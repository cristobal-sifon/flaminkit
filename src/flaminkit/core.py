from argparse import ArgumentParser
import h5py
from multiprocessing.pool import ThreadPool
import numpy as np
import os
import pandas as pd
import swiftsimio as sw
from tqdm import tqdm
import unyt
import warnings

# debugging
from icecream import ic
from time import time

warnings.simplefilter("once", category=RuntimeWarning)

# see https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.merge.html
pd.options.mode.copy_on_write = True


def halofile(args, info=True):
    """Name of file containing halo information, based on command-line arguments

    Parameters
    ----------
    args : output of ``parse_args``

    Returns
    -------
    halofile: ``str``
        filename
    """
    snap = f"{args.snapshot:04d}"
    halofile = os.path.join(args.path["SOAP"], f"halo_properties_{snap}.hdf5")
    if not info:
        return halofile
    # to have some info handy
    ic(halofile)
    with h5py.File(halofile) as f:
        ic(f.keys())
        for key in f.keys():
            ic(key)
            for gr in f[key].items():
                ic(gr)
                if gr[0] in ("HBTplus", "200_crit", "100kpc"):
                    for i, subgr in enumerate(gr[1].items()):
                        ic(subgr)
                        if subgr[0] in ("projx",):
                            for i, subsubgr in enumerate(subgr[1].items()):
                                ic(subsubgr)
            print()
    return halofile


def infalling_groups(
    halofile,
    cluster_mass_min=0,
    cluster_mass_max=np.inf,
    clusters=None,
    distance_max=5,
    group_mass_min=0,
    n=None,
    overdensity="200_crit",
    so_cols=None,
    random_seed=None,
):
    """Find groups falling into a sample of clusters

    This function looks for all clusters in the chosen mass range that are the
    most massive cluster within ``distance_max`` and then identifies all groups
    around them. Here, groups are defined simply by a lower mass cut. The cluster
    center is taken as the center of mass.

    Parameters
    ----------
    halofile : ``str``
        hdf5 file name
    cluster_mass_min, cluster_mass_max : ``float``, optional
        minimum and maximum spherical overdensity cluster mass, in solar masses
    distance_max : ``float``
        maximum 3d distance around which to search for groups, in units of
        the ``SORadius`` specified by ``overdensity``
    group_mass_min : ``float``, optional
        minimum infall group mass
    n : ``int``, optional
        number of clusters to return, chosen randomly given the mass range.
        If not specified all clusters are returned
    overdensity : ``str``
        spherical overdensity as named in ``halofile``
    so_cols : ``list``, optional
        list of spherical overdensity columns to include in addition to
        ``TotalMass`` -- NOT YET IMPLEMENTED
    random_seed : ``int``, optional
        random seed used to select ``n`` clusters, if applicable.

    Returns
    -------
    main_clusters: ``pd.DataFrame``
        most massive clusters within their own ``distance_max * SORadius``.
        If ``clusters`` is provided this will be a subset of it
    infallers : ``pd.DataFrame``
        infalling groups
    """
    with h5py.File(halofile) as file:
        hostid = file.get("InputHalos/HBTplus/HostFOFId")[()]
        mass = file.get(f"SO/{overdensity}/TotalMass")[()]
        # working with central galaxies is enough here
        good = (hostid > -1) & (file.get("InputHalos/IsCentral")[()] == 1)
        df = pd.DataFrame(
            {
                "TrackId": file.get("InputHalos/HBTplus/TrackId")[()][good],
                "HostFOFId": hostid[good],
                "TotalMass": mass[good],
                "SORadius": file.get(f"SO/{overdensity}/SORadius")[()][good],
            }
        )
        com = file.get(f"SO/{overdensity}/CentreOfMass")[()][good]
        for i, x in enumerate("xyz"):
            df[f"CentreOfMass_{x}"] = com[:, i]
        if clusters is None:
            clusters = (mass[good] > cluster_mass_min) & (mass[good] < cluster_mass_max)
            clusters = df.loc[clusters]
        else:
            assert isinstance(clusters, pd.DataFrame)
            assert "TrackId" in clusters.columns
            cols = ["HostFOFId", "TotalMass", "SORadius"] + [
                f"CentreOfMass_{x}" for x in "xyz"
            ]
            addcols = [col for col in cols if col not in clusters.columns]
            if len(addcols) > 0:
                # needed to merge
                addcols = ["TrackId"] + addcols
                clusters = clusters.merge(df[addcols], on="TrackId", how="left")
    # subsample?
    if n is not None:
        rdm = np.random.default_rng(random_seed)
        n = rdm.choice(
            clusters["TrackId"].size,
            n,
            replace=False,
            shuffle=False,
        )
        clusters = clusters.iloc[n]
    groups = df.loc[df["TotalMass"] > group_mass_min]
    del df
    # main_clusters = []
    infallers = {
        "TrackId": [],
        "HostFOFId": [],
        "CentreOfMass_x": [],
        "CentreOfMass_y": [],
        "CentreOfMass_z": [],
        "TotalMass": [],
        "SORadius": [],
        "MainTrackId": [],
        "DistanceToMain": [],
    }
    for cl in clusters.itertuples():
        dist = (
            (cl.CentreOfMass_x - groups.CentreOfMass_x) ** 2
            + (cl.CentreOfMass_y - groups.CentreOfMass_y) ** 2
            + (cl.CentreOfMass_z - groups.CentreOfMass_z) ** 2
        ) ** 0.5
        # we need the main clusters in here too, to do the merging below
        near = dist < distance_max * cl.SORadius
        if not np.any(groups["TotalMass"].loc[near] > cl.TotalMass):
            # main_clusters.append(cl.Index)
            infallers["MainTrackId"].extend(near.sum() * [cl.TrackId])
            for col in infallers.keys():
                if col in ("MainTrackId", "DistanceToMain"):
                    continue
                infallers[col].extend(groups[col].loc[near])
            infallers["DistanceToMain"].extend(dist[near])
    infallers = pd.DataFrame(infallers)
    main_clusters = clusters.loc[np.isin(clusters["TrackId"], infallers["MainTrackId"])]
    # now keep only infallers
    infallers = infallers.loc[infallers["DistanceToMain"] > 0]
    infallers.sort_values(["MainTrackId", "TotalMass"])
    return main_clusters, infallers


def particles_around(
    particle_file,
    coords,
    dmax,
    particle_type,
    dmin=0 * unyt.Mpc,
    squeeze=True,
    progress=True,
    nthreads=1,
):
    """Find particles no further than a given distance from sets of coordinates

    Parameters
    ----------
    particle_file : ``str``
    coords : ``unyt.unyt_array``, `shape: (N, 3)`
        coordinates of interest
    dmax : ``unyt.unyt_quantity``
        maximum 3d distance to include particles
    particle_type : `iterable`
        items must be any subset of ``("dm", "gas", "stars")``
    dmin : ``unyt.unyt_quantity``, optional
        minimum 3d distance
    squeeze : ``bool``
        return squeezed array if only one particle type is provided
    progress : ``bool``
        whether to show a progress bar. Disabled if ``nthreads>1``
    nthreads : ``int``
        number of threads to use to match particles.


    Returns
    -------
    particles : ``swiftsimio.reader.SWIFTDataset`` or list
        Loaded particles within cubical masks
    matching : ``list`` of ``np.ndarray``, each with ``len = len(coords)``
        List of indices pointing to ``particles`` corresponding to all
        particles around each subhalo within the distance constraints
    """
    assert isinstance(dmax, unyt.unyt_quantity)
    assert isinstance(dmin, unyt.unyt_quantity)
    assert isinstance(squeeze, bool)
    assert isinstance(progress, bool)
    # assert particle types
    _valid_types = ["dm", "gas", "stars"]
    if isinstance(particle_type, str):
        particle_type = [particle_type]
    for pt in particle_type:
        if pt not in _valid_types:
            raise ValueError(f"particle type {pt} not recognized")
    if nthreads > 1:
        progress = False
    _tqdm = tqdm if progress else lambda x, **kwargs: x
    if coords.units != dmax.units:
        coords = coords.to(dmax.units)
    mask = sw.mask(particle_file)
    regions = (
        np.transpose(
            [
                [coords[:, i] - dmax, coords[:, i] + dmax]
                for i in range(coords.shape[1])
            ],
            axes=(2, 0, 1),
        )
        * dmax.units
    )
    mask.constrain_spatial(regions[0])
    for region in regions[1:]:
        mask.constrain_spatial(region, intersect=True)
    particles = sw.load(particle_file, mask=mask)
    p = [[]] * len(particle_type)
    matching = [[]] * len(p)
    # convenience
    dmid = dmax - dmin
    for i, pt in enumerate(particle_type):
        if pt == "dm":
            p[i] = particles.dark_matter
        elif pt == "gas":
            p[i] = particles.gas
        elif pt == "stars":
            p[i] = particles.stars
        ic(p[i])
        # using indices instead of boolean masks so I don't have to keep
        # the entire arrays in memory
        rng = np.arange(p[i].masses.size, dtype=int)
        ti = time()
        if nthreads == 1:
            matching[i] = [
                _find_particles_around(p[i].coordinates, xyz, dmid, rng)
                for xyz in _tqdm(coords, total=len(coords))
            ]
        else:
            matching[i] = [[]] * len(coords)
            pool = ThreadPool(nthreads)
            out = [
                pool.apply_async(
                    _find_particles_around,
                    args=(p[i].coordinates, xyz, dmid, rng),
                    kwds=dict(idx=j),
                )
                for j, xyz in enumerate(coords)
            ]
            pool.close()
            pool.join()
            for match in out:
                j, match = match.get()
                matching[i][j] = match
        tf = time()
        ic("matching", nthreads, tf - ti)
    if len(particle_type) == 1 and squeeze:
        p = p[0]
        matching = matching[0]
    return p, matching


def _find_particles_around(particle_coords, xyz, dmid, rng, idx=None):
    """Helper function for ``particles_around``"""
    matches = rng[
        np.abs(((particle_coords - xyz) ** 2).sum(axis=1) ** 0.5 - dmid) < dmid / 2
    ]
    # return indices when multiprocessing
    if idx is None:
        return matches
    return idx, matches


def subhalos_in_clusters(
    halofile,
    cluster_mask=None,
    clusters=None,
    n=None,
    subhalo_mask=None,
    overdensity="200_crit",
    subhalo_cols=None,
    so_cols=None,
    random_seed=None,
):
    """Find subhalos within a given cluster population

    .. note::

        Only works with HBT+ catalogs for now

    Parameters
    ----------
    halofile : ``str``
        hdf5 file name
    cluster_mass_min, cluster_mass_max : ``float``, optional
        minimum and maximum spherical overdensity cluster mass, in Msun
    cluster_mask : ``dict``, optional
        minimum and maximum values for cluster propoerties. Each entry in
        the dictionary should correspond to a ``Dataset`` name in ``halofile``,
        and its value correspond to an iterable with elements ``(vmin,vmax)``.
        The mask is applied as ``vmin <= x < vmax``. Ignored if ``clusters``
        is provided. For instance, ::

            cluster_mask = {"SO/200_crit/TotalMass": (1e15, np.inf)}

    clusters : ``pd.DataFrame``, optional
        cluster sample. Must contain at least the columns
        ``(HostFOFId,CentreOfMass_x,CentreOfMass_y,CentreOfMass_z)``
    n : ``int``, optional
        number of clusters to return, chosen randomly given the mass range.
        If not specified all clusters that comply with ``cluster_mask`` are returned
    subhalo_mask : ``dict``, optional
        minimum and maximum subhalo values, given as an iterable of ``(vmin, vmax)`` for each
        desired column. Dictionary keys can be any subset of columns
        from ``BoundSubhalo``. The mask is applied as ``vmin <= x < vmax``
    overdensity : ``str``
        spherical overdensity as named in ``halofile``
    subhalo_cols : ``list``, optional
        list of columns from ``BoundSubhalo`` to include in
        addition to ``(HostFOFId,TrackId,Depth)``. Because the object returned
        is a ``pd.DataFrame``, columns with more than one dimension are not yet
        supported
    so_cols : ``list``, optional
        list of spherical overdensity columns to include in addition to
        ``TotalMass``. Ignored if ``clusters`` is provided

    Returns
    -------
    cluster_galaxies : ``pd.DataFrame``
        galaxies within clusters
    """
    assert cluster_mask is None or isinstance(
        cluster_mask, dict
    ), f"cluster_mask must be a dictionary, received {type(cluster_mask)}"
    assert (n is None) or (
        isinstance(n, int) & (n > 0)
    ), f"if provided, n must be a positive integer; received {n}"
    if isinstance(cluster_mask, dict):
        if clusters is not None:
            warnings.warn("cluster_mask ignored when clusters is provided")
            cluster_mask = None
    with h5py.File(halofile) as file:
        hostid = file.get("InputHalos/HBTplus/HostFOFId")[()]
        valid = hostid > -1
        hostid = hostid[valid]
        # merge centrals and satellites
        com = file.get("BoundSubhalo/CentreOfMass")[()][valid]
        depth = file.get("InputHalos/HBTplus/Depth")[()][valid]
        galaxies = pd.DataFrame(
            {
                "TrackId": file.get("InputHalos/HBTplus/TrackId")[()][valid],
                "HostFOFId": hostid,
                "Depth": depth,
            }
        )
        for i, coord in enumerate("xyz"):
            galaxies[coord] = com[:, i]
        if subhalo_cols is not None:
            for col in subhalo_cols:
                galaxies[col] = file.get(f"BoundSubhalo/{col}")[()][valid]
        # in case there are constraints on subhalos
        mask = True
        if subhalo_mask is not None:
            for col, (vmin, vmax) in subhalo_mask.items():
                xmask = file.get(f"BoundSubhalo/{col}")[()][valid]
                mask = mask & (xmask >= vmin) & (xmask < vmax)
                ic(col, xmask.min(), xmask.max(), vmin, vmax, mask.sum())
            del xmask
        galaxies = galaxies.loc[mask]
        # load clusters
        if clusters is None:
            cols = []
            clmask = depth == 0
            if isinstance(cluster_mask, dict):
                for col, (vmin, vmax) in cluster_mask.items():
                    v = file.get(col)[()][valid]
                    clmask = clmask & (v >= vmin) & (v < vmax)
                    cols.append((col, v))
                cols = {col: v[clmask] for col, v in cols}
            clusters = {"HostFOFId": hostid[clmask]}
            if isinstance(cluster_mask, dict):
                clusters = {**clusters, **cols}
            clusters = pd.DataFrame(clusters)
            if so_cols is not None:
                if isinstance(so_cols, str):
                    so_cols = (so_cols,)
                for col in so_cols:
                    clusters[f"SO/{col}"] = file.get(f"SO/{overdensity}/{col}")[()][
                        valid
                    ][clmask]
        else:
            # let's just make sure it contains everything we need
            assert isinstance(clusters, pd.DataFrame)
            musthave = ["HostFOFId", "TrackId"]
            assert np.isin(clusters.columns, musthave).sum() == len(
                musthave
            ), f"clusters must contain columns {musthave}"
            bcg = np.isin(galaxies["TrackId"], clusters["TrackId"])
    ic(clusters, clusters["HostFOFId"], clusters["HostFOFId"].size)
    if n is not None:
        rdm = np.random.default_rng(random_seed)
        n = rdm.choice(
            clusters["HostFOFId"].size,
            n,
            replace=False,
            shuffle=False,
        )
        clusters = clusters.iloc[n]
    # we don't need this as it's in the galaxies
    if "TrackId" in clusters.columns:
        clusters.pop("TrackId")
    cluster_galaxies = clusters.merge(
        galaxies, how="inner", on="HostFOFId", suffixes=("_cl", "_gal")
    )
    if "Depth" in cluster_galaxies.columns:
        cluster_galaxies = cluster_galaxies.sort_values("Depth", ignore_index=True)
    return cluster_galaxies.sort_values("HostFOFId", ignore_index=True)


def subhalo_particle_statistic(
    particle_data, subhalo_particles, statistic, weights=None, **kwargs
):
    """Calculate subhalo properties from the particle distribution

    Parameters
    ----------
    particle_data : ``swiftsimio.objects.cosmo_array``
        particle property, i.e., column from a ``swiftsimio.reader.SWIFTDataset``
    subhalo_particles : ``np.ndarray`` or ``list`, `len: N`
        list of indices or masks relating elements of ``particle_data`` to each
        subhalo
    statistic : ``callable``
        function to apply to the particles of each subhalo.
    weights : ``swiftsimio.objects.cosmo_array`` or ``np.ndarray``, optional
        weights applied to obtain a weighted statistic
    kwargs : ``dict``, optional
        additional arguments passed to ``statistic``. Note that if this includes
        another ``swiftsimio.objects.cosmo_array`` (e.g., as weights) it is
        recommended to give the ``ndarray_view`` to ensure proper behaviour

    Returns
    -------
    subhalo_property : ``np.ndarray``, `shape: (N,)`
        subhalo property obtained from the particles

    Notes
    -----
    - If a property has more than one dimension (e.g., a position or velocity)
      it is easiest to pass each array to this function separately as ::

            mean_velocity = [subhalo_particle_statistic(
                particles.velocities[:,i], subhalo_particles, np.mean)
                for i in range(3)]

    """
    # this is the fastest way to get ndarray behaviour
    data = particle_data.ndarray_view()
    if weights is None:
        return np.array([statistic(data[p], **kwargs) for p in subhalo_particles])
    if isinstance(weights, sw.objects.cosmo_array):
        weights = weights.ndarray_view()
    # this avoids ZeroDivisionError if all the weights are zero for a
    # particular subhalo (e.g. if weighting by gas mass)
    f = lambda x, w, **kwargs: statistic(x, weights=w, **kwargs) if np.any(w > 0) else 0
    return np.array([f(data[p], weights[p], **kwargs) for p in subhalo_particles])


def parse_args(args=None):
    """Parse command-line arguments

    Parameters
    ----------
    args : ``list``-like, optional
        additional arguments to include in the parser. Each element in
        ``args`` should contain two elements: the string(s) enabling the
        argument and the kwargs to add it to the parser. For instance, ::

            args=(('--foo', {'type': int, 'default': 1}), ('--bar', {'action': 'store_true'}))

    Returns
    -------
        args : output of `parser.parse_args`
    """
    parser = read_args()
    if args is not None:
        if isinstance(args[0], str):
            args = (args,)
        for argname, kwargs in args:
            parser.add_argument(argname, **kwargs)
    args = parser.parse_args()
    # for now
    args.path = dict(main=os.path.join(os.environ.get("FLAMINGO"), args.box, args.sim))
    args.path["particles"] = os.path.join(
        args.path.get("main"), "snapshots_downsampled"
    )
    args.path["SOAP"] = os.path.join(args.path.get("main"), "SOAP-HBT")
    if args.snapshot is not None:
        args.path["snapshot"] = os.path.join(
            args.path.get("main"), "snapshots", f"flamingo_{args.snapshot:04d}"
        )
        args.snapshot_file = os.path.join(
            args.path.get("snapshot"), f"flamingo_{args.snapshot:04d}.hdf5"
        )
    # is there a cleaner way?
    args.min_mass_cluster /= 1e10
    args.min_mass_group /= 1e10
    args.min_mass_stars /= 1e10
    if args.test:
        args.debug = True
    if not args.debug:
        ic.disable()
    return args


def read_args():
    """Set up the base command-line arguments

    Call this function from any program if there are additional
    command-line arguments in it (to be added manually from that
    program)

    Note that any program might ignore any of these default
    command-line arguments if only a subset are relevant

    Returns
    -------
    parser : ``argparse.ArgumentParser`` object
    """
    parser = ArgumentParser()
    add = parser.add_argument
    add("-b", "--box", default="L1000N1800")
    add("-s", "--sim", default="HYDRO_FIDUCIAL")
    add("-z", "--snapshot", default=77, type=int)
    add("--debug", action="store_true")
    add(
        "--nthreads",
        default=1,
        type=int,
        help="Number of threads (used for specific tasks)",
    )
    add("--seed", default=1, type=int, help="Random seed")
    add("--test", action="store_true")
    add("--min-mass-cluster", default=1e15, type=float, help="Minimum cluster mass")
    add("--min-mass-group", default=1e13, type=float, help="Minimum mass for groups")
    add("--min-mass-stars", default=1e10, type=float, help="Minimum stellar mass")
    add(
        "--ncl", default=None, type=int, help="Number of main clusters (for quick runs)"
    )
    add("--ngal", default=None, type=int, help="Number of galaxies (for quick runs)")
    return parser
