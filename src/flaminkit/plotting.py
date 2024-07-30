"""Custom plotting tools for FLAMINGO

Under development!
"""

# some generic ones to start
axlabels = {
    "StellarMass": "m_{\u2605}",
    "GasTemperature": "T",
}
# these can probably be generalized to also add special cases for Apertures and In/ExclusiveSphere
for over in (
    "2500_crit",
    "500_crit",
    "200_crit",
    "50_crit",
    "5xR_500_crit",
    "200_mean",
):
    *od, ref = over.split("_")
    # for 5xR_500
    if not isinstance(od, str):
        od = "".join(od).replace("x", "")
    od = f"{od}\\mathrm{ref[0]}"
    # masses
    axlabels[f"SO/{over}/TotalMass"] = f"M_{{{od}}}"
    axlabels[f"SO/{over}/DarkMatterMass"] = f"M_{{\mathrm{{DM}},{od}}}"
    axlabels[f"SO/{over}/GasMass"] = f"M_{{\mathrm{{gas}},{od}}}"
    axlabels[f"SO/{over}/HotGasMass"] = f"M_{{\mathrm{{hot\,gas}},{od}}}"
    axlabels[f"SO/{over}/StellarMass"] = f"M_{{\u2605,{od}}}"
    # others (add as needed)
    axlabels[f"SO/{over}/ComptonY"] = f"Y_{{{od}}}"
    axlabels[f"SO/{over}/GasTemperature"] = f"T_{{{od}}}"


def axlabel_exists(column):
    return column in axlabels


def get_axlabel(column):
    """Produce an axis label based on a column name.

    Axis labels are looked up in a hard-coded dictionary. If the
    specific column (e.g., ``SO/200_crit/TotalMass``) exists in this
    dictionary then this specific label is returned. If it does not,
    this function will look for a generic label for ``TotalMass``. If
    that also does not exist it will return simply ``"TotalMass"``

    Columns should be added to this dictionary as needed

    NOTE: Units not implemented
    """
    if column in axlabels:
        label = axlabels.get(column)
    else:
        col = column.split("/")[-1]
        if col in axlabels:
            label = axlabels.get(col)
        else:
            return column
    return f"${label}$"
