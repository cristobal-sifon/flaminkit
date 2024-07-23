from setuptools import setup
import os
import re
from setuptools import find_packages, setup

# folder where pygmos is stored
here = os.path.abspath(os.path.dirname(__file__))


# this function copied from pip's setup.py
# https://github.com/pypa/pip/blob/1.5.6/setup.py
# so that the version is only set in the __init__.py and then read here
# to be consistent
def find_version(fname):
    version_file = read(fname)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


# Taken from the Python docs:
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a
# top level README file and 2) it's easier to type in the README file
# than to put a raw string in below
def read(fname):
    return open(os.path.join(here, fname)).read()


def read_requirements(reqfile):
    return [i for i in open(reqfile).read().split("\n") if i]


setup(
    name="flaminkit",
    version=find_version("src/flaminkit/__init__.py"),
    description="A toolkit for the FLAMINGO simulation suite",
    long_description=read("README.rst"),
    long_description_content_type="text/x-rst",
    author="Cristobal Sifon",
    author_email="cristobal.sifon@pucv.cl",
    url="https://github.com/cristobal-sifon/flaminkit/",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=read_requirements("requirements.txt"),
    python_requires=">=3.7,<3.12",
    zip_safe=False,
)
