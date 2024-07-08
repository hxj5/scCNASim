# setup.py


from codecs import open
from os import path
from setuptools import setup, find_packages


here = path.abspath(path.dirname(__file__))

# load configures
exec(open("./sccnvsim/app.py").read())

# Get the long description from the relevant file
with open(path.join(here, "README.md"), encoding='utf-8') as f:
    long_description = f.read()

# set "numpy<2"?
reqs = ["anndata", "intervaltree", "matplotlib", "numpy", "pandas", "pysam",
        "scipy", "seaborn", "statsmodels"]

# Dependency:
# - samtools: used for sort_bam_by_tag.

setup(
    name = "sccnvsim",

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version = VERSION,

    description = "scCNVSim - Allele-specific CNV Simulator trained on scRNA-seq data",
    long_description = long_description,
    long_description_content_type = "text/markdown",

    # The project's main homepage.
    url = "https://github.com/hxj5/scCNVSim",

    # Author details
    author = 'Xianjie Huang',
    author_email = 'xianjie5@connect.hku.hk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['CNV', 'Simulator', "scRNA-seq"],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages = find_packages(),

    #entry_points={
    #    'console_scripts': [
    #        'sccnvsim = sccnvsim.main:main'
    #    ],
    #},

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires = reqs,

    py_modules = ['sccnvsim']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...
)