import setuptools
import os

from src.pytornado.__version__ import __version__

# See also: https://github.com/kennethreitz/setup.py/blob/master/setup.py

NAME = 'pytornado'
VERSION = __version__
AUTHOR = 'Aaron Dettmann'
EMAIL = 'dettmann@kth.se'
DESCRIPTION = 'A vortex-lattice implementation (VLM) implementation'
URL = 'https://github.com/airinnova/pytornado'
REQUIRES_PYTHON = '>=3.6.0'
REQUIRED = [
    'airfoils>=0.2.0',
    'aeroframe>=0.0.4',
    'ambiance',
    'numpy',
    'scipy',
    'matplotlib>=3.0.2',
    'commonlibs>=0.5.0',
    'model-framework',
]
README = 'README.rst'
PACKAGE_DIR = 'src/'
LICENSE = 'Apache License 2.0'
SCRIPTS = []


here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, README), "r") as fp:
    long_description = fp.read()

setuptools.setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    long_description=long_description,
    url=URL,
    include_package_data=True,
    scripts=SCRIPTS,
    package_dir={'': PACKAGE_DIR},
    license=LICENSE,
    # packages=[NAME],
    packages=setuptools.find_packages(where=PACKAGE_DIR),
    python_requires=REQUIRES_PYTHON,
    install_requires=REQUIRED,
    # See: https://pypi.org/classifiers/
    classifiers=[
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.6',
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
