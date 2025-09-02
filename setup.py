# setup.py

from setuptools import setup, find_packages

setup(
    name='scrin',
    version='1.0.2',
    packages=find_packages(),
    install_requires=[
        'mpi4py',
        'numpy',
        'scipy',
        'pandas',
        'scikit-learn',
        'pyarrow',
        'msgpack',
        'statsmodels',
        'tqdm',
        'rtree',
        'tools'
    ],
    entry_points={
        'console_scripts': [
            'scrin = scrin.main:main',  # Register the main function of the SCRIN package
        ],
    },
)
