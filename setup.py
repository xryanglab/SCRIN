# setup.py

from setuptools import setup, find_packages
from pathlib import Path

this_dir = Path(__file__).parent
readme_path = this_dir / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name='scrin',
    version='1.1.0rc2',
    description="SCRIN is a tool for identifying RNA co-localization networks within subcellular spatial transcriptomics data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="xryanglab",
    url="https://github.com/xryanglab/SCRIN",
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
        'rtree'
    ],
    entry_points={
        'console_scripts': [
            'scrin = scrin.main:main',  # Register the main function of the SCRIN package
            'scrin-colocalization-events = scrin.tools.colocalization_events:main',
            'scrin-false-positive-analysis = scrin.tools.false_positive_analysis:main',
            'scrin-extract-distances = scrin.tools.extract_distances:main',
        ],
    },
)
