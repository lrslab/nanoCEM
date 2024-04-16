from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="nanoCEM",
    version="0.0.5.7",
    author="GUO Zhihao",
    author_email="qhuozhihao@icloud.com",
    description='A simple tool designed to visualize the features that distinguish between two groups of ONT data at the site level.\
                It supports 4 re-squiggle program(tombo resquiggle/f5c resquiggle/f5c eventalign/move_table).',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lrslab/nanoCEM",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7.0,<=3.11.7',
    install_requires=[
        'h5py>=3.8.0',
        'numpy>=1.23.0',
        'pandas>=1.5.0',
        'plotnine==0.12.4',
        'tqdm>=4.62.0',
        "pysam>=0.21.0",
        "pyslow5>=1.0.0",
        "vbz_h5py_plugin>=1.0.1",
        "biopython>=1.80",
        "scikit-learn>=1.2.2",
        'squigualiser==0.5.1'
    ],
    scripts=['nanoCEM/current_events_magnifier','nanoCEM/extract_sub_fast5_from_bam','nanoCEM/alignment_magnifier']
)