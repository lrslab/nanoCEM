from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="nanoCEM",
    version="0.0.2.1",
    author="GUO Zhihao",
    author_email="qhuozhihao@icloud.com",
    description='A simple tool designed to visualize the features that distinguish between two groups of ONT data at the site level.\
                It supports two re-squiggle program(tombo and f5c).',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lrslab/nanoCEM",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7,<3.10',
    install_requires=[
        'h5py==3.8.0',
        'numpy>=1.23.0',
        'pandas>=1.1.5',
        'plotnine>=0.8.0',
        'tqdm>=4.62.0',
        "pysam>=0.21.0",
        "pyslow5>=1.0.0",
        "vbz_h5py_plugin>=1.0.1",
        "biopython>=1.80"
    ],
    scripts=['nanoCEM/current_events_magnifier.py','nanoCEM/extract_sub_fast5_from_bam.py','nanoCEM/extract_sub_fastq_from_bam.py']
)