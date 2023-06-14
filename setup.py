from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="current_events_magnifier",
    version="0.2.3.2",
    author="GUO Zhihao",
    author_email="qhuozhihao@icloud.com",
    description='A sample tool designed to visualize the features that distinguish between two groups of ONT data at the site level.\
                It supports two re-squiggle pipeline(Tombo and f5c).',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lrslab/current_events_magnifier",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
        'h5py==3.8.0',
        'numpy>=1.21.6',
        'pandas>=1.1.5',
        'plotnine>=0.8.0',
        'tqdm>=4.62.0',
        "pyslow5>=1.0.0",
        "pysam>=0.21.0"
    ],
    scripts=['current_events_magnifier/CE_magnifier.py']
)