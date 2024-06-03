# nanoCEM ![logo](docs/logo_tiny.png "nanoCEM")
<a href="https://pypi.python.org/pypi/nanoCEM" rel="pypi">![PyPI](https://img.shields.io/pypi/v/nanoCEM?color=green) </a>
<a href="https://opensource.org/license/mit/" rel="license">![License](https://img.shields.io/pypi/l/nanoCEM?color=orange)</a>

The nanopore current events magnifier (`nanoCEM`) is a python command line to facilitate the analysis of DNA/RNA modification sites by visualizing statistical features of current events. 
NanoCEM can be used to showcase high confidence sites and observe the difference based on the modification sample and the low or no modification sample.

It supports two re-squiggle pipeline(`Tombo` and `f5c`) and support `R9` and `R10`.
If you want to view single read signal or raw signal, [Squigualiser](https://github.com/hiruna72/squigualiser) is recommended.

## Installation

Before pip install, make sure you have installed the `samtools`(>=1.16) , `f5c`(>=1.4), `slow5tools`(>=1.1.0) and `minimap2`(>=2.17),

    conda install samtools=1.16 minimap2 f5c=1.4 slow5tools -c conda-forge -c bioconda 

To install the latest nanoCEM

    pip install nanoCEM

And install from the resource

    git clone https://github.com/lrslab/nanoCEM.git
    cd nanoCEM/
    pip install .
To install nanoCEM from docker,

    docker pull zhihaguo/nanocem_env
    
To check the version of nanoCEM, run:

    pip list | grep nanoCEM


## Solutions for some potential environment problem
Although it does not affect the functionality, the issue of possible missing header files caused by **samtools** installation by conda can be resolved with the following command.

    conda install -c conda-forge ncurses

The potential **vbz** format compression issue when reading **fast5** files.

    conda install -c bioconda ont_vbz_hdf_plugin



## Data release
For the data we used and related commands in our paper, please view our [wiki](https://github.com/lrslab/nanoCEM/wiki/Data-release-and-commands)

## Documents
Full documentation is available at https://nanocem.readthedocs.io/

## Publication
Our paper is online [here](https://academic.oup.com/nargab/article/6/2/lqae052/7676831) now

## Workflow
![workflow](docs/Workflow.png)
