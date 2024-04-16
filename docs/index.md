# Welcome to the nanoCEM's documentation !


 <center>![logo](logo.png "nanoCEM") </center>


<center><span style="color:#0084A9;">nanoCEM (nanopore Current Events Magnifier) </span></center>

A Python command-line tool designed to  visualize and analyse multiple statistical features of current events from nanopore sequencing data.

## Installation
<a href="https://pypi.python.org/pypi/nanoCEM" rel="pypi">![PyPI](https://img.shields.io/pypi/v/nanoCEM?color=green) </a>

Before pip install, make sure you have installed the `samtools`(>=1.16) , `f5c`(>=1.4), `slow5tools`(>=1.1.0), `minimap2`(>=2.17), and python  (>=3.7,<=3.10)

    conda install samtools=1.16 minimap2 f5c=1.4 slow5tools -c conda-forge -c bioconda 

To install the latest nanoCEM,

    pip install nanoCEM

And install from the resource

    git clone https://github.com/lrslab/nanoCEM.git
    cd nanoCEM/
    pip install .
To install nanoCEM from [docker](https://hub.docker.com/r/zhihaguo/nanocem_env),

    docker pull zhihaguo/nanocem_env
    
To check the version of nanoCEM, run:

    pip list | grep nanoCEM

## Solutions for some potential environment problem
Although it does not affect the functionality, the issue of possible missing header files caused by **samtools** installation by conda can be resolved with the following command.

    conda install -c conda-forge ncurses

The potential **vbz** format compression issue when reading **fast5** files.

    conda install -c bioconda ont_vbz_hdf_plugin

## Content

* [Quick start](tutorials.md)
* [Output files description](output_format.md)
* [Data preparation from raw reads](preparation.md)
* [Command line arguments](argument.md)
* [Description of statistical features](statistics.md)
* [Adaption for preprocessing methods](adaption.md)
* [An example for application](example.md)
* [Contacts](contact.md)
