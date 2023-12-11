# Welcome to the nanoCEM's documentation !


 <center>![logo](logo.png "nanoCEM") </center>


<center><span style="color:#0084A9;">nanoCEM (nanopore Current Events Magnifier) </span></center>

A Python command-line tool designed to  visualize and analyse multiple statistical features of current events from nanopore sequencing data.

## Installation

Before pip installation, make sure you have installed the `samtools`(>=1.16) and `minimap2`(>=2.17),

    conda install -c bioconda samtools minimap2

To install the latest nanoCEM

    pip install nanoCEM

See our Installation page for details. 

To check the version of nanoCEM, run:

    pip list | grep nanoCEM


 **Notes:** Additionally, we do not rely on any re-squiggle or eventalign packages. We only need their index files for the sequencing data.




## Content

* [Quick start](tutorials.md)
* [Output files description](output_format.md)
* [Data preparation from raw reads](preparation.md)
* [Command line arguments](argument.md)
* [An example for application](example.md)
* [Contacts](contact.md)
