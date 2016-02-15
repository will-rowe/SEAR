![alt tag](http://i.imgsafe.org/22ed408.jpg)

-----------


SEAR: Search Engine for Antimicrobial Resistance
====

Contents
----
[General updates](https://github.com/wpmr2/sear#general-updates)

[Introduction](https://github.com/wpmr2/sear#introduction)

[Versions](https://github.com/wpmr2/sear#versions)

[References](https://github.com/wpmr2/sear#references)

[Contact](https://github.com/wpmr2/sear#contact)

-----------


General updates
====

    17.12.2015:   The SEAR App has been published on BaseSpace.
    01.12.2015:   A Docker container image for SEAR has been created.

-----------


Introduction
====

The aim of this project is to provide an easy to use tool that can construct full-length, horizontally acquired **Antimicrobial Resistance Genes (ARGs)** from sequencing data.

SEAR has been designed with **environmental metagenomics** / **microbiome** experiments in mind, where the diversity and relative abundance of ARGs would like to be determined both quickly and easily.

SEAR outputs: ARG annotation, relative ARG abundance (%), consensus sequence and search results (from ARDB, RAC, NCBI nucleotide and protein databases).

A paper describing the design and implementation of SEAR has been published in <a href="http://doi.org/10.1371/journal.pone.0133492" target="_blank">PloS ONE</a>.

-----------


Versions
====

There are several versions of SEAR available:

1.	SEAR App
----

+ The SEAR App is available on the <a href="https://basespace.illumina.com/home/sequence" target="_blank">BaseSpace</a> platform (Illumina).

+ This is the recommended version to use as it is actively being developed and has the latest version of the SEAR report output.

+ Visit the SEAR App <a href="https://basespace.illumina.com/apps/2083081/SEAR-Antibiotic-Resistance?preferredversion" target="_blank">page</a> to find out more and try the App!


2.	SEAR.docker_version
----

+ The Docker version of SEAR is the version behind the SEAR App. This version of SEAR runs using a 'Docker container' and will have all the required software and dependencies to run properly. An introduction to using Docker can be found <a href="https://training.docker.com/" target="_blank">here</a>.

+ A user can use the SEAR Docker container on a local machine and run SEAR as they would a command line script. 

+ Alternatively, the SEAR Docker container can be run on a server (as in the App version, which uses Amazon cloud servers).

+ The SEAR Docker container image can be pulled directly from the <a href="https://hub.docker.com/r/wpmr/sear/" target="_blank">Docker Hub</a>:

    `docker pull wpmr/sear`

+ Alternatively, the SEAR Docker container image can be built using the Dockerfile:

    `docker build -t 'docker:sear' .`


3.	SEAR.original_cl_version
----

+ The original version of SEAR is still available. This is the command line version which requires the user to set up all the required dependencies. It also requires USEARCH still, unlike the other SEAR versions that now use VSEARCH (an open source alternative to USEARCH).

+ The README.md file in SEAR.original_cl_version explains how to set up this version of SEAR.


4.	SEAR.original_web_version
----

+ The original of the SEAR and web interface is still available. This is an interface to run SEAR on a personal web server. 

+ Uncompress the SEAR_web_version.tar.gz file and view the README.md file for instructions.

+ Screenshots are available <a href="http://i.imgsafe.org/db9ea41.png" target="_blank">here</a>, <a href="http://i.imgsafe.org/26b32f1.png" target="_blank">here</a> and <a href="http://i.imgsafe.org/3b422ba.png" target="_blank">here</a>.

-----------


References
====

SEAR calls on several pieces of bioinformatic software: 


Number | Reference | Function
------ | --------- | --------
1 | [VSEARCH](https://github.com/torognes/vsearch) | `used to quality check, convert and cluster reads`
2 | [BWA-MEM](http://arxiv.org/abs/1303.3997) | `used to map clustered reads to ARG reference sequences`
3 | [SAMtools](http://www.ncbi.nlm.nih.gov/pubmed/19505943) | `used to create ARG consensus sequences from mapped reads`
4 | [BLAST](http://www.ncbi.nlm.nih.gov/pubmed/2231712) | `used to compare ARG consensus sequences to online databases`
5 | [USEARCH](http://www.drive5.com/usearch/) | `used in place of VSEARCH in older SEAR versions`

SEAR requires a database for the initial clustering of reads. This can be a user supplied database (multifasta file) but a copy of the ARG-ANNOT database is included with SEAR and is used by default:

+ [ARG-ANNOT](http://www.ncbi.nlm.nih.gov/pubmed/24145532)

-----------


Contact
====

Please visit the <a href="https://groups.google.com/forum/#!forum/sear" target="_blank">SEAR forum</a> to post a question, report a bug or start a discussion.

Alternatively, please send an email to Will Rowe at: [will.rowe@nhs.net](mailto:will.rowe@nhs.net?subject=SEAR).
