# SEAR: Search Engine for Antimicrobial Resistance


## General update (1.12.2015)

The server that hosts SEAR is being retired. Consequently, a Docker container image has been created and it is the easiest way to run SEAR yourself.




## Introduction

The aim of this project is to provide an easy to use tool that can construct full-length, horizontally acquired Antimicrobial Resistance Genes (ARGs) from sequencing data.

SEAR has been designed with environmental metagenomics / microbiome experiments in mind, where the diversity and relative abundance of ARGs would like to be determined both quickly and easily.

The paper describing the design and implementation of SEAR has been published in [PLoS ONE](http://doi.org/10.1371/journal.pone.0133492).




## Versions

There are several versions of SEAR available:

### 1.	SEAR App on BaseSpace

As of xx.xx.xxxx, SEAR is available as an App on the Illumina BaseSpace analysis platform.

Link here: [SEAR BaseSpace App](www.google.com).


### 2.	SEAR.docker_version

The docker version of SEAR is the version behind the SEAR App. This version of SEAR runs using a 'Docker container' and will have all the required software and dependencies to run properly.

A user can use the SEAR container on a local machine and run SEAR as they would a command line script. 

Alternatively, the SEAR container can be run on a server (as in the App version, which uses Amazon cloud servers).

An introduction to Docker can be found [here](https://training.docker.com/).

The SEAR container image can be pulled directly from the [Docker Hub](https://hub.docker.com/r/wpmr/sear/):

`docker pull wpmr/sear`

Alternatively, the SEAR container image can be built using the Dockerfile in SEAR.docker_version:

`docker build -t 'docker:sear' .`



### 3.	SEAR.command_line_version

The original version of SEAR is still available. This is the command line version which requires the user to set up all the required dependencies. It also requires USEARCH still, unlike the other SEAR versions that now use VSEARCH (an open source alternative to USEARCH).

The README.md file in SEAR.command_line_version explains how to set up this version of SEAR.




## References




## Contact

Please visit the [SEAR Forum](https://groups.google.com/forum/#!forum/sear) to post a question, report a bug or start a discussion.

Alternatively, please send an email to Will Rowe at: [wpmr2@cam.ac.uk](mailto:wpmr2@cam.ac.uk?subject=SEAR).
