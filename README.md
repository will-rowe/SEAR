# Search Engine for Antimicrobial Resistance (SEAR)


## Introduction

The aim of this project is to provide an easy to use tool that can construct full-length, horizontally acquired Antimicrobial Resistance Genes (ARGs) from sequencing data.

SEAR has been designed with environmental metagenomics / microbiome experiments in mind, where the diversity and relative abundance of ARGs would like to be determined both quickly and easily.

The paper describing the design and implementation of SEAR has been published in [PLoS ONE](http://doi.org/10.1371/journal.pone.0133492).




## Update

The server that hosts SEAR is being retired. Consequently, a Docker container image has been created and it is the easiest way to run SEAR yourself.

This image can be pulled from [Docker](https://hub.docker.com/r/wpmr/sear/):

`docker pull wpmr/sear`

Alternatively, the container image can be built using the Dockerfile in SEAR.docker_version:

`docker build -t 'docker:sear' .`




## Versions

There are several versions of SEAR available:

### SEAR App on BaseSpace

As of xx.xx.xxxx, SEAR is available as an App on the Illumina BaseSpace analysis platform.

Link here: [SEAR BaseSpace App](www.google.com).


### SEAR.docker_version


### SEAR.command_line_version




## Contact

Please visit the [SEAR Forum](https://groups.google.com/forum/#!forum/sear) to post a question, report a bug or start a discussion.

Alternatively, please send an email to Will Rowe at: [wpmr2@cam.ac.uk](mailto:wpmr2@cam.ac.uk?subject=SEAR).
