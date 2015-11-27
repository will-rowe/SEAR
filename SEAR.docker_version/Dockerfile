# This Dockerfile is used to build the SEAR docker image for the BaseSpace App.
# Due to licensing issues, USEARCH must be downloaded by the user as it cannot be distributed with SEAR. The user must name the USEARCH executable ‘usearch’ and add it to the docker image (in /usr/bin).

FROM ubuntu:14.04
MAINTAINER Will Rowe <wpmr2@cam.ac.uk>

RUN apt-get update && apt-get install -y \
git \
curl \
make \
samtools \
ncbi-blast+ \
bwa \
tophat \
seqtk \
r-base \
libJSON-perl

RUN cachebuster=b953b35 git clone git://github.com/wpmr2/SEAR /tmp/SEAR && mv /tmp/SEAR/SEAR.docker_version/SEAR / && cp /tmp/SEAR/SEAR.docker_version/BaseSpace_scripts/perl_parser.pl /usr/bin/ && cp /SEAR/SEAR_bin/vcfutils.pl /usr/bin/ && cp /SEAR/Sear.pm /usr/lib/perl5/

ENV PATH /SEAR:$PATH

