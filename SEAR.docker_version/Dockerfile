# This is file is used to build the sear image
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
r-base

RUN cachebuster=b953b35 git clone git://github.com/wpmr2/SEAR && cp /SEAR/SEAR_bin/vcfutils.pl /usr/bin/ && cp /SEAR/Sear.pm /usr/lib/perl5/

ENV PATH /SEAR:$PATH

