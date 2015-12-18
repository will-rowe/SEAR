################################################################################
########## SEAR_bin SETUP INSTRUCTIONS ##########
################################################################################
1. The SEAR_bin should already contain the following executables:
        cat
        gunzip
        mkdir
        mv
        rm
        split
        tar
        touch
        uname

1. Install the following packages from the 1_src directory (SEAR has not been tested with other versions of these programs):
        bwa-0.7.8.tar.bz2
        curl-7.40.0.tar
        ncbi-blast-2.2.28+-x64-linux.tar.gz
        R-2.12.0.tar.gz
        samtools-0.1.19.tar.bz2
        seqtk-master.tar.gz
        tophat-2.0.9.Linux_x86_64.tar.gz

2. Download and install the USEARCH package (v7.0.959 or later) and R (v2.12.0 or later)

3. Make sure that the following executables have been added to the SEAR_bin (or SEAR path) and can be called from the SEAR_bin:
        bam2fastx* (from tophat-2.0.9)
        bcftools* (from samtools-0.1.19/bcftools/)
        blastn*
        blastx*
        bwa*
        curl*  
        makeblastdb*
        Rscript* (from R-2.12.0/bin/)
        samtools*
        seqtk*
        vcfutils.pl* (from samtools-0.1.19/bcftools/)
        usearch

