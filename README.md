################################################################################
########## SEARCH ENGINE FOR ANTIMICROBIAL RESISTANCE ##########
################################################################################
This is the commandline version of SEAR.
SEAR annotates antimicrobial resistance genes (ARGs) from short read, metagenomic datasets.
SEAR calls on several pieces of software.
The results list the annotated ARGs, the relative abundance of each ARG, the consensus sequence for the annotated ARG and links to BLAST results and further information on several online databases.
It is recommended to leave the SEAR directories and files in situ within the ./SEAR_commandline_version directory.
SEAR uses the remote NCBI BLAST service. To speed up SEAR, you can choose to use a local install of NCBI BLAST. Please refer to the commandline_sear.pl to change from the remote BLAST service.

################################################################################
########## SEAR SETUP INSTRUCTIONS ##########
################################################################################
1.  Set up the SEAR_bin by following the 2_README.md file in ./SEAR_bin to install the required software and setup the correct executables in the SEAR_bin. 

2.  Open commandline_SEAR.pl in an editor and input the SEAR SETTINGS for your system (giving complete paths).
        a.  Set the $ENV{'PATH'} variable to point to the '/SEAR_bin'.
        b.  Set the @INC path to the 'SEAR_commandline_version' directory (that contains the Sear.pm module).
        c.  Set the full paths for the remaining variables as instructed.

3.  Check that the following required Perl modules are installed on your system:
                - Getopt::Long
                - Pod::Usage
                - Time::HiRes
                - File::Basename
                - File::Find
                - List::Util
                - LWP::Simple

4.  Give the SEAR_DATA directory and commandline_SEAR.pl program read, write and execute rights (use chmod).


################################################################################
########## SEAR MANUAL ##########
################################################################################
NAME
    SEAR:   Search Engine for Antimicrobial Resistance   (version 0.1)

SYNOPSIS
    perl commandline_sear.pl -i reads1.fastq reads2.fq

ARGUMENTS
    -i/--infile
        paths to input file(s) (files must be .fastq/.fq although can be gzipped e.g. file.fastq.gz)

OPTIONS
    -ff/--fqformat
        ascii offset for the input fastq files. Accepts either 33 or 64 (default=33)

    -lc/--lengthcutoff
        discard sequences with length < lc (default=70)

    -qc/--qualitycutoff
        quality score cutoff for input fastq files (default=20)

    -f/--filter
        filter reads by mapping to Human Genome reference and discarding mapped reads. Accepts either Y or N (default=N)

    -cc/--coveragecutoff
        percentage reference sequence coverage for gene annotation (default=90)

    -ci/--clusteringident
        identity value for usearch clustering (default=0.99)

    -r/--references
        the reference gene dataset to use (default=arg_annot_database.fa)

    -t/--threads
        number of threads to use in steps that allow multi-threading (default=1)

    -h/--help
        prints usage and exits

    -m/--manual
        prints the manual page and exits

DESCRIPTION
    Command line version of the Search Engine for Antimicrobial Resistance (commandline_SEAR.pl).

DEPENDENCIES
    perl                   required modules: Getopt::Long, Pod::Usage, Time::HiRes, File::Basename, File::Find, List::Util, LWP::Simple
    usearch
    bwa
    samtools               incl. bcftools and vcfutils.pl
    ncbiblast
    seqtk
    tophat                 incl. bam2fastx
    r

AUTHOR
    wpmr2@cam.ac.uk

TESTED
    perl       version 5.12.4
    usearch    version 7.0.959
    bwa        version 0.7.8
    samtools   version 0.1.19
    ncbiblast  version 2.2.28+
    tophat     version 2.0.9
    r          version 2.12.9

CAVEATS
    ***The input fastq file must have a common identifier (min. 5 characters) in the read header line in order for hash tables to work correctly.
    ***Amount of RAM used must be able to hold input file/s.
    ***To view the results files correctly, they must be kept within the SEAR_DATA/RESULTS directory as they rely on the associated html files.

FUTURE UPDATES
    Future updates include: Additional blast databases, mobile element identification, utilising paired end read information in mapping steps.
