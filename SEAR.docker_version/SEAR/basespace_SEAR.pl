#!/usr/bin/perl
# basespace_SEAR.pl is the basespace version of the SEAR website and script that identifies full length antimicrobial resistance genes in fastq sequencing data. It takes the arguments passed from the command line and user specified files and then uses several pieces of open source software.
use warnings;
use strict;
use Sear;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);
use Time::HiRes;
use File::Basename;
use File::Find;
use List::Util 'sum';

############################################################################
#################### SEAR SETTINGS ####################
############################################################################
# The following settings need to be set in order for SEAR to run correctly.
############################################################################
### SET UP ENVIRONMENT ###
############################################################################
# Set ENV variable for the /SEAR_bin (see README.md for how to set up /SEAR_bin correctly).
$ENV{'PATH'} = '/bin:/usr/bin';
# Set path to SEAR_commandline_version directory (containing Sear.pm module)
BEGIN {push @INC, '/SEAR'};

############################################################################
### SET PATHS FOR SERVER ###
############################################################################
# Set webaddress for RAC database;
my $RAC_web_page = "http://rac.aihi.mq.edu.au/rac/feature/list";

# Setup local NCBI blast - not currently used (remote flag is on in NCBI blast commands)
#my $blast_update = "update_blastdb --passive --decompress nr nt";
#print "updating BLAST databases (nr and nt) before starting SEAR\n";
#system ( "$blast_update" ) == 0 or die "Can't download and decrompress BLAST databases: $!\n";

############################################################################
############################################################################
#################### BEGIN SEAR ####################
############################################################################
############################################################################
# The code does not have to be edited from this point on, arguments can be passed from the command line.
############################################################################
### PROCESS COMMAND LINE ARGUMENTS ###
############################################################################
my $start_timestamp = Sear::getDBdate();
# defaults
my $opt_threads = 1;
my $opt_ascii = 33;
my $opt_length = 70;
my $opt_quality = 20;
my $opt_coverage = 90;
my $opt_filter_reads = "N";
my $opt_clustering_identity = "0.99";
my $SEAR_DATA_directory = "./SEAR_DATA";
my $project_output = "0000";
my ( $opt_help, $opt_man, @opt_inputfiles, );
GetOptions(
'help|h'                 =>      \$opt_help,
'manual|m'               =>      \$opt_man,
'threads|t=i'            =>      \$opt_threads,
'infile|i=s{,}'          =>      \@opt_inputfiles,
'filter|f=s'             =>      \$opt_filter_reads,
'fqformat|ff=i'          =>      \$opt_ascii,
'lengthcutoff|lc=i'      =>      \$opt_length,
'qualitycutoff|qc=i'     =>      \$opt_quality,
'coveragecutoff|cc=i'    =>      \$opt_coverage,
'clusteringident|ci=s'	 =>      \$opt_clustering_identity,
'datadirectory|dd=s'	 =>	 \$SEAR_DATA_directory,
'projectoutput|po=s'	 =>	 \$project_output,
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;
if (@opt_inputfiles == 0)
{
    pod2usage("$0: No inputfile specified.\n")
}

#-----------------------------------------------------------------
#----------------  Documentation / Usage / Help ------------------
=head1 NAME
 
=over 8
 
=item B<SEAR:   Search Engine for Antimicrobial Resistance   (version 0.1)>

=back
 
=head1 SYNOPSIS
 
 perl commandline_sear.pl -i reads1.fastq reads2.fq
 
=head1 ARGUMENTS
 
=over 8
 
=item B<-i/--infile>

 paths to input file(s) (files must be .fastq/.fq although can be gzipped e.g. file.fastq.gz)
 
=back

=head1 OPTIONS

=over 8
 
=item B<-ff/--fqformat>
 
 ascii offset for the input fastq files. Accepts either 33 or 64 (default=33)
 
=item B<-lc/--lengthcutoff>
 
 discard sequences with length < lc (default=70)
 
=item B<-qc/--qualitycutoff>
 
 quality score cutoff for input fastq files (default=20)

=item B<-f/--filter>
 
 filter reads by mapping to a reference genome and discarding mapped reads. Y==E.coli genome N=none (default=N)

=item B<-cc/--coveragecutoff>
 
 percentage reference sequence coverage for gene annotation (default=90)
 
=item B<-ci/--clusteringident>
 
 identity value for vsearch clustering (default=0.99 (99%))
 
=item B<-r/--references>
 
 the reference gene dataset to use (default=arg_annot_database.fa)
 
=item B<-t/--threads>
 
 number of threads to use in steps that allow multi-threading (default=1)
 
=item B<-h/--help>
 
 prints usage and exits
 
=item B<-m/--manual>
 
 prints the manual page and exits
 
=back
 
=head1 DESCRIPTION
 
 Command line version of the Search Engine for Antimicrobial Resistance (commandline_SEAR.pl).

=head1 DEPENDENCIES
 
 perl                   required modules: Getopt::Long, Pod::Usage, Time::HiRes, File::Basename, File::Find, List::Util, LWP::Simple
 vsearch
 bwa
 samtools               incl. bcftools and vcfutils.pl
 ncbiblast
 seqtk
 tophat                 incl. bam2fastx
 r                      incl. Rscript
 
=head1 AUTHOR
 
 wpmr2@cam.ac.uk
 
=head1 TESTED
 
 perl       version 5.12.4
 vsearch    version 1.9.3
 bwa        version 0.7.8
 samtools   version 0.1.19
 ncbiblast  version 2.2.28+
 tophat     version 2.0.9
 r          version 2.12.9

=head1 CAVEATS
 
 ***The input fastq file must have a common identifier (min. 5 characters) in the read header line in order for hash tables to work correctly.
 ***Amount of RAM used must be able to hold input file/s.
 ***To view the results files correctly, they must be kept within the SEAR_DATA/RESULTS directory as they rely on the associated html files.
 
=head1 FUTURE UPDATES
 
 Future updates include: Additional blast databases, mobile element identification, utilising paired end read information in mapping steps, more effient fastq lookup
 
=cut


############################################################################
### PROCESS INPUT FILES AND SET UP FILES ###
############################################################################
my $opt_database = "$SEAR_DATA_directory/references/database/arg_annot_database.fa";
my $start_time = [Time::HiRes::gettimeofday()];
print "\n\n################################################\n";
print "Search Engine for Antimicrobial Resistance\n";
print "################################################\nParameters for this SEAR run are as follows:\n";
print "fastq ascii format: $opt_ascii\n";
print "filter against E.coli genome k12: yes\n" if ($opt_filter_reads =~ m/Y/);
print "discard sequences with length < $opt_length bases\n";
print "read quality cutoff: $opt_quality\n";
print "clustering identity: $opt_clustering_identity\n";
print "ARG database for read clustering: $opt_database\n";
print "coverage cuttoff for gene annotation: $opt_coverage%\n";
print "number of threads for multithreading steps: $opt_threads\n";
print "################################################\n";
foreach my $in_file (@opt_inputfiles)
{
    next if ($opt_ascii =~ m/33/);
    next if ($opt_ascii =~ m/64/);
    die "\n\nIncorrect SEAR parameter: fastq encoding must be either 33 or 64\n";
}
foreach my $in_file (@opt_inputfiles)
{
    next if ($opt_filter_reads =~ m/Y/);
    next if ($opt_filter_reads =~ m/N/);
    die "\n\nIncorrect SEAR parameter: filter option must either be Y (e coli K12 genome) or N\n";
}

my $temp_file_name = int ([Time::HiRes::gettimeofday()] * 1000);
my $temp_files_directory = "$SEAR_DATA_directory/$temp_file_name";
my $safe_filename_characters = "a-zA-Z0-9_.-";
print "Supplied input files:\n";
foreach my $in_file (@opt_inputfiles)
{
    print "\t$in_file\n";
    next if ($in_file =~ m/\.fastq$/);
    next if ($in_file =~ m/\.fq$/);
    next if ($in_file =~ m/\.fastq.gz$/);
    next if ($in_file =~ m/\.fq.gz$/);
    die "\n\nIncorrect SEAR parameter: input file is not a fastq file (.fastq/.fq) or gzipped (.fastq.gz/.fq.gz)";
}
print "\n################################################\n";
print "Starting main SEAR program\n";
print "################################################\n";
#   set up files
my $log = "$temp_files_directory/log.txt";

#   check if project id exists in data/output
use File::Path qw(make_path);
my $dir = "/data/output/appresults/$project_output";
eval { make_path($dir) };
if ($@) {
  print "Couldn't create $dir: $@";
}
my $opt_results = $dir . "/SEAR_RESULTS_$temp_file_name";
my $multifastq = "consensus_seqs.fq";
system ( "mkdir -p $temp_files_directory" ) == 0 or die ( "Can't make temporary files directory: $?.\n" );
system ( "touch $log" ) == 0 or die ( "Can't create log file: $?.\n" );
system ( "mkdir $opt_results" ) == 0 or die ("Can't create results directory: $!\n");
system ( "touch $opt_results/$multifastq") == 0 or die ( "Can't create files: $?.\n" );
foreach my $in_file (@opt_inputfiles)
{
    my $gzip_count = 1;
    if ($in_file =~ m/\.gz$/)
    {
	my $uncompressed_file = $temp_files_directory . "SEAR_uncompressed_in_file" . $gzip_count;
        my $gunzip_command = "gunzip -dc $in_file > $uncompressed_file";
        system (" $gunzip_command ") == 0 or die ( "Can't gunzip infile: $gunzip_command $?.\n" );
        $in_file = $uncompressed_file;
	$gzip_count ++;
    }
}


############################################################################
### ADD FASTQ TAGS AND REMOVE WHITE SPACE ###
############################################################################
my $read_tag_counter = 1;
my @parsed_input_filenames;
foreach my $in_file (@opt_inputfiles)
{
    my $input_file_name = fileparse($in_file);
    push(@parsed_input_filenames, $input_file_name);
    open (INFILE, $in_file) or die "Can't open $in_file: $!";
    my $first_read = <INFILE>;
    close INFILE;
    my $read_substr = substr($first_read, 1, 5);
    my $read_tag = "SEARinputFILE$read_tag_counter" . $read_substr;
    my $tagged_in_file = $temp_files_directory . "/SEARinputFILE" . $read_tag_counter . ".fq";
    open (INFILE, $in_file) or die "Can't open $in_file: $!";
    my @in_file_lines = <INFILE>;
    close INFILE;
    foreach my $in_file_line (@in_file_lines)
    {
        if ($in_file_line =~ m/$read_substr/g)
	{
		$in_file_line =~ s/$read_substr/$read_tag/g;
		$in_file_line =~ s/\s//g;
		$in_file_line .= "\n";
    	}
    }
    open OUTFILE, ">$tagged_in_file" or die "Can't open $tagged_in_file for writing tagged infiles: $!";
    print OUTFILE @in_file_lines;
    close OUTFILE;
    $in_file = $tagged_in_file;
    $read_tag_counter ++;
}


############################################################################
### REMOVE CONTAMINATION ###
############################################################################
# Filter reads using Ecoli K12.
if ($opt_filter_reads =~ m/Y/)
{
    #   set up local copy of e coli genome bwa index
    print "\npreparing to filter reads against the E. coli K12 reference genome . . .\n";
    my $setup_hg = "cp /SEAR/SEAR_DATA/references/EcoliK12.tar.gz $temp_files_directory/ && tar -xvf $temp_files_directory/EcoliK12.tar.gz -C $temp_files_directory/";
    system ( "$setup_hg" ) == 0 or die "Can't set up ecoli genome bwa index: $!\n";
    my $hg_path = "$temp_files_directory/EcoliK12/E.coli_K_12.fasta";
    foreach my $in_file (@opt_inputfiles)
    {
        chomp $in_file;
        print "\nfiltering reads in $in_file . . .\n\n";
        my $filtered_reads = $in_file . ".filtered.fq";
        my $bwa_filter_sam_file = "filter_sam_file.sam";
        my $bwa_filter_bam_file = "filter_bam_file.bam";
        my $bwa_filtered_reads_bam_file = "filtered_reads.bam";
        my $bwa_filter_command_1 = "bwa mem -t $opt_threads $hg_path $in_file > $temp_files_directory/$bwa_filter_sam_file";
        my $bwa_filter_command_2 = "samtools view -bS $temp_files_directory/$bwa_filter_sam_file > $temp_files_directory/$bwa_filter_bam_file";
        my $bwa_filter_command_3 = "samtools view -b -f4 $temp_files_directory/$bwa_filter_bam_file > $temp_files_directory/$bwa_filtered_reads_bam_file";
        my $bwa_filter_command_4 = "bam2fastx -q -Q -A -N -o $filtered_reads $temp_files_directory/$bwa_filtered_reads_bam_file";
        my $bwa_filter_cleanup_command = "rm $temp_files_directory/*.sam && rm $temp_files_directory/*.bam";
        system ( $bwa_filter_command_1 ) == 0 or die ( "Can't filter $in_file using bwa mem: $?.\n" );
        system ( $bwa_filter_command_2 ) == 0 or die ( "Can't filter $in_file using bwa mem - can't convert sam > bam: $?.\n" );
        system ( $bwa_filter_command_3 ) == 0 or die ( "Can't filter $in_file using bwa mem: $?.\n" );
        system ( $bwa_filter_command_4 ) == 0 or die ( "Can't convert bam > fastq using bma2fastx (tophat): $?.\n" );
        my $filter_hg_mapped_reads = `samtools view -c -F4 $temp_files_directory/$bwa_filter_bam_file`;
        my $filter_hg_unmapped_reads = `samtools view -c -f4 $temp_files_directory/$bwa_filter_bam_file`;
        print "\nnumber of reads mapped to K12: $filter_hg_mapped_reads\nnumber of reads unmapped to K12: $filter_hg_unmapped_reads\nfiltering complete . . .\n";
        system ( $bwa_filter_cleanup_command ) == 0 or die ( "Can't remove filtering files: $?.\n" );
        $in_file = "$filtered_reads";
        #   check if no reads passed this step
        finddepth(\&filecheck, "$in_file");
    }
}


############################################################################
### LOAD DATA AND CONVERT ###
############################################################################
## Filter the reads and convert to fasta format.
#   for each split fastq file, filter reads based on read length, convert to fasta format and load all reads into a single temp file
print "quality checking and converting to FASTA with vsearch . . .\n";
my $out_file = "reads.fasta";
my $counter = 1;
foreach my $in_file (@opt_inputfiles)
{
    my $fastq_command = "vsearch -fastq_filter $in_file -fastq_minlen $opt_length -fastaout $temp_files_directory/$counter.fasta -notrunclabels -fastq_qmax 100 -threads $opt_threads >> $log 2>&1";
    system (" $fastq_command ") == 0 or die ( "\n\nError in converting fastq to fasta format: check input fastq files\n\n" );
    my $fastareads = "$temp_files_directory/$counter.fasta";
    open (READS_IN, $fastareads) or die ("$!");
    open (OUTFILE, ">>$temp_files_directory/$out_file") or die ("$!");
    print OUTFILE <READS_IN>;
    close (OUTFILE);
    close (READS_IN);
    unlink("$temp_files_directory/$counter.fasta");
    $counter ++;
}

#   check if no reads passed this step
finddepth(\&filecheck, "$temp_files_directory/$out_file");


############################################################################
### CLUSTERING ###
############################################################################
## Each split file of fasta reads is clustered against the reference ARG database.
#   for every split file containing fasta reads in the temp directory, run each as separate vsearch against database of reference sequences
print "clustering with vsearch . . .\n";
my $matched_reads_fa = "matched.reads.fa";
my $vsearch_command = "vsearch --usearch_global $temp_files_directory/$out_file --db $opt_database --uc $temp_files_directory/vsearchfile.uc --matched $temp_files_directory/$matched_reads_fa --id $opt_clustering_identity --strand both --maxhits 1 --threads $opt_threads --notrunclabels --query_cov 0.7 >> $log 2>&1";
system(" $vsearch_command ") == 0 or die ( "Error in vsearch command: $?.\n" );

#   check if no reads passed this step
finddepth(\&filecheck, "$temp_files_directory/$matched_reads_fa");


############################################################################
### PARSING ###
############################################################################
##  Each vsearch output file is searched through and hits are removed and placed in a new file.
#   declarations for PARSING
my @vsearch_hits;
my @rearranged_output;
my $vsearch_hits = "vsearch_hits.txt";

#    for vsearch output file (.uc) in the temp directory, remove any line that begins with H (denotes hit) and push line to an array
print "\nparsing vsearch output . . .\n";
open (VSEARCH_OUTPUT, "<$temp_files_directory/vsearchfile.uc") or die ("$!");
my @vsearch_output_lines = <VSEARCH_OUTPUT>;
close VSEARCH_OUTPUT;
foreach my $vsearch_line (@vsearch_output_lines)
{
    if ($vsearch_line =~ m/^H/)
    {
        push (@vsearch_hits, $vsearch_line);
    }
}
foreach my $vsearch_hit (@vsearch_hits)
{
    $vsearch_hit =~ s/\s/\t/g;
    my $vsearch_line =  (join "\t", reverse split /\t/, $vsearch_hit) . "\n";
    push (@rearranged_output, $vsearch_line);
}

#   the array containing hits is then sorted according to ascending cluster number
my @sorted_vsearch_hits = sort { (split ' ', $a)[0] cmp (split ' ', $b)[0] } @rearranged_output;
print "\ncreating vsearch_hits.txt and cleaning up temp files . . .\n";
open (OUTFILE, ">$temp_files_directory/$vsearch_hits") or die ("$!");
foreach my $sorted_vsearch_hit (@sorted_vsearch_hits)
{
    print OUTFILE ("$sorted_vsearch_hit");
}
close (OUTFILE) or die ("$!");


############################################################################
### FASTQ LOOKUP ###
############################################################################
## The matched fasta reads (matched.reads.fa) from vsearch are used to look up the original fastq reads for use in the mapping stage.
#   the header of each matched read is used to retrieve the original fastq read
#   declarations for FASTQ LOOKUP
my @strippedreads;
my @total_number_of_input_reads;
print ("\nextracting matched reads from original fastq files . . .\n\n");
open (INFILE, "$temp_files_directory/matched.reads.fa") or die ("Can't find $temp_files_directory/matched.reads.fa: $!\n");
my @readstostrip = <INFILE>;
close (INFILE);
foreach my $readtostrip (@readstostrip)
{
    if ($readtostrip =~ /^>/)
    {
        $readtostrip =~ s/^./\@/;
        push (@strippedreads, $readtostrip);
    }
}

#   load input fastq files into hashes -> fastq headers = keys
foreach my $in_file (@opt_inputfiles)
{
    open (INFILE, $in_file) or die ("Can't find input files: $!\n");
    my $header = <INFILE>;
    close (INFILE);
    # have to use first few characters of header line/identifier for keys as "@" can be used as first char on quality score line as well header line
    my $header_first5_chars_identifier;
    ($header_first5_chars_identifier) = $header =~ /(^.{5}).*$/s;
    my %fastq;
    my $fastq_key = 'unlabeled';
    open (INFILE, $in_file) or die ("Can't find input files: $!\n");
    while (defined (my $line = <INFILE>))
    {
        chomp $line;
        if ($line =~ /(^$header_first5_chars_identifier.*)/)
        {
            $fastq_key = $1;
            next;
        }
        push @{$fastq{$fastq_key}}, $line;
    }
    close INFILE;
    
    #   cycle through array of fastq headers, if the header is in hash key table, print key and hash for that read in "matched.reads.fastq"
    foreach my $strippedread (@strippedreads)
    {
        chomp $strippedread;
        if (exists $fastq {$strippedread})
        {
            my $fastqlookup = ("$strippedread\n" . join "\n" => @{$fastq{$strippedread}});
            open (OUTFILE, ">>$temp_files_directory/matched.reads.fastq") or die ("Can't open outfile for fastq lookup: $!\n");
            print OUTFILE ("$fastqlookup\n");
            close OUTFILE;
        }
    }
    
    # count total number of reads from all input files (for use in calculating relative abundance of args)
    my $number_of_reads = keys %fastq;
    push (@total_number_of_input_reads, $number_of_reads);
    
    #   clear hash
    for (keys %fastq)
    {
        delete $fastq{$_};
    }
}
my $total_number_of_input_reads = 0;
for ( @total_number_of_input_reads )
{
    $total_number_of_input_reads += $_;
}

#   load matched.reads.fastq into hash -> fastq headers = keys, subsequent 3 lines = values
my $matchedreads = "matched.reads.fastq";
open (MATCHEDREADS, "$temp_files_directory/$matchedreads") or die ("Can't find ./temp_files/matched.reads.fastq: $!\n");
my $header = <MATCHEDREADS>;
my $header_first5_chars_identifier;
($header_first5_chars_identifier) = $header =~ /(^.{5}).*$/s;
close MATCHEDREADS;
open (MATCHEDREADS, "$temp_files_directory/$matchedreads") or die ("Can't find ./temp_files/matched.reads.fastq: $!\n");
my %matchedreadsfq;
my $matchedread_header = 'unlabeled';
while (defined(my $fqline = <MATCHEDREADS>))
{
    chomp $fqline;
    if ($fqline =~ /(^$header_first5_chars_identifier.*)/)
    {
        $matchedread_header = $1;
        next;
    }
    push @{$matchedreadsfq{$matchedread_header}}, $fqline;
}
close MATCHEDREADS;


############################################################################
### CLUSTER DIRECTORIES ###
############################################################################
## Each reference sequence that received hits during vsearch has a directory created in the temp directory containing the fastq reads that matched the reference sequence.
#   reopen vsearch hits, load columns into an array, if cluster number is duplicated condense the cluster to have all reads on one line
open(vsearch_HITS, "<$temp_files_directory/$vsearch_hits") or die ("$!");
my %info = ();
while (<vsearch_HITS>)
{
    my @columns = split /\t+/;
    if ( exists $info{ $columns[0] } )
    {
        push @{ $info{ $columns[0] }->{C2} }, $columns[1];
    }
    else {
        $info{ $columns[0] } = { C2 =>[ $columns[1] ] }
    }
}
close (vsearch_HITS) or die ("$!");

#   create new directory for every cluster that received hits, then to each directory add text file containing the fastq headings of sequences that hit that cluster
print "creating directories for reference sequence database hits and adding reads to directories . . .\n";
print "reference sequences with hits:";
my $cluster_counter = 1;
foreach my $c1 ( sort {$a cmp $b} keys %info )
{
    print "\n->$c1";
    my $reference_sequence_name = $c1;
    $reference_sequence_name =~ s/\//ForwardSlashHere/;
    my $cluster_directory = "$temp_files_directory/Cluster$cluster_counter";
    system ( "mkdir $cluster_directory" ) == 0 or die ( "Can't make directory $cluster_directory: $?.\n" );
    open ( OUTFILE, ">$cluster_directory/matchingreads.cluster_$reference_sequence_name.txt" ) or die ("$!");
    {
        print OUTFILE ( join ("\n",@{$info{$c1}->{C2}}), "\n", );
    }
    close OUTFILE;

#   for each cluster directory, create multifastq file of all reads in "matchingreads.cluster_no.txt" by comparing to hash keys of the matched.reads.fastq hash
    my $matchingreads = "$cluster_directory/matchingreads.cluster_$reference_sequence_name.txt";
    open (INFILE, $matchingreads) or die ("$!");
    my @array_of_matchingread_headers = <INFILE>;
    close INFILE;
    foreach my $arrayread (@array_of_matchingread_headers)
    {
        chomp $arrayread;
        $arrayread =~ s/^/\@/;
        if (exists $matchedreadsfq {$arrayread})
        {
            my $fastqrecord = ("$arrayread\n" . join "\n" => @{$matchedreadsfq{$arrayread}});
            open (OUTFILE, ">>$cluster_directory/matchingreads.cluster_$reference_sequence_name.fastq") or die ("Can't create outfile for cluster $c1 multifastq: $!\n");
            print OUTFILE ("$fastqrecord\n");
            close (OUTFILE);
        }
        else
        {
            die "cant find matched read '$arrayread' in '$matchedreads' for cluster $c1\n";    
        }
    }
    close (INFILE);
    $cluster_counter ++;
}

#   clear hash
for (keys %matchedreadsfq)
{
    delete $matchedreadsfq{$_};
}


############################################################################
### MAPPING ###
############################################################################
## SEAR cycles through the cluster directories created in the previous step and maps each set of reads within the directory to the reference sequence for that cluster.
#   open temp_files directory and load each cluster into an array
my @cluster_directories;
my $cluster_directory;
my $unmapped_directory_command = "mkdir $temp_files_directory/unsuccessfully_mapped_clusters";
system ( $unmapped_directory_command ) == 0 or die ("$!");
opendir (TEMPDIR, $temp_files_directory) or die "Couldn't open dir '$temp_files_directory': $!";
while ($cluster_directory = readdir (TEMPDIR))
{
    if ($cluster_directory =~ m/^Cluster/)
    {
        push (@cluster_directories, $cluster_directory);
    }
}
close (TEMPDIR);

#   cycle through array, opening each cluster directory and performing mapping
my @reference_args;
my @consensus_seqs;
my %abundance_values;
print "\n\nmapping reads to reference sequences for each cluster . . .\n\n";
foreach my $cluster (@cluster_directories)
{
    my $reference_seq_name;
    print "current cluster: $cluster";
    opendir (DIR, "$temp_files_directory/$cluster") or die "Can't open directory $cluster: $!";
    while (my $cluster_file = readdir(DIR))
    {
       if ($cluster_file =~ m/matchingreads.cluster_\S+\.fastq/)
       {
           $cluster_file =~ m/matchingreads.cluster_(\S+)\.fastq/;
           $reference_seq_name = $1;
       }
    }
    closedir(DIR);
    
    # retrieve reference sequence for cluster
    my $bwa_index;
    my $bwa_index_outfile = "bwa_index.fasta";
    my $data = `cat $opt_database`;
    print ">$reference_seq_name\n\n";
    chomp $reference_seq_name;
    $reference_seq_name =~ s/ForwardSlashHere/\//;
    if ($data =~ /(>\Q$reference_seq_name\E[^>]+)/s)
    {
        $bwa_index = $1;
    }
    else
    {
        print "can't retrieve reference sequence for cluster\n";
        exit;
    }
    open (OUTFILE, ">$temp_files_directory/$cluster/bwa_index.fasta") or die ("Can't find reference ARG for $cluster: $!\n");
    print OUTFILE ("$bwa_index\n");
    close (OUTFILE);
    
    #   run bwa on the cluster
    my $bwa_command_1 = "bwa index $temp_files_directory/$cluster/bwa_index.fasta";
    my $bwa_command_2;
    if ($opt_ascii == 33)
    {
        #$bwa_command_2 = "bwa aln -t $opt_threads $temp_files_directory/$cluster/bwa_index.fasta $temp_files_directory/$cluster/*.fastq > $temp_files_directory/$cluster/reads.sai -q $opt_quality";
        $bwa_command_2 = "bwa mem -t $opt_threads $temp_files_directory/$cluster/bwa_index.fasta $temp_files_directory/$cluster/*.fastq > $temp_files_directory/$cluster/aln.sam";
    }
    if ($opt_ascii == 64)
    {
        $bwa_command_2 = "bwa aln -t $opt_threads $temp_files_directory/$cluster/bwa_index.fasta $temp_files_directory/$cluster/*.fastq > $temp_files_directory/$cluster/reads.sai -q $opt_quality -I";
    }
    #my $bwa_command_3 = "bwa samse $temp_files_directory/$cluster/bwa_index.fasta $temp_files_directory/$cluster/reads.sai $temp_files_directory/$cluster/*.fastq > $temp_files_directory/$cluster/aln.sam";
    system (" $bwa_command_1 ") == 0 or die ( "Error in bwa command: $?.\n" );
    system (" $bwa_command_2 ") == 0 or die ( "Error in bwa command: $?.\n" );
    #system (" $bwa_command_3 ") == 0 or die ( "Error in bwa command: $?.\n" );
    
    #   convert alignment to bam file
    my $samtools_command_1 = "samtools view -bS $temp_files_directory/$cluster/aln.sam > $temp_files_directory/$cluster/aln.bam";
    my $samtools_command_2 = "samtools sort $temp_files_directory/$cluster/aln.bam $temp_files_directory/$cluster/aln.sorted";
    my $samtools_command_3 = "samtools index $temp_files_directory/$cluster/aln.sorted.bam";
    system (" $samtools_command_1 ") == 0 or die ( "Error in samtools command: $?.\n" );
    system (" $samtools_command_2 ") == 0 or die ( "Error in samtools command: $?.\n" );
    system (" $samtools_command_3 ") == 0 or die ( "Error in samtools command: $?.\n" );
    
    #   calculate number of mapped reads
    my $mappedreads = `samtools view -c -F 4 $temp_files_directory/$cluster/aln.sorted.bam`;
    print "\nnumber of reads mapped: $mappedreads\n";
    
    #   if no reads map to reference, move cluster directory to unsuccessful mapping directory
    if ($mappedreads == 0)
    {
        my $move_directory = "mv $temp_files_directory/$cluster $temp_files_directory/unsuccessfully_mapped_clusters/";
        system ( $move_directory ) == 0 or die ("$!");
    }
    
    #   if reads map to reference, create a consensus sequence and annotate as putative ARG if within annotation thresholds
    else
    {
        #   create pileup for the alignment and eactract consensus sequence
        my $samtools_command_4 = "samtools mpileup -EAuf $temp_files_directory/$cluster/bwa_index.fasta $temp_files_directory/$cluster/aln.sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq > $temp_files_directory/$cluster/consensus_$cluster.fastq";
        system ( "$samtools_command_4" )  == 0 or die ( "Error in mpileup command: $?.\n" );
        my $consensusfq = "$temp_files_directory/$cluster/consensus_$cluster.fastq";
        
        #   calculate sequence length of cluster representative sequence
        my $index = "$temp_files_directory/$cluster/bwa_index.fasta";
        open (BWA_INDEX, "<$index") or die "Could not open $cluster representative sequence file: $!";
        my $ref_seq = do { local $/; <BWA_INDEX> };
        close BWA_INDEX;
        $ref_seq =~ s/>$cluster\n//s;
        $ref_seq =~ s/\n//sg;
        my $ref_seq_length = length($ref_seq);
        
        #   calculate breadth of coverage
        open (CONFILE, "<$consensusfq") or die "Could not open consensus sequence file: $!";
        my $consensus_file = do { local $/; <CONFILE> };
        close CONFILE;
        $consensus_file =~ s/\@$cluster\n//;
        $consensus_file =~ s/\n\+\n/$cluster/;
        $consensus_file =~ s/$cluster.*//s;
        $consensus_file =~ s/\n//sg;
        my $ncount = ($consensus_file =~ tr/n//);
        my $consensus_length = length($consensus_file);
        my $basecount = ($consensus_length - $ncount);
        my $ref_coverage = (($basecount/$ref_seq_length) * 100);

        #   if breadth of coverage is below $opt_coverage threshold, omit cluster from results
        if ($ref_coverage < $opt_coverage)
        {
            print "coverage of reference sequence less than $opt_coverage%\n\n\n";
            my $move_directory = "mv $temp_files_directory/$cluster $temp_files_directory/unsuccessfully_mapped_clusters/";
            system ( $move_directory ) == 0 or die ("$!");
        }
        else
        {
            #   create multifastq of consensus sequences
            open (INFILE, $consensusfq);
            open (OUTFILE, ">>$opt_results/$multifastq") or die "can't write to $opt_results/$multifastq ($!)";
            print OUTFILE <INFILE>;
            close INFILE;
            close OUTFILE;

            #   calculate abundance value for cluster -> this is the proportion of reads from input file that aligned/length of the reference sequence
            my $proportion = ($mappedreads/$total_number_of_input_reads);
            my $abundance_value = ($proportion/$ref_seq_length);
            print "abundance value of gene cluster: $abundance_value\n\n";

            #   push abundance value to hash, using reference sequence name as key. Keep an array of reference sequence names for sorting later
            $abundance_values{$reference_seq_name} = $abundance_value;
            push (@consensus_seqs, $reference_seq_name);
        }
    }
}

#   calculate percentage relative abundances of ARGs by dividing each abundance value by the sum of the total abundance values and multiplying by 100
my $total_of_abundance_values = sum values %abundance_values;
map { $_ /= $total_of_abundance_values } values %abundance_values;
map { $_ *= 100 } values %abundance_values;

#   check if no reads passed this step
finddepth(\&filecheck, "$opt_results/$multifastq");


############################################################################
### BLAST PUTATIVE ARGs ###
############################################################################
## For each putative ARG, blast the sequence against multiple databases to confirm annotations
#   if putative ARGs have been found:
if (%abundance_values)
{
    #   sort abundance values numberically
    my @abundance_keys = sort { $abundance_values{$b} <=> $abundance_values{$a} } keys %abundance_values;
    
    #   convert fastq of consensus seqs to multifasta file for blast input
    my $seqtk_command = `seqtk seq -a $opt_results/$multifastq > $opt_results/BLAST_QUERY_SEQUENCES.fasta`;
    my $blast_multifasta = "$opt_results/BLAST_QUERY_SEQUENCES.fasta";
    local $/ = ">";
    open INFILE, $blast_multifasta;
    my @blast_multifasta = <INFILE>;
    close INFILE;
    chomp @blast_multifasta;
    shift @blast_multifasta;
    foreach my $fasta_entry (@blast_multifasta)
    {
        $fasta_entry =~s/^/>/;
    }
    local $/ = "\n";

    print "\n################################################\n";
    print "Parsing output, preparing for database searching\n";
    print "################################################\n\n";
    print "Extracting consensus sequences from read mapping data . . .\n\n";
    print "Searching the following databases for antimicrobial resistance using consensus sequences . . .\n\n#...NCBI nr database\n\nrunning blastn and blastx remotely . . .\n";
    
    my $blast_query_fasta = "BLAST_QUERY_SEQUENCES.fasta";
    my $query_db_blast = "$opt_results/BLASTn_QUERYdb_RESULTS.html";
    my $ncbi_blastn = "$opt_results/BLASTn_NCBI_RESULTS.html";
    my $ncbi_blastx = "$opt_results/BLASTx_NCBI_RESULTS.html";
    my $rac_blast = "$opt_results/BLASTn_RAC_RESULTS.html";
    
    #   remote blastn + blastx at NCBI
    my $blastn_command = "blastn -db nt -query $opt_results/BLAST_QUERY_SEQUENCES.fasta -out $ncbi_blastn -html -remote";
    print "$blastn_command\n";
    system ("$blastn_command") == 0 or die "Can't run blastn: $!";
    my $blastx_command = "blastx -db nr -query $opt_results/BLAST_QUERY_SEQUENCES.fasta -out $ncbi_blastx -html -num_alignments 1 -remote";
    print "$blastx_command\n";
    system ("$blastx_command") == 0 or die "Can't run blastx: $!";

    #   local blast to query database used in SEAR (default= ARGannot)
    print "\n\n#...blastn against query database $opt_database";
    my $makeblastdb1 = "makeblastdb -in $opt_database -dbtype nucl -out $temp_files_directory/USERdb";
    my $USERdb_blastn_command = "blastn -db $temp_files_directory/USERdb -query $opt_results/BLAST_QUERY_SEQUENCES.fasta -out $query_db_blast -num_alignments 1 -html -num_threads $opt_threads";
    system ( "$makeblastdb1" )  == 0 or die ( "Error in $makeblastdb1: $!\n" );
    system ( "$USERdb_blastn_command" )  == 0 or die ( "Error in $USERdb_blastn_command: $!\n" );

    #   local blast to RAC
    print "\n\n#...the Repository of Antibiotic resistance Cassettes (RAC):\n\ndownloading most recent data from the RAC . . .\n\n";
    my $RACfasta = Sear::getRACupdate($temp_files_directory, $RAC_web_page);
    my $RACfastafile = "$temp_files_directory/RAC.fasta";
    open OUTFILE, ">$RACfastafile" or die "$!\n";
    print OUTFILE $RACfasta;
    close OUTFILE or die "$!\n";
    if (-s $RACfastafile)
    {
        print "\ncreating new blast database from RAC data and BLASTing consensus sequences with blastn . . .";
        my $makeblastdb2 = "makeblastdb -in $RACfastafile -dbtype nucl -out $temp_files_directory/RACdb";
        my $RAC_blastn_command = "blastn -db $temp_files_directory/RACdb -query $opt_results/BLAST_QUERY_SEQUENCES.fasta -out $rac_blast -num_alignments 1 -html -num_threads $opt_threads";
        system ( "$makeblastdb2" )  == 0 or die ( "Error in $makeblastdb2: $!\n" );
        system ( "$RAC_blastn_command" )  == 0 or die ( "Error in $RAC_blastn_command: $!\n" );
    }

    
    #   ARDB blast
    print "\n\n#...the Antibiotic Resistance Genes Database (ARDB):\n\nsubmitting sequences for blastn @ http://ardb.cbcb.umd.edu/ . . .\n\n";
    my %ardb_links;
    foreach my $blast_multifasta (@blast_multifasta)
    {
        my ($ARDB, $ARDB_key) = Sear::getARDBdata($blast_multifasta);
        push @{$ardb_links{$ARDB_key}}, $ARDB;
    }
    

    ############################################################################
    ### SEAR OUTPUT ###
    ############################################################################
    ## The SEAR program creates the output files and results page (in html format), complete with links out to online databases.
    #   create raw data file
    my $raw_data_outfile = "$opt_results/SEARdata.xls";
    open (OUTFILE, ">>$raw_data_outfile") or die "can't write to $raw_data_outfile: $!";
    print OUTFILE "$opt_results\n#\tREFERENCE SEQUENCE MATCHED\tRELATIVE ARG ABUNDANCE/%\tblastn ARDB\n";
    my $result_counter = 1;
    foreach my $reference_seq_name (@abundance_keys)
    {
        foreach my $consensus_seq_name (@consensus_seqs)
        {
            if ($consensus_seq_name eq $reference_seq_name)
            {
                print OUTFILE "$result_counter\t$reference_seq_name\t$abundance_values{$reference_seq_name}\t";
            }
        }
        foreach my $ardb_links_key (keys %ardb_links)
        {
            if ($ardb_links_key eq $reference_seq_name)
            {
                print OUTFILE "@{$ardb_links{$ardb_links_key}}\n";
                $result_counter ++;
                next;
            }
        }
    }
    close OUTFILE;
    open (CONSESUS_SEQS, "$opt_results/$blast_query_fasta") or die ("Can't open $opt_results/$blast_query_fasta: $!\n");
    my %consensus_seqs;
    my $matchedread_header = "";
    while (defined(my $line = <CONSESUS_SEQS>))
    {
        if ($line =~ /^>(.*)\n/)
        {
            $matchedread_header = $1;
            next;
        }
        chomp $line;
        push @{$consensus_seqs{$matchedread_header}}, $line;      
    }
    close CONSESUS_SEQS;
    open INFILE,  "$raw_data_outfile" or die "Can't open $raw_data_outfile: $!";
    my @raw_outfile = <INFILE>;
    close INFILE;
    
    #   create graph
    my $graph = Sear::graphGEN(@raw_outfile);
    open OUTFILE, ">$temp_files_directory/r_barchart_script.R";
    print OUTFILE "$graph";
    close OUTFILE;
    my $html_files_directory = "$opt_results/html_files";
    unless (-e $html_files_directory or mkdir $html_files_directory) {
      die "Unable ot make html_files directory\n";
    }
    my $rcommand = "Rscript $temp_files_directory/r_barchart_script.R";
    system ( $rcommand ) == 0 or die ("Can't run r script $!");
    
    #   create html files
    my $results_redirect = "./FILES/html_files/";
    my $html_results_page = `cat $SEAR_DATA_directory/references/html_templates/RESULTS_template.html`;
    my $sear_settings_line = "<td>...end of SEAR settings...</td>";
    my $sear_settings_entry = "<tr>\n<td>" . $opt_ascii . "</td>\n<td>" . $opt_length . "</td>\n<td>" . $opt_quality . "</td>\n<td>" . $opt_filter_reads . "</td>\n<td>" . $opt_clustering_identity . "</td>\n<td>" . $opt_coverage . "</td>" . "\n</tr>\n";
                if ($html_results_page =~ m/\Q$sear_settings_line\E/)
                {
                    $html_results_page =~ s/\<tr\>\s+\<td\>\.\.\.end\sof\sSEAR\ssettings\.\.\.\<\/td\>\s+\<\/tr\>/$sear_settings_entry\<tr\>/;
                }
    my $graph_insert_line = "<img alt=\"results graph\" src=\"";
    my $html_graph_entry = $graph_insert_line . $results_redirect . "/graph.png";
    if ($html_results_page =~ m/\Q$graph_insert_line\E/)
    {
        $html_results_page =~ s/\Q$graph_insert_line\E/$html_graph_entry/;
    }
    open OUTFILE, ">$raw_data_outfile" or die "Can't open $raw_data_outfile: $!";
    foreach my $raw_outfile_line (@raw_outfile)
    {
        foreach my $key(keys %consensus_seqs)
        {
            chomp $raw_outfile_line;
            if ($raw_outfile_line =~ m/\Q$key\E/)
            {
                my $query_db_blastn_result = Sear::blastPARSER($key,$query_db_blast);
                my $blastn_result = Sear::blastPARSER($key,$ncbi_blastn);
                my $blastx_result = Sear::blastPARSER($key,$ncbi_blastx);
                my $rac_blastn_result = Sear::blastPARSER($key,$rac_blast);

                my $key_no_specials = $key;
                $key_no_specials =~ s/[^a-z0-9]/\_/gi;
                
                open BLAST_OUTFILE, ">$opt_results/html_files/$key_no_specials.querydb.blastn.html";
                print BLAST_OUTFILE $query_db_blastn_result;
                close BLAST_OUTFILE;
                
                open BLAST_OUTFILE, ">$opt_results/html_files/$key_no_specials.ncbi.blastn.html";
                print BLAST_OUTFILE $blastn_result;
                close BLAST_OUTFILE;
                
                open BLAST_OUTFILE, ">$opt_results/html_files/$key_no_specials.ncbi.blastx.html";
                print BLAST_OUTFILE $blastx_result;
                close BLAST_OUTFILE;
                
                open BLAST_OUTFILE, ">$opt_results/html_files/$key_no_specials.rac.blastn.html";
                print BLAST_OUTFILE $rac_blastn_result;
                close BLAST_OUTFILE;
                
                open FASTA_OUTFILE, ">$opt_results/html_files/$key_no_specials.consensus.fasta.html";
                print FASTA_OUTFILE "<HTML>\n<TITLE>$key</TITLE>\n" . '<BODY BGCOLOR="#FFFFFF" LINK="#0000FF" VLINK="#660099" ALINK="#660099">' . "\n<PRE>\n" . "@{$consensus_seqs{$key}}" . "\n</PRE>\n</BODY>\n</HTML>";
                close FASTA_OUTFILE;
                
                print OUTFILE "$raw_outfile_line\t@{$consensus_seqs{$key}}\n";
                my ($results_counter, $ref_match, $abundance, $ardb_link) = split /\t/, $raw_outfile_line;
                my $query_db_blastn_result_html = "$results_redirect$key_no_specials.querydb.blastn.html";
                my $blastn_result_html = "$results_redirect$key_no_specials.ncbi.blastn.html";
                my $blastx_result_html = "$results_redirect$key_no_specials.ncbi.blastx.html";
                my $rac_blastn_result_html = "$results_redirect$key_no_specials.rac.blastn.html";
                my $consensus_sequence = "$results_redirect$key_no_specials.consensus.fasta.html";
                my $end_of_results_line = "<td>...end of results...</td>";
                my $html_table_entry = "<tr>\n<td>" . $results_counter . "</td>\n<td>" . $key . "</td>\n<td>" . $abundance . "</td>\n<td>" . '<p><a class="btn btn-info btn-xs" href="' . $consensus_sequence . '" role="button">VIEW SEQUENCE&raquo;</a></p>' . "</td>\n<td>" . '<p><a class="btn btn-info btn-xs" href="' . $query_db_blastn_result_html . '" role="button">VIEW BLASTn RESULT&raquo;</a></p>' . "</td>\n<td>" . '<p><a class="btn btn-info btn-xs" href="' . $ardb_link . '" role="button">VIEW BLASTn RESULT&raquo;</a></p>' . "</td>\n<td>" . '<p><a class="btn btn-info btn-xs" href="' . $rac_blastn_result_html . '" role="button">VIEW BLASTn RESULT&raquo;</a></p>' . "</td>\n<td>" . '<p><a class="btn btn-info btn-xs" href="' . $blastn_result_html . '" role="button">VIEW BLASTn RESULT&raquo;</a></p>' . "</td>\n<td>" . '<p><a class="btn btn-info btn-xs" href="' . $blastx_result_html . '" role="button">VIEW BLASTx RESULT&raquo;</a></p>' . "</td>" . "\n</tr>\n";
                if ($html_results_page =~ m/\Q$end_of_results_line\E/)
                {
                    $html_results_page =~ s/\<tr\>\s+\<td\>\.\.\.end\sof\sresults\.\.\.\<\/td\>\s+\<\/tr\>/$html_table_entry\<tr\>\n$end_of_results_line\n\<\/tr\>/;
                }
            }
        }
    }
    close OUTFILE;
    
    #   make results flatfile (may not be useful...--> leftover function from web version of SEAR)
    my $make_flat_file = "mkdir $opt_results/flat_file_and_blast_files && mv $opt_results/*.html $opt_results/flat_file_and_blast_files/ && mv $opt_results/SEARdata.xls $opt_results/flat_file_and_blast_files/ && tar -cvf $opt_results/flat_file.tar -C $opt_results/flat_file_and_blast_files .";
    system ( "$make_flat_file" ) == 0 or die ("Can't make flat file: $!\n");
    my $flat_file_insert_line = "href=\"download_flat_file\"";
    my $flat_file_entry = "href=\"./FILES/flat_file.tar\"";
    if ($html_results_page =~ m/\Q$flat_file_insert_line\E/)
    {
        $html_results_page =~ s/\Q$flat_file_insert_line\E/$flat_file_entry/;
    }

    #   write results to html file
    open RESULTS_PAGE_OUTFILE, ">$opt_results/SEAR_RESULTS_$temp_file_name.html" or die "Can't generate results page: $!";
    print RESULTS_PAGE_OUTFILE "$html_results_page";
    close RESULTS_PAGE_OUTFILE or die "Can't generate results page: $!";

    #   tidy basespace output files
    my $tidy_dir = "$opt_results/FILES";
    eval { make_path($tidy_dir) };
    if ($@) {
      print "Couldn't create $tidy_dir: $@";
    }
    my $move_files = "rm -r $opt_results/flat_file_and_blast_files && mv $opt_results/html_files $opt_results/BLAST_QUERY_SEQUENCES.fasta $opt_results/consensus_seqs.fq $opt_results/flat_file.tar $tidy_dir && cp -r /SEAR/SEAR_DATA/references $tidy_dir/ && cp /SEAR/README_basespace_results $opt_results/";
    system ( "$move_files" ) == 0 or die ("Can't tidy files: $!\n");

    #   tidy results
    print "\n\n################################################\n";
    print "SEAR Finished\n";
    print "################################################\n\n";
    print "Results can be viewed by downloading the output file directory and opening the SEAR results web page\n";
    my $SEAR_finish_time = Time::HiRes::tv_interval($start_time);
    print "Elapsed time:\t$SEAR_finish_time\n\n";
    system ("rm -r $temp_files_directory") == 0 or die "Can't delete temporary data directory on scratch disk";
    exit(0);
}


############################################################################
### SEAR ALTERNATIVE OUTPUT ###
############################################################################
## If no results have been found.
else
{
    print "\n\n\n\n################################################\n";
    print "SEAR Finished\n";
    print "################################################\n\n";
    print "\nSearch Engine for Antimicrobial Resistance did not find any antimicrobial resistance genes.\n\nTry adjusting the stringency settings.\n\n";
    my $SEAR_finish_time = Time::HiRes::tv_interval($start_time);
    print "Elapsed time:\t$SEAR_finish_time\n\n";
    system ("rm -r $temp_files_directory && rm $opt_results/$multifastq && touch $opt_results/NO_ARGs_FOUND") == 0 or die "Can't move log file and delete temp file folder from server: $!";
    exit(0);
}

sub filecheck
{
    if (-f $_ && -z $_)
    {
        print "\n\n\n\n################################################\n";
        print "SEAR Finished\n";
        print "################################################\n\n";
        print "Search Engine for Antimicrobial Resistance did not find any antimicrobial resistance genes.\n\nTry adjusting the stringency settings.\n\n";
        my $SEAR_PROGRAM_time = Time::HiRes::tv_interval($start_time);
        print "Elapsed time:\t$SEAR_PROGRAM_time\n\n";
        system ("rm -r $temp_files_directory && rm $opt_results/$multifastq && touch $opt_results/NO_ARGs_FOUND") == 0 or die "Can't move log file and delete temp file folder from server: $!";
        exit(0);
    }
    else
    {
        my $filecheck_time = Time::HiRes::tv_interval($start_time);
        print "\n\n\n################################################\nSEAR UPDATE:\nSEAR proceeding as normal\nElapsed time=\t$filecheck_time\n################################################\n";
    }
}
