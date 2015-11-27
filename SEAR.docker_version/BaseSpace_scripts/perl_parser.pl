#!/usr/bin/perl
## Written by Will Rowe (wpmr2@cam.ac.uk)
## This script parsers the App input (json format), creates a copy of the SEAR program on the BaseSpace scratch disk and initates a SEAR run using a AWS machine (cc2.8xlarge) using 32 vCPUs.

use strict;
use warnings;                   
use JSON;
use Getopt::Long;

my $jsonfile = "/data/input/AppSession.json";

open INFILE, $jsonfile or die "couldn't open json file: $jsonfile";
my $json_data = <INFILE>;
close INFILE;
die "Could not open json file!" unless defined $json_data;

# Decode the entire JSON
my $decoded_json = decode_json( $json_data );

# Assign parameters to variables
my $length_cutoff = $decoded_json->{Properties}->{Items}->[4]->{Content};
my $QS_cutoff = $decoded_json->{Properties}->{Items}->[8]->{Content};
my $contamination_removal = $decoded_json->{Properties}->{Items}->[9]->{Content};
my $clustering_stringency = $decoded_json->{Properties}->{Items}->[3]->{Content};
my $annotation_stringency = $decoded_json->{Properties}->{Items}->[1]->{Content};

# Get sample data
my $sample_ID = $decoded_json->{Properties}->{Items}->[12]->{Items}->[0]->{Id};
my $sample_data_DIR = "/data/input/samples/" . $sample_ID . "/";
my @input_files;
opendir (DIR, $sample_data_DIR) or die $!;
while (my $file = readdir(DIR))
{
	next unless (-f "$sample_data_DIR/$file");
	print "$file\n";
	if ($file =~ m/\.gz$/)
	{
		push(@input_files, $file);
	}
        if ($file =~ m/\.fq$/)
        {
                push(@input_files, $file);
        }
        if ($file =~ m/\.fastq$/)
        {
                push(@input_files, $file);
        }
}
closedir(DIR);

# Check just one or two files are present as input files
my $num_files = scalar @input_files;
if ($num_files == 1)
{
	$num_files ++;
}
die "Can't find FASTQ files - needs 1 or 2 fastq files (either .fq/.fastq or .fq.gz/.fastq.gz)!\n" if ($num_files != 2);
for my $x (@input_files)
{
	$x =  "-i " . $sample_data_DIR . $x ;
}
my $input_file_list = join(" ", @input_files);

# App output files
my $project_output = $decoded_json->{Properties}->{Items}->[13]->{Items}->[0]->{Id};

# Move a copy of SEAR to the basespace scratch directory
my $sear_scratch = "cp -r /SEAR /data/scratch/";
system ( $sear_scratch ) == 0 or die $!;

# Initiate SEAR
my $SEAR = "/usr/bin/perl /data/scratch/SEAR/basespace_SEAR.pl $input_file_list -lc $length_cutoff -qc $QS_cutoff -cc $annotation_stringency -ci $clustering_stringency -dd /data/scratch/SEAR/SEAR_DATA -t 32 -po $project_output";
print "SEAR run with the following command:\n\n$SEAR\n\n";
system ( $SEAR ) == 0 or die $!;