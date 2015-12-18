package Sear;

use strict;
use warnings;
use LWP::Simple;

sub getDBdate
{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $DBtimestamp = sprintf ( "%02d.%02d_%02d:%02d:%04d",$hour,$min,$mday,$mon+1,$year+1900);
    return $DBtimestamp;
}

sub getRACupdate
{
my $temp_files_directory = shift;
my $web_page = shift;
my $web_page_grab = "$temp_files_directory/rac_feature_list";
my $curl = "curl $web_page > $web_page_grab";
system ( $curl ) == 0 or die "Can't grab $web_page using curl: $!\n";
my $time = getDBdate();
my $outfile = "$temp_files_directory/updated_RAC_database";
my $infile = $web_page_grab;
my @infile_array;
open INFILE, $infile or die "Can't open $infile\n";
@infile_array = <INFILE>;
close INFILE;
my $dead_link = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nucleotide&list_uids=&from=0&to=0";
my $link;
my $previous_link = "0";
my @link_array;
my $name;
my %link_name_hash;
my @accession_list;
foreach my $line (@infile_array)
{
    if ($line =~ m/\<td\>\<a\s+target=\"_blank\"\s+href\=\"\S+\"\>(\S+)\<\/a\>\<\/td\>/)
    {
        $name = $1;
        $line =~ m/\<td\>\<a\s+target=\"_blank\"\s+href\=\"(\S+)\"/;
        $link = $1;
        if ($link eq $dead_link)
        {
            next;
        }
        if ($link eq $previous_link)
        {
            next;
        }
        push (@link_array, $link);
        push @{$link_name_hash{$link}}, $name;
        $previous_link = $link;
        next;
    }
}
my $accession;
my $rac_data = `cat $infile`;
my %link_accession_hash;
foreach my $array_link (@link_array)
{
    chomp $array_link;
    if ($rac_data =~ m/\<td\>\<a\s+target=\"_blank\"\s+href\=\"\Q$array_link\E\"\>\S+\<\/a\>\<\/td\>\s+\<td\>(\S+)\<\/td\>/s)
    {
        $accession = $1;
        push @{$link_accession_hash{$array_link}}, $accession;
    }
    else
    {
        push @{$link_accession_hash{$array_link}}, "no accession for link";
    }
}
foreach (keys %link_name_hash) {
    if (exists $link_accession_hash {$_})
    {
        my $entry = ("@{$link_accession_hash{$_}}\t" . join "\t" => $_);
        open OUTFILE, ">>$outfile" or die "Can't open $outfile: $!";
        print OUTFILE "@{$link_name_hash{$_}}\t$entry\n";
        close OUTFILE;
        push (@accession_list, @{$link_accession_hash{$_}});
    }
}
system ( "rm $web_page_grab" ) == 0 or die "Can't remove file $web_page_grab: $!";
foreach my $accession_number (@accession_list)
{
    $accession_number .= "[accn]";
}
my $query = join('+OR+',@accession_list);
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=nucleotide&term=$query&usehistory=y";
my $output = get($url);
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $url2 = $base . "efetch.fcgi?db=nucleotide&query_key=$key&WebEnv=$web";
$url2 .= "&rettype=fasta&retmode=text";
my $RACfasta = get($url2);
return $RACfasta;
}

sub getARDBdata
{
    my $blast_multifasta = shift;
    $blast_multifasta =~ m/^>(.+)\n/;
    my $header = $1;
    $blast_multifasta =~ s/^>.+\n//;
    $blast_multifasta =~ s/\n//gs;
    my $ardb_curl_command = 'curl "http://ardb.cbcb.umd.edu/cgi/blast/nblast.cgi?PROGRAM=blastn&DATALIB=ardb_dna&ALIGNMENT_VIEW=0&PERC_IDENT=80&EXPECT=0.0001&SEQUENCE=' . $blast_multifasta . '&MAT_PARAM=BLOSUM62%09+11%09+1&GENETIC_CODE=Standard+%281%29&DB_GENETIC_CODE=Standard+%281%29&OOF_ALIGN=0&OTHER_ADVANCED=&press=Search"';
    my $ardb_data = `$ardb_curl_command`;
    my $ardb_accession;
    my $ardb_top_hit;
    if ($ardb_data =~ m/Sequences\sproducing\ssignificant\salignments:\s+\(bits\)\sValue\s+(\S+)\<\/a\>/m)
    {
        $ardb_accession = $1;
        $ardb_top_hit = "http://ardb.cbcb.umd.edu/cgi/search.cgi?db=L&field=ni&term=" . $ardb_accession;
    }
    else
    {
        $ardb_top_hit = "no ardb result found";
    }
    return ($ardb_top_hit, $header);
}

sub blastPARSER
{
    my ($query,$blast_file) = @_;
    my $entry = "<HTML>\n<TITLE>$query\t$blast_file</TITLE>\n" . '<BODY BGCOLOR="#FFFFFF" LINK="#0000FF" VLINK="#660099" ALINK="#660099">' . "\n<PRE>\n";
    my $load_blast_file = `cat $blast_file`;
    if ($load_blast_file =~ m/(\Q$query\E.*?)Effective\ssearch\sspace\sused:/s)
    {
        $entry .= $1;
    }
    $entry .= "\n</PRE>\n</BODY>\n</HTML>";
    return $entry;
}

sub graphGEN
{
    my @pipelineoutput = @_;
    my $opt_results;
    my @barchart_data;
    my @barchart_names;
    my $other = 0;
    foreach my $line (@pipelineoutput)
    {
        chomp $line;
        if ($line =~ m/^\//)
        {
            $opt_results = $line . "/html_files/graph.png";
        }
        if ($line =~ m/^#/) {next};
        if ($line =~ m/^\d+\s+(\S+)/)
        {
                my $data_chunk = $1;
                $line =~ /^\d+\s+\S+\s+(\S+)/;
                $line = $1;
                if ($line < 1)
                {
                    $other += $line;
                }
                else
                {
                    my $barchart_name = "\"$data_chunk\"";
                    push (@barchart_names, $barchart_name);
                    my $barchart_data = "$line";
                    push (@barchart_data, $barchart_data);
                }
        }
    }
    my $r_barchart_data;
    my $r_barchart_names;
    if ($other > 0)
    {
        $r_barchart_data = join(',', @barchart_data, $other);
        $r_barchart_names = join(',', @barchart_names, '"combined at <1%"')
    }
    else
    {
        $r_barchart_data = join(',', @barchart_data);
        $r_barchart_names = join(', ', @barchart_names);
    }
    
    my $outfile = "options(device = \"png\")\npng(\"$opt_results\",type=\"cairo\",width=800,height=800,units=\"px\",res=NA,pointsize=12)\npar(mar=c(20,4,4,2))\nx=c($r_barchart_data)\nnames(x)=c($r_barchart_names)\nbarplot(x, main=\"RESULTS\", ylab=\"RELATIVE ABUNDANCE (%)\", las=3)\n";
    return $outfile;
}
1;