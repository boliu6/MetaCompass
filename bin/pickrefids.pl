#!/usr/bin/perl

#############################################
#
# Program: Pick reference genomes based on
#          metaphyler output
#
# Author: Bo Liu, boliu@umiacs.umd.edu
#
# Fri Mar  9 23:24:30 EST 2012
#
#############################################

use strict;
use warnings;
use FindBin qw($Bin);


#----------------------------------------#
# read command line options
#----------------------------------------#
my $blast = "";
my $part = 0;
if (scalar @ARGV == 2) {
    ($blast, $part) = @ARGV;
} else {
    Usage();
}
#----------------------------------------#


#----------------------------------------#
# information of reference genomes
#----------------------------------------#
my %seq2tid = ();
my %tid2seqs = ();
my %tid2name = ();
my %tid2sp = ();
open(FH, "$Bin/../refseq/tid2par.tab") or die("Could not open $Bin/tid2par.tab\n");
foreach my $line (<FH>) {
    chomp $line;
    my ($seq, $tid, $sp, $ge, $name) = split("\t", $line);
    $seq2tid{$seq} = $tid;
    $tid2seqs{$tid}{$seq} = 1;
    $tid2name{$tid} = $name;
    if ($sp eq "NA") { $sp = $tid;}
    $tid2sp{$tid} = $sp;
}
close FH;
#----------------------------------------#


#----------------------------------------#
# abundance for each reference genome
#----------------------------------------#
my %tid2num = ();
# read tab delimited file
open(FH, "$blast") or die("Could not open $blast.\n");
foreach my $line (<FH>) {
    chomp $line;
    my ($qid, $rid, $pct, $hspl) = split("\t", $line);
    if ($pct < 90 || $hspl < 60) { next;}
    $rid =~ /^(\S+)\_\d+\_\d+$/;
    $rid = $1;
#    print "$rid\n";
    my $tid = $seq2tid{$rid};
    $tid2num{$tid}++;
}
close FH;
#----------------------------------------#


#----------------------------------------#
# based on the abundance of genomes,
# pick reference
#----------------------------------------#
my $total = 0;
my %sps = ();
foreach my $tid (sort {$tid2num{$b} <=> $tid2num{$a}} keys %tid2num) {

    my $num = $tid2num{$tid}*$part;
    
    if ($num < 60) { last;}    # if abundance < 60 then ignore
    if ($total > 120) { last;} # if more than 120 genomes have been used then stop
    
    my $sp = $tid2sp{$tid};
    if (exists $sps{$sp}) {
	if ($sps{$sp} >= 5) { next;}  # for a single species, allow max 5 reference genomes
    }
    
    $sps{$sp}++;
    $total++;
    foreach my $seq (keys %{$tid2seqs{$tid}}) {
	print "$seq\t$tid\t$num\t$tid2name{$tid}\n";
    }
}



exit;


sub Usage {
    die("
Usage:
       perl pickrefids.pl <blast file> <partition>

Options:
       <blast>        Query sequences in FASTA format to be classified.
       <partition>    Suppose you used a quarter of total data to estimate the
                      taxonomic profile. Then partition number is 4.

Output:
       A list of NCBI genome IDs, which will be used as reference genomes.

Contact:
       Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}

# todo
# pick reference genomes based on the estimation of depth of coverage
# allow users to tune this parameter, right now it's hard coded
