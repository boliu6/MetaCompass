#!/usr/bin/perl

#############################################
#
# Program: Metagenomic comparative assembly
#
# Author: Bo Liu, boliu@umiacs.umd.edu
#
# Thu Mar 22 15:32:47 EDT 2012
#
#############################################

use strict;
use warnings;
use FindBin qw($Bin);


#----------------------------------------#
# read command line options
#----------------------------------------#
my $fastafile = "";
my $reffile = "";
my $prefix = 0;
my $program = "";
my $nthreads = 0;
if (scalar @ARGV == 5) {

    ($fastafile, $reffile, $prefix, $program, $nthreads) = @ARGV;

} else {
    Usage();
}
#----------------------------------------#


if ($program eq "mummer-map") {

    for (my $i = 1; $i <= 3; ++$i) {

	print STDERR "### iteration ",  $i, " ###\n";

	if ($i != 1) { my $j = $i - 1; $reffile = "$prefix." . $j . ".newref"; }
	
	my $cmd = "$Bin/mummer-map/mummer-map -o $prefix.$i -p $nthreads $reffile $fastafile";
	print STDERR "Run mummer-map read mapping\n";
	print STDERR "$cmd\n";
	system($cmd);

	
	my $pickref = "breadth";
	if ($i != 1) {$pickref = "all";}
	$cmd = "$Bin/buildcontig/buildcontig -m $prefix.$i.map -r $reffile -o $prefix.$i -c 1 -l 150 -n T -b F -k $pickref";
	print STDERR "# Build contigs\n";
	print STDERR "$cmd\n";
	system($cmd);
	print "\n";
    }
}
else {
    
    for (my $i = 1; $i <= 3; ++$i) {

	print STDERR "### iteration ",  $i, " ###\n";

	if ($i != 1) { my $j = $i - 1; $reffile = "$prefix." . $j . ".newref"; }
	
	my $cmd = "$program/bowtie2-build -q $reffile $prefix.refseq.$i";
	print STDERR "# Build bowtie2 index\n";
	print STDERR "$cmd\n";
	system($cmd);
	
	$cmd = "$program/bowtie2 --sam-nohead --sam-nosq --end-to-end --quiet -k 30 -p $nthreads -x $prefix.refseq.$i -f $fastafile > $prefix.$i.bowtie2";
	print STDERR "# Run bowtie2 read mapping\n";
	print STDERR "$cmd\n";
	system($cmd);

	my $pickref = "breadth";
	if ($i != 1) {$pickref = "all";}
	$cmd = "$Bin/buildcontig/buildcontig -s $prefix.$i.bowtie2 -r $reffile -o $prefix.$i -c 1 -l 150 -n T -b F -k $pickref";
	print STDERR "# Build contigs\n";
	print STDERR "$cmd\n";
	system($cmd);
	print "\n";
    }
}


my $cmd = "";

exit;


sub Usage {
    die("
Usage:
       perl metaCompass.pl <FASTA> <reference> <prefix> <program> <# threads>

Options:
       <FASTA>        DNA reads in FASTA format.
       <reference>    Reference sequences in FASTA format.
       <prefix>       Output file prefix.
       <program>      mummer-map or bowtie2 (specify the DIRECTORY of bowtie2 package).
       <# threads>    # threads to run read mapping.

Output:
       

Contact:
        Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
