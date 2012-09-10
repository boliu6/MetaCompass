#!/usr/bin/perl

#############################################
#
# Program: Metagenomic comparative assembly
#
# Author: Bo Liu, boliu@umiacs.umd.edu
#
# Sun Aug 12 16:12:16 EDT 2012
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
my $niter = 3;
if (scalar @ARGV == 6) {

    ($fastafile, $reffile, $prefix, $program, $niter, $nthreads) = @ARGV;

} else {
    Usage();
}
#----------------------------------------#


if ($program eq "mummer-map") {

    for (my $i = 1; $i <= $niter; ++$i) {

	print STDERR "### iteration ",  $i, " ###\n";

	if ($i != 1) { my $j = $i - 1; $reffile = "$prefix." . $j . ".newref"; }
	
	my $cmd = "$Bin/mummer-map -o $prefix.$i -p $nthreads $reffile $fastafile";
	print STDERR "Run mummer-map read mapping\n";
	print STDERR "$cmd\n";
	system($cmd);

	
	my $pickref = "breadth";
	if ($i != 1) {$pickref = "all";}
	$cmd = "$Bin/buildcontig -m $prefix.$i.map -r $reffile -o $prefix.$i -c 1 -l 150 -n T -b F -k $pickref";
	print STDERR "# Build contigs\n";
	print STDERR "$cmd\n";
	system($cmd);
	print "\n";
    }
}
else {
    
    for (my $i = 1; $i <= $niter; ++$i) {

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
	$cmd = "$Bin/buildcontig -s $prefix.$i.bowtie2 -r $reffile -o $prefix.$i -c 1 -l 150 -n T -b F -k $pickref";
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
       perl compass.pl <FASTA> <reference> <prefix> <mapping> <# iterations> <# threads>

Options:
       <FASTA>        DNA reads in FASTA format.
       <reference>    Reference sequences in FASTA format.
       <prefix>       Output file prefix.
       <mapping>      mummer-map or bowtie2 (specify the DIRECTORY of bowtie2 package).
       <# iterations>  # iterations to run comparative assembly. 3 is recommended.
       <# threads>    # threads to run read mapping.

Output:
       

Contact:
        Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
