#!/usr/bin/perl

#############################################
#
# Program: Given a set of metagenomic sequences,
#          pick reference genomes for
#          comparative assembly.
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
my $prefix = 0;
my $nthreads = 0;
if (scalar @ARGV == 3) {
    ($fastafile, $prefix, $nthreads) = @ARGV;
} else {
    Usage();
}
#----------------------------------------#


my $cmd = "";

# run metaphyler
print STDERR "# Run MetaPhyler\n";
$cmd = "perl $Bin/metaphyler/runMetaphyler.pl $fastafile $prefix $nthreads";
print STDERR "$cmd\n";
system($cmd);
print STDERR "\n";

# select reference genomes
print STDERR "Pick reference genomes based on MetaPhyler output\n";
$cmd = "perl $Bin/refseq/pickrefids.pl $prefix.blastn 1 > $prefix.refseq.ids";
print STDERR "$cmd\n";
system($cmd);
print STDERR "\n";


print STDERR "Extract reference genome sequences\n";
$cmd = "$Bin/refseq/extractSeq $Bin/refseq/bacgeno.fna $prefix.refseq.ids > $prefix.refseq.fna";
print STDERR "$cmd\n";
system($cmd);
print STDERR "\n";


exit;


sub Usage {
    die("
Usage:
       perl pickrefseqs.pl <FASTA> <output prefix> <# threads>

Options:
       <FASTA>        DNA reads in FASTA format.
       <prefix>       Output file prefix.
       <# threads>    Used to run BLAST.

Output:
       MetaPhyler output:
             prefix.blast[n/x]
                   Raw blast output.
             prefix.classification
                   Classification results.
             prefix.<genus|family|order|class|phylum>.taxprof
                   Taxonomy profiles at each level.

       Reference genomes:
             prefix.refseq.fna

Description:
       This program first runs MetaPhyler to estimate the taxonomic composition.
       Then extract corresponding genomes for comparative assembly.

Contact:
        Have problems? Contact Bo Liu - boliu\@umiacs.umd.edu

");
}
