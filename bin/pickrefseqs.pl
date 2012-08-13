#!/usr/bin/perl

#############################################
#
# Program: Pick reference genomes using MetaPhyler.
#
# Author: Bo Liu, boliu@umiacs.umd.edu
#
# Tue Aug  7 21:19:05 EDT 2012
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
$cmd = "perl $Bin/../src/metaphyler/metaphyler.pl $fastafile $prefix $nthreads";
print STDERR "$cmd\n";
system($cmd);
print STDERR "\n";

# select reference genomes
print STDERR "# Pick reference genomes based on MetaPhyler output\n";
$cmd = "perl $Bin/pickrefids.pl $prefix.blastn 1 > $prefix.refseq.ids";
print STDERR "$cmd\n";
system($cmd);
print STDERR "\n";


print STDERR "# Extract reference genome sequences\n";
$cmd = "$Bin/extractSeq $Bin/../refseq/bacgeno.fna $prefix.refseq.ids > $prefix.refseq.fna";
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
