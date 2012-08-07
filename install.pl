#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);

my $cmd = "";
 
print "\n";
print "# compile libbase for mummer-map\n";
chdir("$Bin/src/mummer-map/libbasedir/");
my $compile = "gcc";
my $opt = "-c -o";
my @libbasepros = ("cleanMUMcand", "clock", "mapfile", "multiseq", "procopt", "safescpy", "seterror", "space");
foreach my $pro (@libbasepros) {
    $cmd = "$compile $opt $pro.o $pro.c";
    print "$cmd\n";
    system($cmd);
}
$cmd = "ar sru libbase.a cleanMUMcand.o clock.o mapfile.o multiseq.o procopt.o safescpy.o seterror.o space.o";
print "$cmd\n";
system($cmd);
chdir($Bin);
print "\n";


print "# compile stree for mummer-map\n";
chdir("$Bin/src/mummer-map/streesrc/");
$opt = "-I../libbasedir -DSTREEHUGE -c -o";
my @streepros = ("construct", "access", "scanpref", "linkloc", "depthtab", "ex2leav", "dfs", "overmax", "oversucc", "addleafcount", "iterator");
foreach my $pro (@streepros) {
    $cmd = "$compile $opt $pro.o $pro.c";
    print "$cmd\n";
    system($cmd);
}
$cmd = "ar sru libstree.a construct.o access.o scanpref.o linkloc.o depthtab.o ex2leav.o dfs.o overmax.o oversucc.o addleafcount.o iterator.o";
print "$cmd\n";
system($cmd);
chdir($Bin);
print "\n";


print "# compile mummer-map\n";
chdir("$Bin/src/mummer-map/mainsrc/");
$opt = "-I../libbasedir -I../streesrc -c -o";
my @pros = ("mummer-map", "maxmatopt", "fasta", "findmaxmat", "procmaxmat");
foreach my $pro (@pros) {
    $cmd = "$compile $opt $pro.o $pro.c";
    print "$cmd\n";
    system($cmd);
}

$cmd = "$compile -lpthread mummer-map.o maxmatopt.o fasta.o findmaxmat.o procmaxmat.o ../streesrc/libstree.a ../libbasedir/libbase.a -o $Bin/bin/mummer-map";
print "$cmd\n";
system($cmd);
chdir($Bin);
print "\n";


$cmd = "./src/metaphyler/installMetaphyler.pl";
print "$cmd\n\n";
system($cmd);

$cmd = "g++ -Wall -W -O2 -o ./bin/extractSeq ./refseq/extractSeq.cpp";
print "$cmd\n\n";
system($cmd);



$cmd = "g++ -Wall -W -O2 -o ./bin/buildcontig src/buildcontig/buildcontig.cpp src/buildcontig/cmdoptions.cpp src/buildcontig/memory.cpp src/buildcontig/procmaps.cpp src/buildcontig/outputfiles.cpp";
print "$cmd\n\n";
system($cmd);

