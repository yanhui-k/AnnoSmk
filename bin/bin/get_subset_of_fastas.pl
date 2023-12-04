#!/usr/bin/perl -w
# download from https://github.com/Childs-Lab/GC_specific_MAKER
# get_subset_of_fastas.pl

# 22 November 2008
# Kevin Childs

# This list reads a list of gene ids from a file.  The script also reads a file of fastas.
# The sequences of the genes that are on the list are output to a third file.

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Std;
use strict;

our ($opt_l, $opt_f, $opt_o, $opt_h);
getopts("f:l:o:h") or die usage();

if ($opt_h) {
    usage();
}
my $fasta_file = $opt_f;
my $list_file = $opt_l;
my $output_file = $opt_o;
if (!defined($fasta_file) && !defined($list_file)) {
    die usage();
}
if (!defined($output_file)) {
    die usage();
}
if (!(-e $fasta_file) || !(-e $list_file)) {
    print "\nThere is a problem with one of the input files.\n";
    die;
}
if (-e $output_file) {
    print "\nThe file, $output_file, already exists.\n";
    die;
}

my %list_ids;

open IN, "$list_file" || die "\nUnable to open $list_file for reading.\n\n";
while (my $line = <IN>) {
    chomp $line;
    $list_ids{$line} = 1;
}

my $fastas = Bio::SeqIO->new(-file => $fasta_file, -format => "fasta");
my $output = Bio::SeqIO->new(-file => ">$output_file", -format => "fasta");
while (my $seq_obj = $fastas->next_seq()) {
    my $id = $seq_obj->id();
    #$id =~ /([^\.]+)/;
    #$id = $1;
    #$id =~ tr/TG/tg/;

    if (exists($list_ids{$id})) {
	$output->write_seq($seq_obj);
    }
}


exit;

sub usage {
    print "\nUsage: $0 -l list_gene_ids  -f multifasta_file  -o output_fasta_file\n\n";

    exit(1);
}
    