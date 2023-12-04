#!/usr/bin/perl -w
# download from https://github.com/Childs-Lab/GC_specific_MAKER
# fathom_to_genbank.pl

# Kevin Childs

# Converts the uni.ann and uni.dna output from SNAP's fathom to genbank format.

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeatureI;
use Bio::Location::Split;
use Getopt::Long;
use strict;
use warnings;

my $usage = "\n$0\n    --annotation_file  uni.ann from fathom\n" .
                  "    --dna_file  uni.dna from fathom\n" .
                  "    --genbank_file  output file in genbank format\n" .
                  "    --number  number of annotations to convert\n" .
                  "    [--help]\n\n";

my ($annotation_file, $dna_file, $genbank_file, $convert_num);
my $help;

Getopt::Long::GetOptions( "annotation_file=s" => \$annotation_file,
                          "dna_file=s" => \$dna_file,
                          "genbank_file=s" => \$genbank_file,
                          "number=i" => \$convert_num,
                          "help" => \$help) || die;

if (defined($help)) {
    print $usage;
    exit;
}

if (!defined($annotation_file) ||  !(-e $annotation_file) ||
    !defined($dna_file) ||  !(-e $dna_file) ||
    !defined($genbank_file) ||  (-e $genbank_file) ||
    !defined($convert_num)) {
    die $usage;
}

# Read the uni.dna file and put all sequences into a hash.
my %dna_seqs;
my $fastas = Bio::SeqIO->new(-file => $dna_file, -format => "fasta");
while (my $seq_obj = $fastas->next_seq()) {
    my $new_seq_obj = Bio::Seq->new( -seq => $seq_obj->seq(),
				     -id  => $seq_obj->display_id());
    $dna_seqs{$seq_obj->display_id()} = $new_seq_obj;
}

# Read the uni.ann file and store exon coordinates and strand information in hashes.
# >region-420 (5046649 0..500)
# Einit  	262    	285    	MODEL1025
# Exon   	301    	335    	MODEL1025
# Eterm  	432   	452   	MODEL1025

my (%annot_strands, %annot_coords, $annot_id);
my $annot_fh;
open $annot_fh, $annotation_file or die "\nUnable to open $annotation_file for reading.\n\n";
while (my $line = <$annot_fh>) {

    chomp $line;

    if ($line =~ /^>(\S+)/) {
	$annot_id = $1;
    }
    else {
	my ($exon, $coord1, $coord2, $model) = split '\t', $line;	
	if ($coord1 < $coord2) {
	    $annot_strands{$annot_id} = '+';
	    push @{$annot_coords{$annot_id}}, $coord1 . "\t" . $coord2;
	}
	else {
	    $annot_strands{$annot_id} = '-';
	    push @{$annot_coords{$annot_id}}, $coord2 . "\t" . $coord1;
	}
    }
}
close $annot_fh;


my $num_processed = 0;
my $output = Bio::SeqIO->new(-file => ">$genbank_file", -format => "genbank");
foreach my $seq_id (keys(%dna_seqs)) {
    my $seq_obj = $dna_seqs{$seq_id};

    my $strand;
    if ($annot_strands{$seq_id} eq '+') {
	$strand = 1;
    }
    else {
	$strand = -1;
    }

    my $source_feat =  Bio::SeqFeature::Generic->new( 
	-start        => 1, 
	-end          => $seq_obj->length(),
	-primary      => 'source'
    );
    $seq_obj->add_SeqFeature($source_feat);

    my $split_location = Bio::Location::Split->new();
    foreach my $coord_string (@{$annot_coords{$seq_id}}) {
	my ($start, $end) = split "\t", $coord_string;
	#print "full $coord_string\n";
	#print "start $start\n";
	#print "end $end\n";

	$split_location->add_sub_Location(Bio::Location::Simple->new(-start=> $start,
								     -end=> $end,
								     -strand=> $strand));
    }

    my $source_cds =  Bio::SeqFeature::Generic->new( 
	-strand       => $annot_strands{$seq_id},
	-primary      => 'CDS',
	-location     => $split_location
    );
    $seq_obj->add_SeqFeature($source_cds);

    $output->write_seq($seq_obj);

    ++$num_processed;
    if ($num_processed == $convert_num) {
	last;
    }
}

