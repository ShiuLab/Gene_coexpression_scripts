#! /usr/bin/perl 

# 6 September 2010

# This script will use the GSEXXXX_mas5_normalized_averaged_per_gene.txt files.  
# The data in these files are gene-based and have only been averaged across probes.  
# The data are already log2 transformed.  The COV filter script will untransform 
# the data to determine expression levels.  The data will not be filtered based on 
# the number of samples that have no expression levels.  Only genes that have zero 
# expression in all treatments will be discarded.  The output of this script will 
# be an abbreviated matrix of genes and their expression values identical in format 
# to the input file.

# Kevin Childs
#########
# 7 September 2010

# Modify script to output expression data that is averaged across replicates.

# Kevin Childs
#########
# 13 March 2011

# This version of the script does not perform average across replicates before performing
# the CV test.
# The existing headers will be used.

# Kevin Childs

use strict;
use warnings;
use Statistics::Descriptive;
use Getopt::Std;

my $usage = "\n$0  -i input_file  -o output_file  -c cov_low_end_cutoff \n\n";

our ($opt_o, $opt_i, $opt_c, $opt_h);
getopts("o:i:c:h") or die usage();

if ($opt_h) {
    die $usage;
}

my $output_file = $opt_o;
my $input_file = $opt_i;
my $low_cov_cutoff = $opt_c;

if (!defined($output_file) || -e $output_file || !defined($input_file) || !-e $input_file ||
    !defined($low_cov_cutoff)) {
    die $usage;
}

open IN, $input_file or die "Unable to open $input_file\n";
open OUT, ">$output_file" or die "Unable to open $output_file for writing.\n";

my $line = <IN>;
print OUT "$line";  # Print out the header.

while (my $line = <IN>) {
    chomp $line;
    my @elems = split "\t", $line;

    # Check that at least one array value is non-zero.
    my $all_zeroes = 1;
    for (my $i = 1; $i < scalar(@elems); ++$i) {
	if ($elems[$i] != 0) {
	    $all_zeroes = 0;
	}
    }
    if ($all_zeroes == 1) {
	# Skip this gene.
	next;
    }

    # Grab the expression data, ignore the gene name.
    my @this_data = @elems[1 .. (scalar(@elems) - 1)];

    my $stat = Statistics::Descriptive::Sparse->new();
    $stat->add_data(@this_data);
    my $mean = $stat->mean();
    my $std_dev = $stat->standard_deviation();
    my $cov = $std_dev / $mean;

    if ($cov >= $low_cov_cutoff && $cov <= 1000) {
	print OUT "$line\n";
    }
	

}
close IN;
close OUT;

print "\nFinished $0\n\n";
exit;

sub log2 {
    my $n = shift @_;
    return log($n)/log(2);
}


