#!/usr/bin/perl -w

# Script to extract expression values for individual modules.

# 31 August 2010
# Kevin Childs

use strict;
use warnings;
use Getopt::Std;

my $usage = "\n$0  -e expression_matrix  -d directory_with_modules.txt_files   \n\n";

our ($opt_e, $opt_d, $opt_h);
getopts("e:d:h") or die usage();

if ($opt_h) {
    die $usage;
}

my $matrix_file = $opt_e;
my $module_directory = $opt_d;

if (!defined($matrix_file) || 
    !defined($module_directory)) {
    die $usage;
}
if (! -e $matrix_file || 
    ! -d $module_directory) {
    die $usage;
}


my %gene_expression;
open IN, $matrix_file or die "Unable to open $matrix_file";
my $header = <IN>;
while (my $line = <IN>) {
    chomp $line;
    my @elems = split "\t", $line;
    my $gene = shift @elems;
    if (defined($gene)) {
	$gene_expression{$gene} = join "\t", @elems;
    }
}
close IN;

# Get the list of module files.
opendir DIR, $module_directory || die "\nUnable to open $module_directory for reading.\n\n";
my @files= readdir(DIR);
closedir DIR;

foreach my $file (@files) {
    if ($file eq ".." || $file eq ".") {
        next;
    }

    if ($file !~ /module\.txt$/) {
        next;
    }

    $file = $module_directory . "/" . $file;

    $file =~ /([\.\w]+module)\.txt$/;
    my $module_name = $1;
    my $output = $module_directory . "/" . $module_name . "_expression_matrix.txt";

    open OUT, ">$output" or die "Unable to open $output\n";
    print OUT "$header";
    open IN, $file or die "Unable to open $file\n";
    while (my $gene = <IN>) {
	chomp $gene;
	if (exists($gene_expression{$gene})) {
	    print OUT "$gene\t$gene_expression{$gene}\n";
	}
	else {
	    print "Missing expression value for $gene\n";
	}
    }
    close IN;

}
close OUT;

print "\nFinished $0\n\n";
exit;

