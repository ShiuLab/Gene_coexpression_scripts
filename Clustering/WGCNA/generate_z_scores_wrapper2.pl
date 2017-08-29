#! /usr/bin/perl -w

# A simple wrapper that will run generate_z_scores_for_matrix.pl on all 
# appropriate module files in a directory.

# Kevin Childs
# 14 September 2010

use strict;
use warnings;
use Getopt::Std;

my $usage = "\n$0   -d directory_with_files\n\n";

our ($opt_d,  $opt_h);
getopts("d:h") or die usage();

if ($opt_h) {
    die $usage;
}

my $input_directory = $opt_d;

if (!defined($input_directory) || 
    !(-d $input_directory)) {
    die $usage;
}

opendir DIR, $input_directory || die "\nUnable to open $input_directory for reading.\n\n";
my @files = readdir(DIR);
closedir DIR;
foreach my $file (@files) {
    if ($file eq ".." || $file eq ".") {
        next;
    }

    my $base;
    if ($file =~ /([\w\.]+)_module_expression_matrix\.txt$/) {
        $base = $1;
    }
    else {
	next;
    }
    my $full_file = $input_directory . "" . $file;
    my $results_file = $input_directory . "" . $base . "_module_expression_matrix_z_scores.txt";

    my $z_call = "/mnt/home/uygunsah/WGCNA/codes/generate_z_scores_wrapper.pl  -i  $full_file  -o  $results_file";
    print "$z_call\n";
    system($z_call);
    
}

print "\n\nFinished: $0\n\n";
exit;





