#!/usr/bin/perl -w

# This script takes the template pdb and the output of the tcl script pdb to produce a new pdb in which all of the atoms designated in a file called "missing_residues.tmp" or similar.
# It uses a simplistic approach, and replaces the whole residue
# Script will fail if residue numbering in the files is not consecutive
# Script also might give trouble with comments

# USAGE: ./this_file pdb_template pdb_target missing_residues_file

$template	= $ARGV[0];
$target         = $ARGV[1];
$missing_file   = $ARGV[2];

# reading all files
open (TEMPLATE,$template);
open (TARGET,$target);
open (MISSING,$missing_file);

#output pdb file
open (OUT,">${template}_corrected");

#goes over all the residues that need replacing
#sets initial conditions, reading in first line of both files
$template_line = <TEMPLATE>;
$template_residue = 0;

$target_line = <TARGET>;
$target_residue = 0;

while (<MISSING>){
	$current_residue=$_;
	#as long as current residue hasn't been matched, keep printing original file
	while ($current_residue>$template_residue){
		print OUT $template_line;
		# read in next line
		$template_line = <TEMPLATE>;
		$template_line =~ /^ATOM\s*\d+\s*\w+\s*\w+\s*\w+\s*(\d+)/;
		$template_residue = $1;
	} # as soon as there is a match, loop carries on here
	# now it needs to skip to the next residue in the template
	while ($current_residue==$template_residue){
                # read in next line
                $template_line = <TEMPLATE>;
                $template_line =~ /^ATOM\s*\d+\s*\w+\s*\w+\s*\w+\s*(\d+)/;
                $template_residue = $1;
        }
	# now loop runs over target pdb file, looking for the residue that needs to be placed
	while ($current_residue>$target_residue){
		$target_line = <TARGET>;
		$target_line =~ /^ATOM\s*\d+\s*\w+\s*\w+\s*\w+\s*(\d+)/;
		$target_residue = $1;	
	}
	# now it writes down the current residue
        while ($current_residue==$target_residue){
		print OUT $target_line;
                $target_line = <TARGET>;
                $target_line =~ /^ATOM\s*\d+\s*\w+\s*\w+\s*\w+\s*(\d+)/;
                $target_residue = $1;
        }
}

#finish up the rest of the pdb template
print OUT $template_line;
while (<TEMPLATE>){
	print OUT $_;
}

close (OUT);
close (TARGET);
close (MISSING);
close (TEMPLATE);

