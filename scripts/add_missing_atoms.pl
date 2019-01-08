#!/usr/bin/env perl

# This script will do the work that the two other scripts require
# for this script to work, you'll need one pdb of the original coordinates and a second pdb with the missing residues, and they have to be numbered the same!
# The rest is done by itself.
# pay attention that some files might be overwritten if you happen to have files with the same name on your computer

# Usage: perl <this script> original.pdb missing.pdb indent

$original	= $ARGV[0];
$target         = $ARGV[1];

# makes sure that the user is doing everything correctly:

print "\n********************\nPlease make sure that $original and $target have the same number of residues and that they are numbered the same. You can use http://dicsoft2.physics.iisc.ernet.in/pdbgoodies/inputpage.html for renumbering\n*********************\nType anything to continue or no to abort\n";

# last chance to cancel
chomp ( $reply = <STDIN> );
if (($reply eq "n") or ($reply eq "no")) {exit;}

# get the number of the first residue
open (GETFIRSTRESIDUE,$original);
$n=-1;
while($n==-1){
        $_ = <GETFIRSTRESIDUE>;
        /^ATOM\s*\d+\s*\w+\s*\w+\s*\w+\s*(\d+)/ ;
        $n=$1-1;
}
close (GETFIRSTRESIDUE);


# first step, create the missing_residues.tmp file
# create an input file for molaris and run it
# change the molaris path here if it's different
open (MOLINP,">molaris_input_file_for_script.inp");
print MOLINP "$original\nexit\n";
close (MOLINP);

system('molaris molaris_input_file_for_script.inp output.out');
system("cat ./molaris_input_file_for_script/output.out | grep 'is trimmed to GLY' | awk '{print \$3+$n}' > missing_residues_for_script.tmp");
# delete temp rubbish
system("rm -rf ./molaris_input_file_for_script/");
system("rm ./molaris_input_file_for_script.inp");

# create and run the vmd script
open (TCLINP,">tcl_for_script.tcl");
print TCLINP "
set template_path [lindex \$argv 0] 
set target_path [lindex \$argv 1]
set missing_residues_file [lindex \$argv 2]

proc read_missing_residues { file_path template target } {
	set infile [open \$file_path r]
	set cur_residue [gets \$infile]
	while {\$cur_residue > 0 } {
		align_sidechain \$template \$target \$cur_residue
		set cur_residue [gets \$infile]
		if {[eof \$infile]} {
    			close \$infile
    			set cur_residue 0
    			break
		}
                align_sidechain \$template \$target \$cur_residue
	}
}

proc align_sidechain { template target residue } {
	set original [atomselect \$template \"resid \$residue and (backbone or type CB)\"]
	set corrected_core [atomselect \$target \"resid \$residue and (backbone or type CB)\"]
	set corrected [atomselect \$target \"resid \$residue\"]
	set M [measure fit \$corrected_core \$original]
	\$corrected move \$M
}


        set template [mol load pdb \$template_path]
        set target [mol load pdb \$target_path]
        read_missing_residues \$missing_residues_file \$template \$target
        set output_pdb [atomselect \$target \"all\"]       
        \$output_pdb writepdb template_pdb_for_script.pdb
	quit
";
close (TCLINP);

system("vmd -dispdev none -e tcl_for_script.tcl -args $original $target ./missing_residues_for_script.tmp");
# delete the script
system("rm ./tcl_for_script.tcl");

# create the new pdb

open (TEMPLATE,$original);
open (TARGET,"template_pdb_for_script.pdb");
open (MISSING,"./missing_residues_for_script.tmp");

#output pdb file
open (OUT,">new_corrected_$original");

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
        while (($current_residue==$target_residue) && (defined $target_line) ){
		print OUT $target_line;
		print $target_line;
                $target_line = <TARGET>;
                $target_line =~ /^ATOM\s*\d+\s*\w+\s*\w+\s*\w+\s*(\d+)/;
                $target_residue = $1;
        }
}

#finish up the rest of the pdb template
print OUT $template_line;
print $template_line;
while (<TEMPLATE>){
	print OUT $_;
	print $_;
}

close (OUT);
close (TARGET);
close (MISSING);
close (TEMPLATE);

# delete the rest of the rubbish
system("rm ./missing_residues_for_script.tmp");
system("rm ./template_pdb_for_script.pdb");

