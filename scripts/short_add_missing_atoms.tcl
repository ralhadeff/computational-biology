# This script is supposed to help you add missing atoms to a pdb file.
# It assumes that all backbone and CB are present. Otherwise, it will crash.
# To override this, manually delete from the missing_residues.tmp file any residue you wish to skip.

# Using MaxSprout or a similar server, create a pdb with all the residues complete. This is the template file, then run it through PDBGoogies to renumber the residues.

# In terminal, on a molaris output file (where you loaded the original pdb), type:
# > cat output.out  | grep "is trimmed to GLY" | awk '{print $3+X}' > missing_residues.tmp
# where X is the value to correct the renumbering of molaris

set template_path [lindex $argv 0] 
set target_path [lindex $argv 1]
set missing_residues_file [lindex $argv 2]

proc read_missing_residues { file_path template target } {
	set infile [open $file_path r]
	set cur_residue [gets $infile]
	while {$cur_residue > 0 } {
		align_sidechain $template $target $cur_residue
		puts "moved resid $cur_residue"
		set cur_residue [gets $infile]
		if {[eof $infile]} {
    			close $infile
    			set cur_residue 0
    			break
		}
                align_sidechain $template $target $cur_residue
	}
}

proc align_sidechain { template target residue } {
	set original [atomselect $template "resid $residue and (backbone or type CB)"]
	set corrected_core [atomselect $target "resid $residue and (backbone or type CB)"]
	set corrected [atomselect $target "resid $residue"]
	set M [measure fit $corrected_core $original]
	$corrected move $M
}


        set template [mol load pdb $template_path]
        set target [mol load pdb $target_path]
        read_missing_residues $missing_residues_file $template $target
        puts "test"
        set output_pdb [atomselect $target "all"]       
        $output_pdb writepdb ${target_path}_corrected
	quit
