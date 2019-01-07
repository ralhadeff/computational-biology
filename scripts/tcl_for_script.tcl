
set template_path [lindex $argv 0] 
set target_path [lindex $argv 1]
set missing_residues_file [lindex $argv 2]

proc read_missing_residues { file_path template target } {
	set infile [open $file_path r]
	set cur_residue [gets $infile]
	while {$cur_residue > 0 } {
		align_sidechain $template $target $cur_residue
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
        set output_pdb [atomselect $target "all"]       
        $output_pdb writepdb template_pdb_for_script.pdb
	quit
