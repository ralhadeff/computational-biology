proc build_dihedral_vector_pdb {molA molB alpha} {
# this rotates all dihedrals of the FIRST molID provides by the delta of the dihdrals (second-first) times an alpha factor
# both pdbs MUST be of the exact same atom number and order

	set numerOfResidues [[atomselect $molA "type CA"] num]
	set firstResidue [[atomselect $molA "index 0"] get resid]

	set dihedrals {}

	for {set i [expr $firstResidue + 1]} {$i < [expr $firstResidue+$numerOfResidues -1]} {set i [expr $i + 1]} {
	# residues loop starts from 2 to avoid null for the i-1 C atom
	# TODO residue 1 is actually skipped, consider adding manual treatment above. Not important probably
	# this loops just puts the dihedrals in the list
		set j [expr $i - 1]
		set k [expr $i + 1]
		# start with CA-C-N-CA
		mol top $molA
		set dihA [measure dihed [[atomselect $molA "(resid $j and type CA C) or (resid $i and type N CA)"] get index]]
		mol top $molB
		set dihB [measure dihed [[atomselect $molB "(resid $j and type CA C) or (resid $i and type N CA)"] get index]]
		set diff [expr ($dihB - $dihA)]
		if {$diff > 180} {
			set diff [expr $diff - 360]
		}
		if {$diff < -180} {
			set diff [expr $diff + 360]
		}
		set diff [expr $diff * $alpha]
		lappend dihedrals $diff
                puts "$j $i $dihA $dihB $diff"
		# then C-N-CA-C
                mol top $molA
                set dihA [measure dihed [[atomselect $molA "(resid $j and type C) or (resid $i and type N CA C)"] get index]]
                mol top $molB
                set dihB [measure dihed [[atomselect $molB "(resid $j and type C) or (resid $i and type N CA C)"] get index]]
                set diff [expr ($dihB - $dihA)]
                if {$diff > 180} {
                        set diff [expr $diff - 360]
                }
                if {$diff < -180} {
                        set diff [expr $diff + 360]
                }
                set diff [expr $diff * $alpha]
               lappend dihedrals $diff
                puts "$j $i $dihA $dihB $diff"
		# finally N-CA-C-O
                mol top $molA
                set dihA [measure dihed [[atomselect $molA "(resid $i and type N CA C) or (resid $k and type N)"] get index]]
                mol top $molB
                set dihB [measure dihed [[atomselect $molB "(resid $i and type N CA C) or (resid $k and type N)"] get index]]
                set diff [expr ($dihB - $dihA)]
                if {$diff > 180} {
                        set diff [expr $diff - 360]
                }
                if {$diff < -180} {
                        set diff [expr $diff + 360]
                }
                set diff [expr $diff * $alpha]
                lappend dihedrals $diff
                puts "$j $i $dihA $dihB $diff"
	}

	set count 0
	mol top $molA
	for {set i [expr $firstResidue + 1]} {$i < [expr $firstResidue+$numerOfResidues-1]} {set i [expr $i + 1]} {
	# this loops does the rotations
		set j [expr $i - 1]
		# start with CA-C-N-CA
		set indexS [[atomselect $molA "resid $j and type C"] get index]
		set indexE [[atomselect $molA "resid $i and type N"] get index] 
		rotate_around_bond $molA $indexS $indexE [lindex $dihedrals $count]
		set count [expr $count +1]
		# then C-N-CA-C
                set indexS [[atomselect $molA "resid $i and type N"] get index]
                set indexE [[atomselect $molA "resid $i and type CA"] get index]   
                rotate_around_bond $molA $indexS $indexE [lindex $dihedrals $count]
                set count [expr $count +1]
		# finally N-CA-C-O
                set indexS [[atomselect $molA "resid $i and type CA"] get index]
                set indexE [[atomselect $molA "resid $i and type C"] get index]   
		puts "$i $indexS $indexE"
                rotate_around_bond $molA $indexS $indexE [lindex $dihedrals $count]
                set count [expr $count +1]
	}
}


proc rotate_around_bond {moleculeID atomA atomB degrees} {
# rotate all residues AFTER an index around the selected bond 

	set a1 [atomselect $moleculeID "index $atomA"]
        set a2 [atomselect $moleculeID "index $atomB"]

	set c1 [lindex [$a1 get {x y z}] 0]
	set c2 [lindex [$a2 get {x y z}] 0]

	set sel [atomselect $moleculeID "index > $atomB"]

	set rot_mat [trans bond $c1 $c2 $degrees deg]

	$sel move $rot_mat
}
