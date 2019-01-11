# tcl script for vmd rigid body docking
# written by Raphael Alhadeff Oct 2018 (v.1)
# the purpose is to produce non-overlapping ocnfigurations, moving one protein with respect to the other by small perturbations, then using the MCPT CG energy of molaris (will be submitted externally using a perl script) to run an overarching MC simulation that will refine the docking

# rotate a protein around its own geometric COM
proc rotate_com {selection axis value} {
        set com [measure center $selection]
        $selection moveby [vecscale -1 $com]
        set matrix [transaxis $axis $value]
        $selection move $matrix
        $selection moveby [vecscale 1 $com]
}

# get number of atoms that collide between two selections (given as strings)
proc areWithin {selectionA selectionB ID cutoff} {
        set sel [atomselect $ID "$selectionA and within $cutoff of $selectionB"]
	$sel num
}

# generate configurations
proc runIt {receptor ligand ID distanceCutoff overlapCutoff translation rotation cycles} {
	for {set i 0} {$i<$cycles} {incr i} {
		# generate random translation step, isometric with size [-translation/2,translation/2]
		set trans "[expr $translation*(rand()-0.5)] [expr $translation*(rand()-0.5)] [expr $translation*(rand()-0.5)]"
		# generate random rotation vector, isometric with size [-rotation/2,rotation/2
		set rot "[expr $rotation*(rand()-0.5)] [expr $rotation*(rand()-0.5)] [expr $rotation*(rand()-0.5)]"
		# move and rotate
		rotate_com [atomselect $ID $ligand] x [lindex $rot 0]
		rotate_com [atomselect $ID $ligand] y [lindex $rot 1]
		rotate_com [atomselect $ID $ligand] z [lindex $rot 2]
		[atomselect $ID $ligand] moveby $trans
		# check if new configuration overlaps
		if {[areWithin $receptor $ligand $ID $distanceCutoff ] <= $overlapCutoff} {
			# success, present and save
			display update
			[atomselect $ID "all"] writepdb configuration_$i.pdb
		} else {
			# failed, try again
			set i [expr $i-1]
		}		
		# revert
		[atomselect $ID $ligand] moveby [vecscale -1 $trans]
		rotate_com [atomselect $ID $ligand] z [expr -[lindex $rot 2]]
		rotate_com [atomselect $ID $ligand] y [expr -[lindex $rot 1]]
		rotate_com [atomselect $ID $ligand] x [expr -[lindex $rot 0]]
	}
}

