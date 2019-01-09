# a bunch of tcl scripts for manipulating and analyzing proteins in VMD

# Move a selection so that the geometric center is at 0 0 0 
proc move_to_origin {selection} {

	set gc [veczero] ;# make a vector
	set sel $selection ;# make a selection
	foreach coord [$sel get {x y z}] {set gc [vecadd $gc $coord]} ;# add all coordinates to vector
	set gc [vecscale [expr 1.0 /[$sel num]] $gc] ;# normalize vector
	$sel moveby [vecscale -1 $gc] ;# move to origin
}

proc align_to_vector {selection vecX vecY vecZ} {
;# aligns the given vector on the molecule to the Z axis

	set gc "$vecX $vecY $vecZ" ;# asign vector
	set sel $selection ;# make a selection
	set M [transvecinv $gc] ;# a transformation
	$sel move $M ;# move
	set M [transaxis y -90] ;# correction to the transformation 
	$sel move $M ;# move again
}


# This script is supposed to help you add missing atoms to a pdb file.
# It assumes that all backbone and CB are present. Otherwise, it will crash.
# To override this, manually delete from the missing_residues.tmp file any residue you wish to skip.

# Using MaxSprout or a similar server, create a pdb with all the residues complete. This is the template file, then run it through PDBGoogies to renumber the residues.

# In terminal, on a molaris output file (where you loaded the original pdb), type:
# > cat output.out  | grep "is trimmed to GLY" | awk '{print $3+X}' > missing_residues.tmp
# where X is the value to correct the renumbering of molaris

proc correct_my_pdb {template_path target_path missing_residues_file } {
	set template [mol load pdb $template_path]
	set target [mol load pdb $target_path]
	read_missing_residues $missing_residues_file $template $target
	set output_pdb [atomselect $target "all"]	
	$output_pdb writepdb ${target_path}_aligned

}

proc correct_my_pdb_bb {template_path target_path missing_residues_file } {
        set template [mol load pdb $template_path]
        set target [mol load pdb $target_path]
        read_missing_residues_bb_only $missing_residues_file $template $target
        set output_pdb [atomselect $target "all"]
        $output_pdb writepdb ${target_path}_aligned

}

proc read_missing_residues_bb_only { file_path template target } {
        set infile [open $file_path r]
        set cur_residue [gets $infile]
        while {$cur_residue > 0 } {
                align_backbone $template $target $cur_residue
                echo "moved resid $cur_residue"
                set cur_residue [gets $infile]
                if {[eof $infile]} {
                        close $infile
                        set cur_residue 0
                        break
                }
                align_backbone $template $target $cur_residue
        }
}


proc read_missing_residues { file_path template target } {
	set infile [open $file_path r]
	set cur_residue [gets $infile]
	while {$cur_residue > 0 } {
		align_sidechain $template $target $cur_residue
		echo "moved resid $cur_residue"
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

proc align_backbone { template target residue } {
        set original [atomselect $template "resid $residue and (backbone)"]
        set corrected_core [atomselect $target "resid $residue and (backbone)"]
        set corrected [atomselect $target "resid $residue"]
        set M [measure fit $corrected_core $original]
        $corrected move $M
}

# end of add missing chunk

# this rotates a selection around its own COM, specifying x y or z axes, value in degrees
proc rotate_com {selection axis value} {
	set com [measure center $selection]
	$selection moveby [vecscale -1 $com]
	set matrix [transaxis $axis $value]
	$selection move $matrix
        $selection moveby [vecscale 1 $com]
}

proc areWithin {selectionA selectionB cutoff} {
	set sel [atomselect top "$selectionA and within $cutoff of $selectionB"]
	$sel num
	# [atomselect top "chain B"] moveby {1 1 1}
}

# the following scripts have been copied from http://www.ks.uiuc.edu/Research/vmd/vmd-1.7/ug/node181.html
proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}

proc geom_center {selection} {
        # set the geometrical center to 0
        set gc [veczero]
        # [$selection get {x y z}] returns a list of {x y z} 
        #    values (one per atoms) so get each term one by one
        foreach coord [$selection get {x y z}] {
           # sum up the coordinates
           set gc [vecadd $gc $coord]
        }
        # and scale by the inverse of the number of atoms
        return [vecscale [expr 1.0 /[$selection num]] $gc]
}

proc gyr_radius {sel} {
  # make sure this is a proper selection and has atoms
  if {[$sel num] <= 0} {
    error "gyr_radius: must have at least one atom in selection"
  }
  # gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
  set com [center_of_mass $sel]
  set sum 0
  foreach coord [$sel get {x y z}] {
    set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
  }
  return [expr sqrt($sum / ([$sel num] + 0.0))]
}
# end of copied scripts

# calculates largest radisu
proc largest_radius {sel} {

 set com [geom_center $sel]
 set largest 0
 foreach coord [$sel get {x y z}] {
	set distance [veclength [vecsub $coord $com]]
	if {$distance > $largest} {
		set largest $distance
	}
 }
 return $largest
}

# rotates a molecule to the current viewpoint
proc rotate_to_this_viewpoint {molecule_index} {

	[atomselect $molecule_index "all"] move [lindex [molinfo $molecule_index get rotate_matrix] 0 ]

}

# load vesso outputs (loads several files with the same name, but non consecutive numbers)
# this is useful because if there are too many files the bash can't handle it
proc load_vesso_files {} {
	for {set i 0} {$i<50000} {incr i} {
		set name [format "%012i_proc.pdb" $i]
		if {[file exists $name]} {animate read pdb $name}
	}
}

# rotates a selection to the current viewpoint
proc rotate_to_this_viewpoint_sel {selection molecule_index} {

        $selection move [lindex [molinfo $molecule_index get rotate_matrix] 0 ]

}

# This is suppose to move an atom around until it is at least threshold from anything else.
proc relax_atom {selection threshold step dimensions} {

# start by QA
if {[$selection num] == 0} {
	puts "no selection"
	return
}
if {$threshold < $step} {
        puts "steps must be smaller than the threshold"
        return
}

puts "This only find the local minima"

if {$dimensions=="xy"} {
	puts "will only perform XY relaxation (z constant)"
}


# get selection index
set indices [$selection get index]

# get selection center
set center [veczero] 
foreach atom [$selection get {x y z}] {set center [vecadd $center $atom]}
set center [vecscale [expr 1.0 /[$selection num]] $center] 

# get surrounding within threshold
set surroundings [atomselect top "(within $threshold of index $indices) and (not index $indices)"]

# at this point, it could get trapped in an endless loop
while {[$surroundings num] != 0} {

	#first, get the closest atom
	set shell [atomselect top "(within $step of index $indices) and (not index $indices)"]

	for {set i $step} {[$shell num]==0} {set i [expr $i + $step]} {
		set shell [atomselect top "(within $i of index $indices) and (not index $indices)"]
	}

	puts $i
	puts [$shell get index]

	set vector [veczero]
	foreach coordinate [$shell get {x y z}] {
		set vector [vecadd $vector $center]
		set vector [vecsub $vector $coordinate]
	}

	set vector [vecscale [expr $step/[veclength $vector]] $vector]
	puts $vector
	if {$dimensions=="xy"} {
		set test [list [lindex $vector 0] [lindex $vector 1] 0]
	} else {
		set test [list [lindex $vector 0] [lindex $vector 1] [lindex $vector 2]]
	}
	$selection moveby $test
	after 10

set surroundings [atomselect top "(within $threshold of index $indices) and (not index $indices)"]

set center [veczero]
foreach atom [$selection get {x y z}] {set center [vecadd $center $atom]}
set center [vecscale [expr 1.0 /[$selection num]] $center]

display update

}

}
