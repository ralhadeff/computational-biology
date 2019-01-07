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

puts "at this point, dimensions doesn't work and it only does xy relaxation."
puts "also, this only find the local minima"

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
	set test [list [lindex $vector 0] [lindex $vector 1] 0]
	$selection moveby $test
	after 10

set surroundings [atomselect top "(within $threshold of index $indices) and (not index $indices)"]

set center [veczero]
foreach atom [$selection get {x y z}] {set center [vecadd $center $atom]}
set center [vecscale [expr 1.0 /[$selection num]] $center]

display update

}

}

