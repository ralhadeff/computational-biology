# prints out the resid number of each continuous segment in a pdb
# used to find missing residues, and to compare different structures of the same protein where different parts have been resolved

# read in PDB
$pdb = $ARGV[0];

open(IN,$pdb);

$previous = "-1";

while(<IN>){
	$line  = $_;
#	get type, chain id and resid
	$type  = substr $line, 13, 3;
	$chain = substr $line, 21, 1;
	$ID    = substr $line, 22, 4;
#	only CA 
	if ($type eq "CA "){
#	if they are consecutive, keep going
		if ($ID == $previous+1){
			$previous = $ID;
		} else {
#	if discontinueous, print the start and end of the previous segment
			# trim whitespace
			$previous =~ s/^\s+|\s+$//g ; 
			# skip the first -1 
			if ($previous=="-1"){
				print "$chain $ID - ";
			} else {
				print "$previous\n$chain $ID - " ;
			}
	                $previous = $ID;			
		}
	}
}

$previous =~ s/^\s+|\s+$//g ;
print "$previous\n";


#ATOM      1  N   THR A   9     -21.342  38.912 -27.670  1.00171.15           N

