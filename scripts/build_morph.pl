# This script takes two pdb files, and creates a third one where all coordinates are the result of applying the delta vector of the two, on one of them, n times. Notice that the two pdbs must be IDENTICAL in terms of atom order , number and type etc.

# Usage: perl carteian_vector_continue.pl A.pdb B.pdb n
# notice that n doesn't have to be an integer. 1 n = 1 B-A vector, it's not normalized to any unit.

$pdbA	= $ARGV[0];
$pdbB   = $ARGV[1];
$n	= $ARGV[2];

# makes sure that the user is doing everything correctly:

print "\n********************\nPlease make sure that $pdbA and $pdbB have the same number of atoms and that they are aligned\n*********************\nType anything to continue or no to abort\n";

# last chance to cancel
#chomp ( $reply = <STDIN> );
#if (($reply eq "n") or ($reply eq "no")) {exit;}

# open both pdbs
open (APDB,$pdbA);
open (BPDB,$pdbB);

# create output file name
$pdbA =~ /(.*)\.pdb/;
$name = $1;
$pdbB =~ /(.*)\.pdb/;
$name = $name."_to_".$1."_".$n."_times.pdb";



# create the new pdb
open (RESULT,">$name");

# for each atom, prints out the new coordinates and otherwise same as A
while (<APDB>){
		# get A's part
		# notice that the second \w is with a * because the previous one and the current might be joint
		# the one with ? is for chain, which might not be present
		$lineA = $_;
                $lineB = <BPDB>;
		if ($lineA =~ /(^ATOM\s+\d+\s+\w+\s*\w*\s+\w?\s+\d+)\s+(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)/){
			$xA	=	$2;
			$yA     =       $3;
			$zA     =       $4;
			print RESULT $1;
			# get B's part
			# not using $1 but kept it there for consistency
	                $lineB =~ /(^ATOM\s+\d+\s+\w+\s*\w*\s+\w?\s+\d+)\s+(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)/;
	                $xB     =       $2;
	                $yB     =       $3;
	                $zB     =       $4;
			# calculate coordinates
			$x = $xA + ($xB - $xA)*$n;
	                $y = $yA + ($yB - $yA)*$n;
	                $z = $zA + ($zB - $zA)*$n;
			# print output
			printf RESULT ("    %8.3f%8.3f%8.3f\n", $x, $y, $z);
		} else {
			print RESULT $lineA;
		}
}

close (APDB);
close (BPDB);
close (RESULT);
