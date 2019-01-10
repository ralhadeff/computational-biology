#!/usr/bin/env perl

# input files with xyz coordinates for the center of the residues
$input1	= $ARGV[0];
$input2	= $ARGV[1];
# number of distances to calculate per 'frame' (map)
$length = $ARGV[2];
# distance of the starting state (of the pdb, not the simulation)
$startD = $ARGV[3];
# distance at the ending state (of the pdb, not the simulation)
$endD   = $ARGV[4];

# predict number of frames
$n = `wc $input1 | awk '{print \$1}'`;
$n /= $length ;

# open files
open (ONE,$input1);
open (TWO,$input2);

# counter to keep track of frame
$counter = 0;
# frame number for file naming
$frame = 1;
# reference
$ref = $startD;
# write to first file
open (OUT,">pmf_${frame}.dat");
print "pmf_${frame}.dat $ref 10\n";

while (<ONE>){
	/(\S+)\s+(\S+)\s+(\S+)/;
	$x1 = $1;
	$y1 = $2;
       	$z1 = $3;
	$line = <TWO>;
	$line =~ /(\S+)\s+(\S+)\s+(\S+)/;
	$x2 = $1;
        $y2 = $2;
        $z2 = $3;
        $distance = sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
        if ($counter<$length){
		$counter++;
		printf OUT "%3.3f %3.3f\n", $counter, $distance;
	} else {
		$counter = 0;
		close (OUT);
		$frame++;
		$ref+= ($endD-$startD)/($n-1);
		open (OUT,">pmf_${frame}.dat");
		print "pmf_${frame}.dat $ref 10\n";
	}

}

close (ONE);
close (TWO);
