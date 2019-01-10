#!/usr/bin/env perl

use Switch;

# This script takes a txt file that I created using grep to analyze the big data from the large pdld analysis
# it is a modification that does all of the trimming as well.
$txt	= $ARGV[0];

# checks that file exists
if (!(-e $txt)){
	print "\nfile does not exist!\n";
	exit;
}

# read in file
open (IN,"$txt");
# find all labels that are present in file
# counter for number of labels
$k = 0;
# get the size of the dataset
$minMax = 0;
while (<IN>){
        /..\/frame_(\-?\d+)\/conf_(\d+)\/(\w+)\/bind_pdld_s_lra_beta.out: (\S+)/;
	if (abs($1) > $minMax){
		$minMax = abs($1);
	}
	$found = -1;
	for ($i=0;$i<$k;$i++){
		if ($3 eq $labels[$i]){
			# labels already marked
			$found = 1;
		}
	}
	if ($found == -1){
		# add labels, update counter
		$labels[$k] = $3;
		$k++;
	}	
}
close (IN);

# create the data matrix with all numbers set to 0. Matrix is data[frame(-15 to 15)][protonation state][value]
# because indices have to be positive, I'm gonna add +20 to the number for the indices
# also, I'm gonna have 154 = 0, 155 = 1, both = 2, none = 3
# and finally, value = 0 is the sum, and value = 1 is the count. value = 2 will be the average and value = 3 will be the variance
for ($i=-$minMax;$i<=$minMax;$i++){
	for ($j=0;$j<=$k;$j++){
		# sum
		$data[$i+$minMax][$j][0]=0;
	        # count
		$data[$i+$minMax][$j][1]=0;
		# average
                $data[$i+$minMax][$j][2]=0;
		# SS / variance
                $data[$i+$minMax][$j][3]=0;
	}
}

open (IN,"$txt");
# first run, computes the total and count so that the average can be created.
while(<IN>){
	/..\/frame_(\-?\d+)\/conf_(\d+)\/(\w+)\/bind_pdld_s_lra_beta.out: (\S+)/;
	# translate protonation
        $found = 0;
        for ($i=0;$i<$k;$i++){
		if ($3 eq $labels[$i]){
			$data[$1+$minMax][$i][0]+=$4;
			$data[$1+$minMax][$i][1]++;
			$found = 1;
			last;
		}
	}
	if ($found==0){
		$data[$1+$minMax][$k][0]+=$4;
                $data[$1+$minMax][$k][1]++;

	}
}

close(IN);

# print out results so far and calculate averages
for ($i=-$minMax;$i<=$minMax;$i++){
        for ($j=0;$j<=$k;$j++){
		if ($data[$i+$minMax][$j][1]!=0){
			$data[$i+$minMax][$j][2] = ($data[$i+$minMax][$j][0]+0.0) / $data[$i+$minMax][$j][1];
		} else {
			$data[$i+$minMax][$j][2] = 99;
		}
	}
}

# run over file again and calculate variance (SoS)
open (IN,"$txt");

while(<IN>){
        /..\/frame_(\-?\d+)\/conf_(\d+)\/(\w+)\/bind_pdld_s_lra_beta.out: (\S+)/;
        # translate protonation
        $found = 0;
        for ($i=0;$i<$k;$i++){
                if ($3 eq $labels[$i]){
			$data[$1+$minMax][$i][3]+=( $data[$1+$minMax][$i][2]- $4)**2;
                        $found = 1;
                        last;
                }
        }
        if ($found==0){
		$data[$1+$minMax][$k][3]+=( $data[$1+$minMax][$k][2]- $4)**2;
        }
}

printf "frame\tprot\tcount\tsum\taverage\tvariance\n";

# calculate variance and print out results
for ($j=0;$j<=$k;$j++){
	for ($i=-$minMax;$i<=$minMax;$i++){
		if ($data[$i+$minMax][$j][1]!=0){
	                $data[$i+$minMax][$j][3] = (($data[$i+$minMax][$j][3]+0.0) / $data[$i+$minMax][$j][1])**(0.5);
		} else {
			$data[$i+$minMax][$j][3] = 99;
		}
                if (!($data[$i+$minMax][$j][2]==99 && $data[$i+$minMax][$j][3]==99)){
			printf "$i\t%s\t$data[$i+$minMax][$j][1]\t$data[$i+$minMax][$j][0]\t%.3f\t%.3f\n", $labels[$j], $data[$i+$minMax][$j][2], $data[$i+$minMax][$j][3] ;
		}
        }
}

close(IN);
