#!/usr/bin/env perl

# Raphael Alhadeff June 2016
# this script is made to rescale an xy data file so that the x values are on a different scale
# the purpose is to allow subtraction of two datasets that originally have different sizes and are not aligned in the x-axis
# it will dumbly run over min to max numbers at interval intervals, simply averaging over all the values it could grab from the data set provided. If there are no values in that range it will interpolate using linear regression

# input dataset needs to be a simple list of 2 numbers, x and y
# input command line will be min, max, interval and it will print out the dataset (so user can copy paste or > into a file

# the code makes no quality control, USER HAS TO OPEN BEFORE AND AFTER AND VISUALLY CHECK IT LOOKS GOOD
# also user has to make sure that the scale makes sense (i.e. in the range of the original x values)

#

# Usage: perl <this script> <dataset file> <min> <max> <interval>

$input		= $ARGV[0];
$min		= $ARGV[1];
$max		= $ARGV[2];
$interval 	= $ARGV[3];

if (!(-e $input)){
	print "\ninput file does not exist!\n";
	exit;
}


# open file, but sort it first
@rawData = `sort -n $input`;
# split data into x and y

for ($i = 0;$i<scalar @rawData; $i++){
	$rawData[$i]=~/(\S+)\s+(\S+)/;
	$dataX[$i] = $1;
	$dataY[$i] = $2;
}

# x runs on the original scale
$x = 0;
$oldX = $dataX[$x];

# j runs on the new scale
$j = 0;

# bin the original data
for ($newX=$min;$newX<=$max;$newX+=$interval){
	# grab all values that are below this 'bin'
	$sum = 0;
	$count = 0;
	while ($oldX<=$newX && $x<scalar @dataX){
		$count++;
		# add up y's and increment i
		$sum+=$dataY[$x];
		$x++;
		# update x
		$oldX = $dataX[$x];		
	}
	# print out
	if ($count!=0){
		$newY = $sum / $count;
		$newX[$j] = $newX;
		$newY[$j] = $newY;
		$j++;
	} else {
		$newY = "A";
                $newX[$j] = $newX;
                $newY[$j] = $newY;
                $j++;
	}
}

# change all A's into average values

# start with all leading A's, arbitrarily set them as the first value
$first = 0;
$y = $newY[$first];
while ($y == "A"){
	$y = $newY[$first];
	$first++;
}
# correction
$first--; 
# update
for ($i=0;$i<$first;$i++){
	$newY[$i]=$newY[$first];
}

# proceed with all final A's, arbitrarily set them as the last
$last = scalar @newX;
$y = $newY[$last];
while ($y == "A"){
        $y = $newY[$last];
        $last--;
}
# correction
$last++;
# update
for ($i=scalar @newX;$i>$last;$i--){
        $newY[$i]=$newY[$last];
}

# end with in betweens
for ($i=$first;$i<=$last;$i++){
	# as long as Y is not A, keep record of the value
	if ($newY[$i]=="A"){
		# count how many A's and keep record of the next value that is not A
		for ($j=$i+1;$j<=$last;$j++){
			if ($newY[$j]!="A"){
				# update all until this point
				# dY
				$gap = $newY[$j]-$newY[$i-1];
				# dY interval
				$gap/= $j-($i-1);
				for ($k=$i;$k<$j;$k++){
					$newY[$k]=$newY[$i-1]+($k-($i-1))*$gap;
				}
				# return to i loop and exit j loop
				$i = $j-1;
				$j = $last;
			}
		}
	}
}

# print output
for ($i=0;$i<scalar @newX ; $i++){
        if ($newY[$i]!="A"){
                printf "%7.4f\t%7.4f\n", $newX[$i],$newY[$i];
        } else {
                printf "%7.4f\tA\n", $newX[$i];
        }
}
