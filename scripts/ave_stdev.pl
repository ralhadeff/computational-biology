# this script takes a list of numbers and prints out the average and the stdev
# usage: list | perl thisScript
# output is first line = averages, second line = stdev

@sum = 0;
@count = 0;

# first pass
while(<STDIN>){
        @words = split;
	$rows = scalar @words;
        for ($i=0;$i<$rows;$i++){
                if ($words[$i]=~/-?\d+\.?\d*/){
                        $sum[$i]+=$words[$i];
                        $value[$i][$count[$i]++]=$words[$i];
                }
	}
}

# second pass
@ss = 0;
for ($j=0;$j<scalar @value;$j++) {
	$ave[$j]=$sum[$j]/$count[$j];
	for ($i=0;$i<$count[$j];$i++){
	        $ss[$j]+=($ave[$j]-$value[$j][$i])**2;
	}
	$stdev[$j] = sqrt($ss[$j]/($count[$j]-1));
}

for ($i=0;$i<scalar @value;$i++) {
	printf "%f\t", $ave[$i]
}
print "\n";
for ($i=0;$i<scalar @value;$i++) {
        printf "%f\t", $stdev[$i]
}
print "\n";
