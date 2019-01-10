# Raphael Alhadeff May 2016 v.1

# simple markov chain matrix multiplier (for manual inspection of the process)

# usage: perl script matrix (matrix as simple text file)

# read in input file
open(IN,$ARGV[0]);

$line = 0;
while(<IN>){
        @{$matrix[$line]} = split , $_;
	$line++;
}
close(IN);

print "Matrix:\n";

# print matrix once
for($x=0; $x<scalar @matrix; $x++){
	for($y=0; $y<scalar @matrix; $y++){
		printf "%4.4f\t" , $matrix[$x][$y];
	}
	print "\n";
}

#keep printing multiplications as long as input is not "quit"
while ($input ne "quit"){
	chomp ($input = <STDIN>);
	# create new matrix
	for($x=0; $x<scalar @matrix; $x++){
	        for($y=0; $y<scalar @matrix; $y++){
			# sum the vector's multiplication
			$value = 0;
		        for($i=0; $i<scalar @matrix; $i++){
				$value+=$matrix[$x][$i]*$matrix[$i][$y];
			}
			$newMatrix[$x][$y]=$value;
	        }
	}
	# now move newMatrix to matrix and print it;
	for($x=0; $x<scalar @matrix; $x++){
        	for($y=0; $y<scalar @matrix; $y++){
			$matrix[$x][$y] = $newMatrix[$x][$y];
                	printf "%4.4f\t" , $matrix[$x][$y];
	        }
	        print "\n";
	}

}
