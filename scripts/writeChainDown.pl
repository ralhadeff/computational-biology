# this script takes a pdb where the resid or the chains are incorrect (usually due to editting by a server of software such as MOLARIS and modifies it based on user input
# user should edit the pdb file (vim etc) and renumber the FIRST atom of each chain and provide a chain identifier. Then the script will use those numbers to fill in all the 'below' gaps

open(IN,$ARGV[0]);

$chain = "DUMMY";
$id = -1;
$current = -1;
$first = 1;

while(<IN>){
	if (/(ATOM.{17})(.)(....)(.*)/){
		if ($chain eq $2 || $2 eq " "){
			# check if ID needs to be increased
			if ($3>$current && $first!=1){
				$id++;	
			}
			$current = $3;
			$first = 0 ;
		} else {
			# update chain
			$chain = $2;
			$id = $3;
			$first = 1;
		}
		printf "$1${chain}%4d$4\n",$id;
	} else {
		print $_;
	}
}

close(IN);
