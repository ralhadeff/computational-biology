#!/usr/bin/env perl

# This script takes a pdb file, and renumbers all the chain to whatever is given
# New chain name must be a single character
# use lower case s for space

# Usage: perl <this script> input.pdb new_chain_name

$pdb	= $ARGV[0];
$chain	= $ARGV[1];

# checks that file exists
if (!(-e $pdb)){
	print "\npdb file does not exist!\n";
	exit;
}

# checks $chain
if ($chain eq "s"){
	$chain = " ";
} elsif ((length($chain)!=1) || !($chain=~/[A-Z]/)) {
	print "\nChain identifier must be a single capital letter\n";
	exit;
}

# prompt
print "\n********************\n$pdb is going to be OVERWRITTEN! are you sure you want to continue? (type c for copy)\n";

# last chance to cancel
chomp ($reply = <STDIN>);
if (!($reply eq "yes") && !($reply eq "y") && !($reply eq "c")) {
	print "\nExit requested by user\n";
	exit;
}

# creates the new pdb
system ("cat $pdb | awk '/^ATOM/{\$0=substr(\$0,1,21)\"${chain}\"substr(\$0,23,length(\$0))}1' > ${pdb}_copy.pdb");
# overwrite the old one
if (!($reply eq "c")){ 
	system ("mv ${pdb}_copy.pdb $pdb");
}

print "\nSuccessful\n";

