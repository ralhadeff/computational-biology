#!/usr/bin/env perl

# Raphael Alhadeff Sept 2018 v.34


# usage: perl script inputFileName
# input file syntax is complex.

# based on http://web.chem.ucsb.edu/~kalju/MonteCarlo_3.html (monte carlo metropolis)

# $positionEnergy[$site][$type][$state] = energy in kcal/mol
# barriers for the MC will be calculated by calculating delta energies

# $positionTransition[$site][$type][$transitionIndex] = destinationSite
# transition index = 0 is the number of allowed transfers, so the code knows how many to count, and 1...n are the indices.

# $particle[$index][$type] = position

# $maxOccupancy[$site][$type]= maxParticles for that site/type.
# $maxGeneralOccupancy[$site]= maxParticles total for that site.

# $states[$state][$site][$type]= numberOfAtoms.

# read in input file

print "\nReading input file:\n###########\n";
open(IN,$ARGV[0]);
# set all flags to false
$cyclic = -1;
$occupancy = -1;
$starting = -1;
$energy = -1;
$alting = -1;
$altEnergies = -1;
$specificAlts = 0;
$differentAlts = 0;
$printRate = -1;
while(<IN>){
	print $_;
	if ($occupancy == -1 && $starting == -1 && $energy==-1 && $alting ==-1){
		if (/^basicMap\s+(.+)/){
			@basicMap = split /=/,$1;
		} elsif (/^siteNames\s+(.+)/){
			@siteNames = split /=/,$1;
	        } elsif (/^particleTypes\s+(.+)/){
			@particleTypes = split /,/,$1;
	        } elsif (/^occupancies\s+(.+)/){
			@generalOccupancies = split /=/,$1;
			$occupancy=0;
		} elsif (/^numberOfSteps\s+(.+)/){
			$totalSteps = $1+0;
		} elsif (/^printFrequency\s+(.+)/){
			$printRate = $1+0;
		} elsif (/^startingConditions/){
			$starting = 0;
		} elsif (/^energies/){
			$energy = 0;
			$totalStates = 0;
		} elsif (/^cyclic/){
			$cyclic = 0;
		} elsif (/^barriers\s+(.+)/){
                        @barriers = split /=/,$1;
		} elsif (/^alternate\s+(\d+)/){
			$alting = 0;
			$alternateFreq = $1;
                } elsif (/^alternateEnergies\s+(.+)/){
                        @generalAltEnergies = split /=/,$1;
			$altEnergies=0;
			# specific energies provided
			$specificAlts=1;
                } elsif (/^temperature\s+(.+)/){
                        $beta = -1.0 / (0.002 * $1);
		}
	} elsif ($occupancy == 0){
		if(/^(\w+):\s+(.+)/){
			# trim spaces
			$1=~s/^\s+|\s+$//g;
			$specificOccupancies{$1} = $2;
		} else {
			$occupancy = -1;
		}
	} elsif ($starting == 0){
		if(/^(\w+):\s+(.+)/){
                        # trim spaces
                        $1=~s/^\s+|\s+$//g;
                        $startingConditions{$1} = $2;
                } else {
                        $starting = -1;
                }
	} elsif ($energy == 0){
		if (/=/){
			$_=~s/^\s+|\s+$//g;
			# first initiation of states matrix
			$rawStates[$totalStates]=$_;
			$energy = 1;
		} else {
			$energy = -1;
		}
	} elsif ($alting > -1){
		if (/(\d+\.?\d*)\s+(.+)/){
			@{$alternates[$alting]} = split /=/, $2;
			$alternateBarrier[$alting] = $1;
			$alting++;
		} else {
			$alting = -1;
		}
	} elsif ($energy == 1){
		if (/^(\w+):\s+(.+)/){
			# trim spaces
                        $1=~s/^\s+|\s+$//g;
                        $energies[$totalStates]{$1} = $2;
		} else {
			$totalStates++;
			$energy = 0;
		}
	}
}
close(IN);
# fix default
if ($printRate<=0){
	$printRate = $totalSteps+1;
}

print "\n###########\nInput file read. Starting to analyze your system\n";
print "###########\n";

if ($beta==0){
	print "No temperature was selected, using 300 as default\n";
	print "###########\n";
	$beta = -1.0 / (0.002 * 300);
}

print "Types detected:\n";
for ($i=0;$i<scalar @particleTypes;$i++){
	$particleTypes[$i]=~s/^\s+|\s+$//g;
	$types[$i]=$particleTypes[$i];
        print "$types[$i]\n";
}

# make some quality control
if (scalar @basicMap != scalar @siteNames || scalar @basicMap != scalar @generalOccupancies){
	print "\nThere is an error in your input! Inconsistency in the number of sites.\n";
	exit;
}

# create the sites and transitions as well as fill the max occupancies. Also populate sites.
$site=0;
$particle=0;
for($i=0; $i<scalar @basicMap; $i++){
	@localNames = split /,/, $siteNames[$i];
	@localGeneralOccupancies = split /,/, $generalOccupancies[$i];
	@localAltEnergies = split /,/, $generalAltEnergies[$i];
	@localBarriers = split /,/, $barriers[$i];
	# stuff for the multiple lines of alternating energy
	for($q=0; $q<scalar @alternates;$q++){
		@{$localAlternates[$q]} = split /,/, $alternates[$q][$i];
	}
	for($j=0; $j<$basicMap[$i];$j++){
		# do for each particle type
		for ($type=0; $type<scalar @types; $type++){
			# create the states matrix
			for ($s=0; $s<$totalStates; $s++){
				@globalStates = split /=/, $rawStates[$s];
				@localStates = split /,/, $globalStates[$i];
				if ($localStates[$j]==-1){
					$states[$s][$site][$type]=-1;
				} else {
					@pattern = $localStates[$j] =~ /$types[$type]/g;
					$count = scalar @pattern;
					$states[$s][$site][$type]=scalar @pattern;
				}
				# create energies array
				@globalEnergies = split /=/, $energies[$s]{$types[$type]};
				@localEnergies = split /,/, $globalEnergies[$i];
				$positionEnergy[$site][$type][$s]= $localEnergies[$j];
			}
                        @starting =  split /=/, $startingConditions{$types[$type]};
			@localStarting = split /,/, $starting[$i];
			@occupancies =  split /=/, $specificOccupancies{$types[$type]};
			# first add transitions to the 'local' sites
			@localOccupancy = split /,/, $occupancies[$i];
			$maxOccupancy[$site][$type] = $localOccupancy[$j];
			# mark barriers;
			$barrier[$site] = $localBarriers[$j];
			# alternation stuff
			$altEnergy[$site] = $localAltEnergies[$j];
			# mark alternating sites
                        for ($q = 0; $q<scalar @alternates;$q++){
		                $alternate[$q][$site] = $localAlternates[$q][$j];
		        }
			# populate site, making sure it will not lead to over population
                        if ($localStarting[$j]>$localOccupancy[$j] && $localOccupancy[$j]!=-1){
                                print "\nYou are trying to over-populate one of your site. Please review your input file.\nAborting\n";
                                exit;
                        } else {
                                for ($k=0; $k<$localStarting[$j]; $k++){
					# write down population array
					$population[$site][$type]++;
                                        $particleType[$particle]=$type;
                                        $particlePosition[$particle]=$site;
                                        $particle++;
                                }
                        }
			# skip sites that cannot be occupied in the first place
			if($localOccupancy[$j]!=0){
				for($local=$j+1; $local<$basicMap[$i];$local++){
					if ($localOccupancy[$local]!=0){
						# add both forward and backward transitions
						$positionTransition[$site][$type][0]++;
						$positionTransition[$site][$type][$positionTransition[$site][$type][0]]=$site+$local-$j;
						$positionTransition[$site+$local-$j][$type][0]++;
						$positionTransition[$site+$local-$j][$type][$positionTransition[$site+$local-$j][$type][0]]=$site;
					}
				}
				# next add transitions to the next sites
				@localOccupancy = split /,/, $occupancies[$i+1];
		                for($local=0; $local<$basicMap[$i+1];$local++){
		                        if ($localOccupancy[$local]!=0){
		                                $positionTransition[$site][$type][0]++;
		                                $positionTransition[$site][$type][$positionTransition[$site][$type][0]]=$site+$local+$basicMap[$i]-$j;
						$positionTransition[$site+$local+$basicMap[$i]-$j][$type][0]++;
						$positionTransition[$site+$local+$basicMap[$i]-$j][$type][$positionTransition[$site+$local+$basicMap[$i]-$j][$type][0]]=$site;
		                        }
		                }
			}
		}
	$siteName[$site]=$localNames[$j];
	$maxGeneralOccupancy[$site]=$localGeneralOccupancies[$j];
	$site++;
	}
}

# print out all transition (quality control)
print "\nDone creating transition array. Please review carefully.\n";
for ($type=0; $type<scalar @types; $type++){
	print "###########\nTransitions for $types[$type]:\n";
	for ($i=0; $i<scalar @siteName; $i++){
		for ($j=1; $j<=$positionTransition[$i][$type][0]; $j++){
			print "From $siteName[$i] to $siteName[$positionTransition[$i][$type][$j]]\n";
		}
	}
}
print "###########\n";

# print out barriers
print "\nBarriers that the code will consider:\n";
for ($i=0; $i<scalar @siteName; $i++){
	if ($barrier[$i] eq "*"){
		print "$siteName[$i] is a barrier\n";
	}
}
print "###########\n";

# alternating code
if ($alternateFreq !=0){
	print "\nCode will run with alternating conformations, with a barrier of $alternateFreq kcal.\n";
	$alternateChance = exp($alternateFreq*$beta);
        for($q=0; $q<scalar @alternates;$q++){
		for ($i=0; $i<scalar @siteName; $i++){
			if ($alternate[$q][$i] eq "+"){
				$plus.= "$siteName[$i]\t";
			} elsif ($alternate[$q][$i] eq "-"){
				$minus.= "$siteName[$i]\t";
			}
		}
	}
	print "Sites on one side: $plus\nSites on the other side: $minus\n";
	print "Increasing the barriers by $alternateBarrier\n";
		$alternateTurn = int(rand(2));
		if($alternateTurn%2== 0){
			$directionOfStart = "peri-Facing";
		} else {
			$directionOfStart = "cyto-Facing";
		}
	print "Will start at $directionOfStart\n\n";

	for ($i=0; $i<scalar @siteName; $i++){
		print "$siteName[$i] increases crossing energy by $altEnergy[$i]\n";
	}

	#check for specific state energies
	for ($s=0; $s<$totalStates; $s++){
		# if there is no + sign in the Conformation energy line, it is a global energy
		if (!($energies[$s]{"C"}=~/(\S+)\+(\S+)/)){
			$stateAltEnergy[$s] = $energies[$s]{"C"};
			if ($stateAltEnergy[$s]!=0){
				$specificAlts=1;
			}
			print "state:$s:\ten:$stateAltEnergy[$s]:\n";
		# otherwise, the energy is different for IF and OF
		# having the first number be the energy for the - conformation and the sum is the energy for both
		} else {
			$stateAltEnergy[$s] = $1;
			$stateAltEnergyDifference[$s] = $2;
			$differentAlts=1;
			$num = $1+$2;
                        print "state:$s:\ten-:$stateAltEnergy[$s]:en+:$num:\n";
		}
	}
print "###########\n";
}

# print out occupancies for quality control
print "\nStarting and Max occupancies considered:\n";
for ($i=0; $i<scalar @siteName; $i++){
	printf "Site %-16s- ", $siteName[$i];
	$sum = 0;
	for ($type=0; $type<scalar @types; $type++){
		$max = $maxOccupancy[$i][$type];
		if ($max==-1){
			$max = "unl";
		}
		$sum+=$population[$i][$type];
		printf "$types[$type]: %-2s/ %-2s\t",$population[$i][$type], $max;
	}
	$max = $maxGeneralOccupancy[$i];
	if ($max==-1){
		$max = "unl";
        } elsif ($sum>$max) {
		print "\nYou are trying to over-populate one of your site. Please review your input file.\nAborting\n";
		exit;
	}
	printf "total: %-2s/ %-2s\n", $sum, $max;
}
print "###########\n";

# print out states and energies for quality control
print "\nDifferent states to be considered and corresponding energies. Review very carefully.\n";
for ($s=0; $s<$totalStates; $s++){
	print "for state $s:\t       state conditions:\tstate energies:\n";
	for ($i=0; $i<scalar @siteName; $i++){
		printf "Site %-16s- ", $siteName[$i];
		for ($type=0; $type<scalar @types; $type++){
			print "$types[$type]: $states[$s][$i][$type]\t";
		}
		print "\t";
                for ($type=0; $type<scalar @types; $type++){
                        print "$types[$type]: $positionEnergy[$i][$type][$s]\t";
                }
		print "\n";
	}
	print "###########\n";
}

# make sure that all values have been initiated
for ($s=0; $s<$totalStates; $s++){
        for ($i=0; $i<scalar @siteName; $i++){
                for ($type=0; $type<scalar @types; $type++){
			$population[$i][$type]+=0;
		}
	}
}

# assess the starting state
$state = -1;
MISMATCH:
for ($s=0; $s<$totalStates; $s++){
        for ($i=0; $i<scalar @siteName; $i++){
		for ($type=0; $type<scalar @types; $type++){
			if ($population[$i][$type]!=$states[$s][$i][$type] && $states[$s][$i][$type]!=-1){
				next MISMATCH;
			}
		}
	}
	$state=$s;
	last;
}

if ($state==-1){
	print "\nYou have a state with no energies defined. Please review your input or add a default state.\nAborting\n";
	exit;
}
print "System has been found to be in state $state\n";

# assess the starting energy
$totalEnergy = 0;
for ($i=0; $i<scalar @siteName; $i++){
	for ($type=0; $type<scalar @types; $type++){
		$totalEnergy+= $population[$i][$type] * $positionEnergy[$i][$type][$state];
	}
}
print "System has an initial total energy of $totalEnergy Kcal\n";

if ($cyclic==0){
	print "System will be run in a cyclic way\n";
        for ($type=0; $type<scalar @types; $type++){
                $movedFirstToLast[$type]=0;
                $transferCost[$type]=$positionEnergy[0][$type][0]-$positionEnergy[scalar @siteName -1][$type][0];
		# make sure these are equal for all states
		for ($s=1; $s<$totalStates; $s++){
			$stateGap = $positionEnergy[0][$type][$s]-$positionEnergy[scalar @siteName -1][$type][$s];
			if ($stateGap != $transferCost[$type]){
				print "\nYour gradients are not equal for all states.\nAborting\n";
				exit;
			}
		}
        }
}
print "###########\n";

print "\nMC process begins and will run for $totalSteps.\n";

# alternating intiation
if ($alternateFreq !=0){
	if ($alternateTurn%2==0){
	        for ($i=0; $i<scalar @siteName; $i++){
		         for($q=0; $q<scalar @alternates;$q++){
		        	if ($alternate[$q][$i] eq "-"){
					for ($type=0; $type<scalar @types; $type++){
						for ($s=0; $s<$totalStates; $s++){
						# increment barriers for each type, states and appropriate site
							$positionEnergy[$i][$type][$s]+=$alternateBarrier[$q];
						}
					}
				}
			}
		}
	} else {
                for ($i=0; $i<scalar @siteName; $i++){
		         for($q=0; $q<scalar @alternates;$q++){	                        
				if ($alternate[$q][$i] eq "+"){
	                                for ($type=0; $type<scalar @types; $type++){
	                                        for ($s=0; $s<$totalStates; $s++){
	                                        # increment barriers for each type, states and appropriate site
	                                               $positionEnergy[$i][$type][$s]+=$alternateBarrier[$q];
	                                        }
	                                }
	                        }
			}
                }

	}
}
$altRoll = 1;

# make a duplicate of the population array
for ($i=0; $i<scalar @siteName; $i++){
	for ($type=0; $type<scalar @types; $type++){
        	$nextPopulation[$i][$type]=$population[$i][$type];
        }
}

# step loop
$onBarrier = 0;
for ($step=1; $step<=$totalSteps; $step++){
	if ($step%$printRate==0){
		print "step $step\n";
		# print so-far results
		print "Final transfers $movedFirstToLast[0] $movedFirstToLast[1]\n";
		# take into account the assymetry of the peri and cytoplasm
		$dummy[0] = $movedFirstToLast[0] + $population[0][0] - $population[scalar @siteName-1][0];
		$dummy[1] = $movedFirstToLast[1] + $population[0][1] - $population[scalar @siteName-1][1];
		print "Clean transfers $dummy[0] $dummy[1]\n";
		print "Total conf changes $alternateTurn\n";
	}
	# handle alternating
	# do not alternate if on a barrier
	if (($alternateFreq !=0) && ($onBarrier ne "*" )){
		# it's time to alternate if:
		$altTotalE = $alternateFreq;
		# did the user provide specific energies for the alternation?
		# otherwise the Freq is the totalE
		if ($specificAlts!=0){
			if ($stateAltEnergy[$state]== 0){
				for ($i=0; $i<scalar @siteName; $i++){
					for ($type=0; $type<scalar @types; $type++){
						$altTotalE+= $population[$i][$type]*$altEnergy[$i];
					}
				}
			} else {
				if ($differentAlts==0){
					# alt energy normal
					$altTotalE = $stateAltEnergy[$state];
				} else {
					# alt energy for - conformation
					$altTotalE = $stateAltEnergy[$state];
					# addition for + conformation
					if ($alternateTurn%2==1){
						$altTotalE+= $stateAltEnergyDifference[$state];
					}
				}
			}
		}
		#print "alternate energy: $altTotalE\n";
		$alternateChance = exp($beta * $altTotalE);
                if ($alternateChance > $altRoll){
			# this prevents changing again if an illegal step was skipped
			$altRoll = 1;
                        $alternateTurn++;
	                #if($alternateTurn%2== 0){
	                #        $directionOfStart = "peri-Facing";
        	        #} else {
                	#        $directionOfStart = "cyto-Facing";
	                #}
                        #print "##\nswitching conformations to $directionOfStart\n##\n";
			if ($alternateTurn%2==0){
			#alternate to --'s
				for ($i=0; $i<scalar @siteName; $i++){
					for($q=0; $q<scalar @alternates;$q++){
						# decrease + side
				                if ($alternate[$q][$i] eq "+"){
	                        			for ($type=0; $type<scalar @types; $type++){
				                                for ($s=0; $s<$totalStates; $s++){
				                                        $positionEnergy[$i][$type][$s]-=$alternateBarrier[$q];
	                        			        }
				                        }
				                }
						# increase - side
	                                        if ($alternate[$q][$i] eq "-"){
	                                                for ($type=0; $type<scalar @types; $type++){
	                                                        for ($s=0; $s<$totalStates; $s++){
	                                                                $positionEnergy[$i][$type][$s]+=$alternateBarrier[$q];
	                                                        }
	                                                }
	                                        }
					}
			        }
				#correct energy for bound particles (if relevant)
				if ($differentAlts=1){
					$totalEnergy+= $stateAltEnergyDifference[$state];
				}
			} else {
			#alternate to ++'s
                                for ($i=0; $i<scalar @siteName; $i++){
					for($q=0; $q<scalar @alternates;$q++){
	                                        # increase + side
       	                                	if ($alternate[$q][$i] eq "+"){
       	                                        	for ($type=0; $type<scalar @types; $type++){
       	                                                	for ($s=0; $s<$totalStates; $s++){
       	                                                        	$positionEnergy[$i][$type][$s]+=$alternateBarrier[$q];
       	                                                 	}
                                                	}
                                        	}
                                        	# decrease - side
	                                        if ($alternate[$q][$i] eq "-"){
        	                                        for ($type=0; $type<scalar @types; $type++){
                	                                        for ($s=0; $s<$totalStates; $s++){
                        	                                        $positionEnergy[$i][$type][$s]-=$alternateBarrier[$q];
                                	                        }
                                        	        }
                                       		}
					}
                                }
				#correct energy for bound particles (if relevant)
				if ($differentAlts=1){
					$totalEnergy-= $stateAltEnergyDifference[$state];
				}
			}
		}
	}
	# pick one particle at random, starting from index 0 or the previous particle, if it's on a barrier
	if ($onBarrier ne "*" ){
		$random = int(rand($particle));
	}
	# where is the particle?
	$currentPosition = $particlePosition[$random];
	#print "Trying to move particle $random ($types[$particleType[$random]]) from $siteName[$currentPosition] ";
	# what are the transition options?
	$transitionOptions = $positionTransition[$currentPosition][$particleType[$random]][0];
	# select one transition at random, note that it starts from 1
	$transitionIndex = int(rand($transitionOptions)+1);
	$transition = $positionTransition[$currentPosition][$particleType[$random]][$transitionIndex];
	#print "to $siteName[$transition] ";
	$nextPopulation[$currentPosition][$particleType[$random]]--;
	$nextPopulation[$transition][$particleType[$random]]++;
	# make sure that new population does not lead to overpopulation
	if ($nextPopulation[$transition][$particleType[$random]]>$maxOccupancy[$transition][$particleType[$random]] && $maxOccupancy[$transition][$particleType[$random]]!=-1){
		#print "\nIllegal move\n";
		$step--;
		# undo changes in nextPopulation
	        $nextPopulation[$currentPosition][$particleType[$random]]++;
	        $nextPopulation[$transition][$particleType[$random]]--;
		next;
	}
	# only bother if 
	if ($maxGeneralOccupancy[$transition]!=-1){
		$overPop=0;
		for ($type=0; $type<scalar @types; $type++){
			$overPop+=$nextPopulation[$transition][$type];
		}
		if ($overPop > $maxGeneralOccupancy[$transition]){
			#print "\nIllegal move\n";
			$step--;
	                # undo changes in nextPopulation
	                $nextPopulation[$currentPosition][$particleType[$random]]++;
        	        $nextPopulation[$transition][$particleType[$random]]--;
			next;
		}
	}
	# get the state after this transition
	$nextState = -1;
	NEXTMISMATCH:
	for ($s=0; $s<$totalStates; $s++){
	        for ($i=0; $i<scalar @siteName; $i++){
	                for ($type=0; $type<scalar @types; $type++){
				# have to emulate a transfter
				if ($nextPopulation[$i][$type]!=$states[$s][$i][$type] && $states[$s][$i][$type]!=-1){
		                	next NEXTMISMATCH;
		                }
	                }
	        }
	        $nextState=$s;
	        last;
	}
	if ($nextState==-1){
       		print "\nYou have a state with no energies defined. Please review your input or add a default state.\nAborting\n";
	        exit;
	}
	#print "(from state $state to $nextState )\n";
	# calculate energy after this transition
	$nextTotalEnergy = 0;
	for ($i=0; $i<scalar @siteName; $i++){
	        for ($type=0; $type<scalar @types; $type++){
	                $nextTotalEnergy+= $nextPopulation[$i][$type] * $positionEnergy[$i][$type][$nextState];
	        }
	}
	# this is the correction to prevent the barriers from changing as the particles leave the sites.
	if ($barrier[$transition] eq "*" && $state!=$nextState){
		$correction = ($positionEnergy[$transition][$particleType[$random]][$state]-$positionEnergy[$transition][$particleType[$random]][$nextState]);
		if($correction!=0){
			print "\nYour input file contains barriers that do not correspond to detailed balance. Please review and correct.\nAborting\n";
			exit;
		}
	}
	# add correction for the ions that were manually moved
	if ($cyclic==0){
		for ($type=0; $type<scalar @types; $type++){
			$nextTotalEnergy+=$movedFirstToLast[$type] * $transferCost[$type];
		}
	}
	#print "Energies are $nextTotalEnergy - $totalEnergy ";
	$barrier = $nextTotalEnergy-$totalEnergy;
	#print "= $barrier Kcal barrier\n";
	# if the site is available
		# get chance of success
		$chance = exp($barrier*$beta);
		# try and get over the barrier, Metropolis algorithm
		$roll = rand();
		#print "rolled $roll, and probability is $chance\n";
		if ($chance > $roll){
			#print "particle moved!\n";
	                $population[$currentPosition][$particleType[$random]]--;
	                $population[$transition][$particleType[$random]]++;
			# no need to update nextPop because these values are already stored
			$state = $nextState;
			$totalEnergy = $nextTotalEnergy;
			$particlePosition[$random] = $transition;
			$onBarrier = $barrier[$transition];
			# equilibrate sides
			if ($cyclic==0){
				# if an atom has been added to one of the sides
				if ($transition==0 && $population[$transition][$particleType[$random]]-$population[scalar @siteName-1][$particleType[$random]]>=2){
					#print "Manually moving atom from $siteName[$transition] to $siteName[scalar @siteName-1]\n";
					$population[$transition][$particleType[$random]]--;
					$population[scalar @siteName-1][$particleType[$random]]++;
					# also update nextPop
                                        $nextPopulation[$transition][$particleType[$random]]--;
                                        $nextPopulation[scalar @siteName-1][$particleType[$random]]++;
					$particlePosition[$random] = scalar @siteName - 1 ;
					$movedFirstToLast[$particleType[$random]]++;
				} elsif ($transition==scalar @siteName-1 && $population[$transition][$particleType[$random]]-$population[0][$particleType[$random]]>=2){ 
                                        #print "Manually moving atom from $siteName[$transition] to $siteName[0]\n";
                                        $population[$transition][$particleType[$random]]--;
                                        $population[0][$particleType[$random]]++;
                                        # also update nextPop
                                        $nextPopulation[$transition][$particleType[$random]]--;
                                        $nextPopulation[0][$particleType[$random]]++;
                                        $particlePosition[$random] = 0 ;
                                        $movedFirstToLast[$particleType[$random]]--;
				}
			}
		} else {
			#print "barrier too high\n";
                	# undo changes in nextPopulation
	                $nextPopulation[$currentPosition][$particleType[$random]]++;
	                $nextPopulation[$transition][$particleType[$random]]--;
		}
	# print out populations
	#for ($type=0; $type<scalar @types; $type++){
	#	printf "%-2s : ", $types[$type];
	#	for ($i=0; $i<scalar @siteName; $i++){
	#		print "$population[$i][$type],";	
	#	}
	#	print "\n";
	#}
	#print "Transfered $movedFirstToLast[0] $movedFirstToLast[1]\n";
        $altRoll = rand();
}

print "\nAll done successfully\n";
print "Final transfers $movedFirstToLast[0] $movedFirstToLast[1]\n";

# take into account the assymetry of the peri and cytoplasm
$movedFirstToLast[0] = $movedFirstToLast[0] + $population[0][0] - $population[scalar @siteName-1][0];
$movedFirstToLast[1] = $movedFirstToLast[1] + $population[0][1] - $population[scalar @siteName-1][1];
print "Clean transfers $movedFirstToLast[0] $movedFirstToLast[1]\n";
print "Total conf changes $alternateTurn\n";
