#!/bin/bash 

# specific script for building a working folder for generating an energy profile with Molaris pdld/S-LRA

line=$(pwd)

frames=$(ls pdbs/start_*|wc|awk '{print $1}')
let frames=$frames-1

cd $line
cd pdbs
if [ -e "start_00.pdb" ]; then
	for i in {0..9}
		do
			mv start_0$i.pdb start_$i.pdb
		done
fi

for (( i=0; i<=$frames; i++ ))
	do
		echo "$i"
		cd $line
		if ! [ -e "frame_${i}" ]; then
			mkdir frame_${i}
		fi
		for j in {1..10}
			do
				for k in {1..3}
					do
				                if [ "$k" -eq "1" ]; then
				                        text="charged"
				                elif [ "$k" -eq "2" ]; then
				                        text="asp_only"
                                                elif [ "$k" -eq "3" ]; then
                                                        text="neutral"
				                else
				                        text="bla"
				                fi
						let m=j-1
			                	cd $line
						cd frame_${i}/
                                                if [ "$k" -eq "1" ]; then
					                if ! [ -e "conf_${j}" ]; then
		        		        	        mkdir conf_${j}
							fi
                                                fi
                		        	cd conf_${j}
						if ! [ -e "p_$text" ]; then
	                	        		mkdir p_$text
						fi
		                	        cd p_$text
						cp $line/../../lib/lib_rc .
        		                	cp ../../../../build/${text}.inp .
                                                sed -i "s/np 11000/np 1${m}000/" ${text}.inp
						if [ -f ../../../pdbs/em_${i}/em.pdb ]; then
                	                        	n=$(grep '8809  M61 ZD4' ../../../pdbs/em_${i}/em.pdb | awk '{print substr($0,47,8)}')
				                	sed -i "s/start_NN/em_${i}\/em/" ${text}.inp
						else
                	                	        n=$(grep '8809  M61 ZD4' ../../../pdbs/start_${i}.pdb | awk '{print substr($0,47,8)}')
							sed -i "s/start_NN/start_${i}/" ${text}.inp	
						fi
		                                sed -i "s/TT/$n/" ${text}.inp
					done
			done
	done
