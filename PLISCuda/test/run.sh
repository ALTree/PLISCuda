#!/bin/sh

PLISCUDA="/galileo/home/userinternal/adonizet/progetto-PLISCuda/PLISCuda/PLISCuda/pliscuda"

# run tests on system contained in folder $1
function run_tests {
	echo -e "====" $1 "====\n"

	echo -n "  Running simulation.. "
	cd $1
	$PLISCUDA conf.txt > /dev/null
	echo -e "done.\n"

	echo -n "  == populations test.. "
	pop_test

	# cleanup
	echo -e ""
	rm sim*
	cd ..
}

# count global popolation of each specie in the simulation output
# file, then compare with the expected populations from pops.txt
function pop_test {
	AWKSUM='{for (i=1;i<=NF;i++) a[i]+=$i} END {for (i=1;i<=NF;i++) printf a[i]; printf "\n"}'
	P1=$(tail -n +5 sim* | awk "$AWKSUM")
	P2=$(cat golden/pops.txt)

	if [ "$P1" == "$P2" ]; then
		echo "PASS."
	else
		echo "FAIL!"
		echo "    want: |$P2|"
		echo "    got:  |$P1|"
	fi
}


## tests driver ##

echo -e "\n### Running PLISCuda tests ###\n"

for D in *; do
    if [ -d "${D}" ]; then
        s=`date +%s%N`
		run_tests ${D}
		e=`date +%s%N`
		echo -e "  Done ("$( echo -e "scale=4; ($e - $s)/1000000000" | bc -l ) "s)\n"
    fi
done



