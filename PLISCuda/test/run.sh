#!/bin/sh

PLISCUDA="/galileo/home/userinternal/adonizet/progetto-PLISCuda/PLISCuda/PLISCuda/pliscuda"

# run tests on folder $1
function run_tests {
	echo -e "====" $1 "====\n"

	rm -f sim* # cleanup 

	echo -n "  Running simulation.. "
	cd $1
	$PLISCUDA conf.txt > /dev/null
	echo -e "done.\n"

	echo -n "  == populations test.. "
	pop_test

	echo -n "  == memcheck test.. "
	memcheck_test

	# cleanup
	echo -e ""
	rm -f sim*
	cd ..
}

# count global popolation of each specie in the simulation output
# file, then compare with the expected populations from pops.txt
function pop_test {
	AWKSUM='{for (i=1;i<=NF;i++) a[i]+=$i} END {for (i=1;i<=NF;i++) {printf a[i]; printf " "} printf "\n"}'
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

# run simulation with cuda-memcheck
function memcheck_test {
	# change entTime to 0.01 or it will take ages
	cat conf.txt | sed 's/endTime=.*/endTime=0.01/' > conf_t.txt

	# run memcheck and capture output
	MCHECKOUT=`cuda-memcheck $PLISCUDA conf_t.txt`

	if [[ $MCHECKOUT == *"ERROR SUMMARY: 0 errors"* ]]
	then
		echo "PASS.";
	else
		echo "FAIL!"
		echo "     cuda-memcheck output was written to mcheck-fail.txt"
		echo "$MCHECKOUT" > mcheck-fail.txt
	fi

	rm -f conf_t.txt
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

