#!/bin/sh

PLISCUDA="/galileo/home/userinternal/adonizet/progetto-PLISCuda/PLISCuda/PLISCuda/pliscuda"
FAILURES=""

# run tests on folder $1
function run_tests {
	echo -e "====" $1 "====\n"

	cd $1
	rm -f sim* # cleanup 

	echo -n "  Running simulation.. "
	CUDA_VISIBLE_DEVICES="1" $PLISCUDA conf.txt > /dev/null
	echo -e "done.\n"

	echo -n "  == populations test.. "
	pop_test

	echo -n "  == endstate test.. "
	endstate_test

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
	if [ ! -e golden/pops.txt ]; then
		echo "SKIPPED."
		return
	fi

	AWKSUM='{for (i=1;i<=NF;i++) a[i]+=$i} END {for (i=1;i<=NF;i++) {printf a[i]; printf " "} printf "\n"}'
	P1=$(tail -n +5 sim* | awk "$AWKSUM")
	P2=$(cat golden/pops.txt)

	if [ "$P1" == "$P2" ]; then
		echo "PASS."
	else
		echo "FAIL!"
		echo "    want: |$P2|"
		echo "    got:  |$P1|"
		FAILURES="FAIL" 
	fi
}

# check if the end state is equal to the expected one in
# golden/endstate.txt
function endstate_test {
	if [ ! -e golden/endstate.txt ]; then
		echo "SKIPPED."
		return
	fi

	F1=$(cat sim*)
	F2=$(cat golden/endstate.txt)
	
	if [ "$F1" == "$F2" ]; then
		echo "PASS."
	else
		echo "FAIL!"
		echo "    want: |$F2|" | head -n4
		echo "    got:  |$F1|" | head -n4
		FAILURES="FAIL" 
	fi
}

# run simulation with cuda-memcheck
function memcheck_test {
	# change entTime to 0.01 or it will take ages
	cat conf.txt | sed 's/endTime=.*/endTime=0.01/' > conf_t.txt

	# run memcheck and capture output
	MCHECKOUT=`CUDA_VISIBLE_DEVICES="1" cuda-memcheck $PLISCUDA conf_t.txt`

	if [[ $MCHECKOUT == *"ERROR SUMMARY: 0 errors"* ]]
	then
		echo "PASS.";
	else
		echo "FAIL!"
		echo "     cuda-memcheck output was written to mcheck-fail.txt"
		echo "$MCHECKOUT" > mcheck-fail.txt
		FAILURES="FAIL" 
	fi

	rm -f conf_t.txt
}


## tests driver ##

echo -e "\n### Running PLISCuda tests ###\n"

sg=`date +%s%N`
for D in *; do
    if [ -d "${D}" ]; then
		s=`date +%s%N`
		run_tests ${D}
		e=`date +%s%N`
		echo -e "  Done ("$( echo -e "scale=4; ($e - $s)/1000000000" | bc -l ) "s)\n"
    fi
done
eg=`date +%s%N`

echo -e "================================\n"
if [[ -z $FAILURES ]]; then
	echo -n "    ALL TESTS PASSED "
else
	echo -n "    FUCK! SOME TEST FAILED! "
fi

echo -e "("$( echo -e "scale=4; ($eg - $sg)/1000000000" | bc -l ) "s)\n"


