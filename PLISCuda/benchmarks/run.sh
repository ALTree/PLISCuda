#!/bin/sh

PLISCUDA="/galileo/home/userinternal/adonizet/progetto-PLISCuda/PLISCuda/PLISCuda/pliscuda"

# run benchmark on folder $1
function run_benchmark {
	echo -e "==== bench" $1 "====\n"

	cd $1
	rm -f sim* # cleanup 

	NUM=8

	echo -n "  Running benchmark.. "
	s=`date +%s%N`
	for ((i=1;i<=NUM;i++)); do
		CUDA_VISIBLE_DEVICES="1" $PLISCUDA conf.txt > ./results/out${i}.txt
	done
	e=`date +%s%N`
 	echo -e "done ("$( echo -e "scale=4; ($e - $s)/1000000000" | bc -l ) "s).\n"

	echo "  Average of $NUM runs:"
	
	STEPS=$(cat ./results/out* | grep "steps/" | awk '{sum += $2} END {printf sum/NR "\n"}')
	echo -e "   " $STEPS "steps/s"

	TIMEO=$(cat ./results/out* | grep "elapsed")
	TIMEH=$(echo "$TIMEO" | awk '{sum += $3} END {printf sum/NR "h"}')
	TIMEM=$(echo "$TIMEO" | awk '{sum += $4} END {printf sum/NR "m"}')
	TIMES=$(echo "$TIMEO" | awk '{sum += $5} END {printf sum/NR "s"}')
	echo -e "   " $TIMEH $TIMEM $TIMES

	# cleanup
	echo -e ""
	rm -f sim*
	rm -f ./results/out*
	cd ..
}

## benchmarks driver ##

echo -e "\n### Running PLISCuda benchmarks ###\n"

sg=`date +%s%N`
for D in *; do
	if [ -d "${D}" ]; then
		run_benchmark ${D}
	fi
done
eg=`date +%s%N`

echo -e "=======================\n"
echo -n "    DONE "
echo -e "("$( echo -e "scale=4; ($eg - $sg)/1000000000" | bc -l ) "s)\n"




