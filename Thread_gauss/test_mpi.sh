#!/bin/bash
dir="./small_tests"

data_s4="\
a.txt \
a20.txt \
a40.txt \
b.txt \
f.txt \
"

data_s6="\
c.txt \
d.txt \
e.txt \
g.txt \
"

function parse {
	local fl=$1
	local result=$2

	if [ $fl -ge 1 ]
	then
		local error=`echo "${result}" | grep -a 'ERROR'`
		local residual=`echo "${result}" | grep -a 'residual'`

		echo "${residual}"
		echo "${error}"
	else
		echo "${result}"
		echo "-----------------------------------------------------------------------"
	fi
}

if [ -n "$1" ]
then
	fl=$1
	if [ $fl -eq 0 ]
	then echo -e "---------------------SMALL TESTS - ALL RESULTS-----------------------\n "
	fi

	if [ $fl -eq 1 ]
	then echo -e "----------------SMALL TESTS - ONLY RESIDUAL AND ERROR----------------\n "
	fi

	if [ $fl -gt 1 ] && [ $fl -lt 0 ]
	then 
		fl=-1
		echo "Incorrect arg"
	fi
else
	fl=-1
	echo "USAGE:"
	echo "0 - SMALL TESTS - ALL RESULTS"
	echo "1 - SMALL TESTS - ONLY RESIDUAL AND ERROR"
fi


if [ $fl -ne -1 ]
then
	for data in ${data_s4}
	do \
		data_file="${dir}/${data}"

		echo "=================== ${data_file} k = 1 ==================="
		result=`mpirun -np 1 ./a.out 4 3 4 0 ${data_file} | tr -d '\0'`
		parse $fl "${result}"

		echo "=================== ${data_file} k = 2 ==================="
		result=`mpirun -np 2 ./a.out 4 3 4 0 ${data_file} | tr -d '\0'`
		parse $fl "${result}"		
	done

	for data in ${data_s6}
	do \
		data_file="${dir}/${data}"

		echo "=================== ${data_file} k = 1 ==================="
		result=`mpirun -np 2 ./a.out 6 3 4 0 ${data_file} | tr -d '\0'`
		parse $fl "${result}"

		echo "=================== ${data_file} k = 2 ==================="
		result=`mpirun -np 2 ./a.out 6 3 4 0 ${data_file} | tr -d '\0'`
		parse $fl "${result}"	
	done


	for (( s = 1; s <= 4; s++ ))
	do \
		for (( k = 1; k <= 4; k++ ))
		do \
			for (( n = 1; n <= 30; n++ ))
			do \
				if [ "$s" -ne 4 ] || [ "$n" -le 15 ];
				then
					for (( m = 3; m <= 30; m+=3 ))
					do \
						if [ "$m" -le "$n" ]
						then
							echo "===================== n = ${n} m = ${m} s = ${s} k = ${k} ===================== "
							result=`mpirun -np ${k} ./a.out ${n} ${m} 4 ${s} | tr -d '\0'`		
							parse $fl "${result}"
						fi
					done
				fi
			done
		done
	done
fi








