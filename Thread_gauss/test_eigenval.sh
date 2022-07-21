#!/bin/bash
dir="../tests/small_tests"
const_param="6 1e-14"

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
		local residual=`echo "${result}" | grep -a 'Residual1'`

		echo "${error}"
		echo "${residual}"
	else
		echo "${result}"
	fi
	echo "-------------------------------------------------------------------"
}

if [ -n "$1" ]
then
	fl=$1
	if [ $fl -eq 0 ]
	then echo "---------------------SMALL TESTS - ALL RESULTS-----------------------"
	fi

	if [ $fl -eq 1 ]
	then echo "----------------SMALL TESTS - ONLY RESIDUAL AND ERROR----------------"
	fi

	if [ $fl -eq 2 ]
	then echo "-----------------BIG TESTS - ONLY RESIDUAL AND ERROR-----------------"
	fi

	if [ $fl -gt 2 ] && [ $fl -lt 0 ]
	then 
		fl=-1
		echo "Incorrect arg"
	fi
else
	fl=-1
	echo "USAGE:"
	echo "0 - SMALL TESTS - ALL RESULTS"
	echo "1 - SMALL TESTS - ONLY RESIDUAL AND ERROR"
	echo "2 - BIG TESTS   - ONLY RESIDUAL AND ERROR"
fi


if [ $fl -eq 0 ] || [ $fl -eq 1 ]
then
	for data in ${data_s4}
	do \
		data_file="${dir}/${data}"
		echo "====================== ${data_file} ======================"
		result=`./a.out 4 ${const_param} 0 ${data_file} | tr -d '\0'`
		parse $fl "${result}"
	done

	for data in ${data_s6}
	do \
		data_file="${dir}/${data}"
		echo "====================== ${data_file} ======================"
		result=`./a.out 6 ${const_param} 0 ${data_file} | tr -d '\0'`
		parse $fl "${result}"
	done


	for (( s = 1; s <= 4; s++ ))
	do \
		for (( n = 1; n <= 30; n++ ))
		do \
			echo "======================== n = ${n} s = ${s} ======================== "
			result=`./a.out ${n} ${const_param} ${s}| tr -d '\0'`
			parse $fl "${result}"
		done
	done
fi

if [ $fl -eq 2 ]
then
	for (( s = 1; s <= 4; s++ ))
	do \
		for (( n = 1000; n <= 4000; n+=1000 ))
		do \
			if [ $n -ne 3000 ]
			then 
				echo "======================== n = ${n} s = ${s} ======================== "
				result=`./a.out ${n} ${const_param} ${s}| tr -d '\0'`
				parse $fl "${result}"
			fi
		done
	done
fi






