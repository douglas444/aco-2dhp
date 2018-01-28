runs_per_instance=20

for i in `seq 1 8`
do

	echo "Sequence $i" 

	output=""
	best_energy=0
	best_energy_occurrences=0
	time_sum=0
	solution=""
	best_solution=""
	sequence=""

	for j in `seq 1 $runs_per_instance`
	do
		output=$(./bin/Release/aco-2dhp ./input sequence"$i")

		if echo "$output" | grep -q "Error"
		then
			echo "$output"
			exit 1
		fi

		output=$(printf "%s %s %s %s %s %s" $output)
		energy="$(cut -d' ' -f1 <<<$output)"
		time="$(cut -d' ' -f2 <<<$output)"
		solution="$(cut -d' ' -f3 <<<$output)"
		time_sum=$(echo "$time_sum + $time" | bc)


		if [ $energy -eq $best_energy ]
		then
			best_energy_occurrences=$((1+$best_energy_occurrences))
		fi

		if [ $energy -lt $best_energy ]
		then
			best_energy=$energy
			best_solution=$solution
			best_energy_occurrences=1
		fi

	done

	time_avg=$(echo "scale=2; $time_sum / $j" | bc)

	echo " * Best energy = $best_energy"
	echo " * Best energy occurrences = $best_energy_occurrences"
	echo " * Time avg = $time_avg (s)"
	echo " * Solution = $best_solution"

	sequence="$(cut -d' ' -f4 <<<$output)"

	./2dhp-plot "$sequence" "$best_solution" ./sequence"$i"

done
