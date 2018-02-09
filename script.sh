runs_per_instance=20

for i in `seq 1 8`
do

	echo "************************************** Sequence $i" 
	echo ""

	output=""
	best_energy=0
	best_energy_occurrences=0
	time_sum=0
	solution=""
	best_solution=""
	sequence=""

	echo "\begin{tabular}{|l|l|l|l|l|l|l|} \hline"
	echo "EX & E & I & T(s) & AVG & STDDEV & RATE \\\ \hline"

	for j in `seq 1 $runs_per_instance`
	do
		output=$(./bin/Release/aco-2dhp ./input sequence"$i")

		if echo "$output" | grep -q "Error"
		then
			echo "$output"
			exit 1
		fi

		sequence="$(cut -d'|' -f1 <<<$output)"
		directions="$(cut -d'|' -f2 <<<$output)"
		energy="$(cut -d'|' -f3 <<<$output)"
		final_particles_avg="$(cut -d'|' -f4 <<<$output)"
		final_particles_solution_rate="$(cut -d'|' -f5 <<<$output)"
		final_particles_stddev="$(cut -d'|' -f6 <<<$output)"
		found_on_iteration="$(cut -d'|' -f7 <<<$output)"
		time="$(cut -d'|' -f8 <<<$output)"

		echo "$j & $energy & $found_on_iteration & $time & $final_particles_avg & $final_particles_stddev & $final_particles_solution_rate\% \\\ \hline"

		time_sum=$(echo "$time_sum + $time" | bc)

		if [ $energy -eq $best_energy ]
		then
			best_energy_occurrences=$((1+$best_energy_occurrences))
		fi

		if [ $energy -lt $best_energy ]
		then
			best_energy=$energy
			best_solution=$directions
			best_energy_occurrences=1
		fi
	done

	echo "\end{tabular}"
	echo ""

	time_avg=$(echo "scale=2; $time_sum / $j" | bc)

	echo "best energy: $best_energy"
	echo "best energy occurrences: $best_energy_occurrences"
	echo "time avg: $time_avg (s)"
	echo "solution: $best_solution"

	sequence="$(cut -d' ' -f4 <<<$output)"
	./2dhp-plot "$sequence" "$best_solution" ./sequence"$i"

	echo ""

done
