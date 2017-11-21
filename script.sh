for i in `seq 1 8`
do

	echo "Sequence $i" 

	best_energy=0
	best_energy_occurrences=0
	time_sum=0

	for j in `seq 1 20`
	do

		output=$(./aco-2dhp configs/"$i" plotfile)
		energy="$(cut -d' ' -f1 <<<"$output")"
		time="$(cut -d' ' -f2 <<<"$output")"

		time_sum=$(echo "$time_sum + $time" | bc)

		if [ $energy -eq $best_energy ]
		then
			best_energy_occurrences=$((1+$best_energy_occurrences))
		fi

		if [ $energy -lt $best_energy ]
		then
			best_energy=$energy
			best_energy_occurrences=1
		fi

	done

	time_avg=$(echo "scale=2; $time_sum / $j" | bc)

	echo " * Time avg = $time_avg (s)"
	echo " * Best energy = $best_energy"
	echo " * Best energy occurrences = $best_energy_occurrences"

done

