#!/bin/bash
runs_per_instance=20

[ -e ./results/summary ] && rm ./results/summary
[ -e ./results/outputlog ] && rm ./results/outputlog
mkdir -p ./results/tables
mkdir -p ./results/figures
mkdir -p ./results/energy_evolution

for i in `seq 1 8`
do
	iteration=0
	output=""
	best_energy=0
	best_energy_occurrences=0
	time_sum=0
	solution=""
	best_solution=""
	sequence=""
	energy_evolution=""
	energy_pairs=()
	energy_tuples=()
	pairs_number=0
	found_first_iteration=0
	x=0
	y=0

	[ -e ./results/tables/sequence$i ] && rm ./results/tables/sequence$i
	mkdir -p ./results/energy_evolution/sequence$i

	#Table file
	echo "\documentclass[" >> ./results/tables/sequence$i
	echo "    border=2mm," >> ./results/tables/sequence$i
	echo "    convert={" >> ./results/tables/sequence$i
	echo "        density=300 -alpha deactivate," >> ./results/tables/sequence$i
	echo "        size=1000x1000," >> ./results/tables/sequence$i
	echo "        outext=.png" >> ./results/tables/sequence$i
	echo "    }" >> ./results/tables/sequence$i
	echo "]{standalone}" >> ./results/tables/sequence$i
	echo "\usepackage{array}" >> ./results/tables/sequence$i
	echo "\begin{document}" >> ./results/tables/sequence$i
	echo "\begin{tabular}{|>" >> ./results/tables/sequence$i
	echo "{\centering\arraybackslash}m{0.3in}|>" >> ./results/tables/sequence$i
	echo "{\centering\arraybackslash}m{0.3in}|>" >> ./results/tables/sequence$i
	echo "{\centering\arraybackslash}m{0.3in}|>" >> ./results/tables/sequence$i
	echo "{\centering\arraybackslash}m{0.3in}|>" >> ./results/tables/sequence$i
	echo "{\centering\arraybackslash}m{0.7in}|>" >> ./results/tables/sequence$i
	echo "{\centering\arraybackslash}m{0.7in}|>" >> ./results/tables/sequence$i
	echo "{\centering\arraybackslash}m{0.7in}|} \hline" >> ./results/tables/sequence$i
	echo "Run & \$E\$ & \$i\$ & T(s) & Final Population Avg & Final Population Std Dev & Final Population \$E\$ Rate \\\ \hline" >> ./results/tables/sequence$i

	#Console
	echo "Sequence $i"

	#Output log file
	echo "Sequence $i" >> ./results/outputlog

	for j in `seq 1 $runs_per_instance`
	do
		#Console
		echo "Run $j"

		#Execute program
		output=$(sudo ./bin/Release/aco-2dhp ./input sequence"$i")
		if echo "$output" | grep -q "Error"
		then
			echo "$output"
			exit 1
		fi

		#Output log file
		echo " - Run $j: $output" >> ./results/outputlog

		#Extract results
		iterations="$(cut -d'|' -f1 <<<$output)"
		sequence="$(cut -d'|' -f2 <<<$output)"
		directions="$(cut -d'|' -f3 <<<$output)"
		energy="$(cut -d'|' -f4 <<<$output)"
		final_particles_avg="$(cut -d'|' -f5 <<<$output)"
		final_particles_stddev="$(cut -d'|' -f6 <<<$output)"
		final_particles_solution_rate="$(cut -d'|' -f7 <<<$output)"
		found_on_iteration="$(cut -d'|' -f8 <<<$output)"
		time="$(cut -d'|' -f9 <<<$output)"
		energy_evolution="$(cut -d'|' -f10 <<<$output)"

		#Energy evolution function file
		pairs_number=$(grep -o "," <<<"$energy_evolution" | wc -l)
		[ -e ./results/energy_evolution/sequence$i/run$j ] && rm ./results/energy_evolution/sequence$i/run$j
		echo "\documentclass[" >> ./results/energy_evolution/sequence$i/run$j
		echo "    border=2mm," >> ./results/energy_evolution/sequence$i/run$j
		echo "    convert={" >> ./results/energy_evolution/sequence$i/run$j
		echo "        density=300 -alpha deactivate," >> ./results/energy_evolution/sequence$i/run$j
		echo "        size=500x500," >> ./results/energy_evolution/sequence$i/run$j
		echo "        outext=.png" >> ./results/energy_evolution/sequence$i/run$j
		echo "    }" >> ./results/energy_evolution/sequence$i/run$j
		echo "]{standalone}" >> ./results/energy_evolution/sequence$i/run$j
		echo "\usepackage{pgfplots}" >> ./results/energy_evolution/sequence$i/run$j
		echo "\pgfplotsset{compat=1.12}" >> ./results/energy_evolution/sequence$i/run$j
		echo "\begin{document}" >> ./results/energy_evolution/sequence$i/run$j
		echo "\begin{tikzpicture}" >> ./results/energy_evolution/sequence$i/run$j
		echo "\begin{axis}[" >> ./results/energy_evolution/sequence$i/run$j
		echo "    xlabel=\$Iteration$,ylabel=\$Energy$,grid=both," >> ./results/energy_evolution/sequence$i/run$j
		if [ $pairs_number -eq 1 ]; then
			echo "    ytick={0,$energy}," >> ./results/energy_evolution/sequence$i/run$j
    			echo "    ymax=0," >> ./results/energy_evolution/sequence$i/run$j
		else
			echo "    ytick={0,-1,...,$energy}," >> ./results/energy_evolution/sequence$i/run$j
		fi
		echo "    xtick={0,20,...,$iterations}," >> ./results/energy_evolution/sequence$i/run$j
		echo "    xticklabels={" >> ./results/energy_evolution/sequence$i/run$j
		for k in `seq 1 $pairs_number`
		do
			if [ "$k" -gt 1 ]; then
				echo "        ,," >> ./results/energy_evolution/sequence$i/run$j
			fi
			energy_pairs[$k]="$(cut -d'/' -f$k <<<$energy_evolution)"
			x="$(cut -d',' -f1 <<<${energy_pairs[$k]})"
			y="$(cut -d',' -f2 <<<${energy_pairs[$k]})"
			energy_tuples[$k]="("$(($(($k - 1)) * 20 * 2))",$y)"
			echo "        $x" >> ./results/energy_evolution/sequence$i/run$j
		done
		if [ $pairs_number -eq 1 ]; then
			found_first_iteration=1
			pairs_number=$(($pairs_number+1))
			energy_tuples[2]="("$((1 * 20 * 2))",$y)"
			echo "        ,," >> ./results/energy_evolution/sequence$i/run$j
			echo "        $iterations" >> ./results/energy_evolution/sequence$i/run$j
		fi
		echo "    }" >> ./results/energy_evolution/sequence$i/run$j
		echo "]" >> ./results/energy_evolution/sequence$i/run$j
		echo "\addplot+[const plot, no marks, thick]" >> ./results/energy_evolution/sequence$i/run$j
		echo "coordinates {" >> ./results/energy_evolution/sequence$i/run$j
		for k in `seq 1 $pairs_number`
		do
			echo "    ${energy_tuples[$k]}" >> ./results/energy_evolution/sequence$i/run$j
		done
		echo "};" >> ./results/energy_evolution/sequence$i/run$j
		if [ $found_first_iteration -eq 1 ]; then
			found_first_iteration=0
			echo "\draw [dashed] (20,0) -- (20,"$((${energy_pairs[1]} - 1))");" >> ./results/energy_evolution/sequence$i/run$j
		else 
			for k in `seq 0 $(($pairs_number - 2))`
			do
				echo "\draw [dashed] ("$(($k * 20 * 2 + 20))","$(($energy - 1))") -- ("$(($k * 20 * 2 + 20))","$((${energy_pairs[1]} + 1))");" >> ./results/energy_evolution/sequence$i/run$j
			done
		fi
		echo "\end{axis}" >> ./results/energy_evolution/sequence$i/run$j
		echo "\end{tikzpicture}" >> ./results/energy_evolution/sequence$i/run$j
		echo "\end{document}" >> ./results/energy_evolution/sequence$i/run$j

		#Table file
		echo "$j & $energy & $found_on_iteration & $time & $final_particles_avg & $final_particles_stddev & $final_particles_solution_rate\% \\\ \hline" >> ./results/tables/sequence$i

		#Summary file
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
	


	#Table file
	echo ""
	echo "\end{tabular}" >> ./results/tables/sequence$i
	echo "\end{document}" >> ./results/tables/sequence$i


	#Summary file
	time_avg=$(echo "scale=2; $time_sum / $j" | bc)
	echo "Sequence $i" >> ./results/summary
	echo "Best energy: $best_energy;" >> ./results/summary
	echo "Best energy occurrences: $best_energy_occurrences;" >> ./results/summary
	echo "Time average: $time_avg (s);" >> ./results/summary
	echo "Solution: $best_solution;" >> ./results/summary
	echo "" >> ./results/summary

	#Figure file
	./2dhp-plot "$sequence" "$best_solution" ./results/figures/sequence"$i"


done
