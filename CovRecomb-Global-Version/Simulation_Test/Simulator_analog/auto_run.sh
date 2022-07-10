cd /home/soniali/Desktop/CovRecomb-Global-Version/Simulation_Test/Simulator_analog


echo "\n Start for simulation datasets generation and Comparison between the 3SEQ and CovRecomb method.\n"

# Only one 'sr_candidate value' could be input at one time, and the recommended value is from 50 to 100

for sr_candidate in $*100
do
	python3 Simulator_CovRecombTest_sampler.py -ld [15,21,38,46] -gen [8] -sr $sr_candidate
	wait
	python3 3seq_run.py -ld [15,21,38,46] -gen [8]
	wait
	python3 Compare_CoverageRate.py -ld [15,21,38,46] -gen [8] -top 6 -sr $sr_candidate
	wait
	echo "\n End of test. \n"
done

echo `date`

python3 Compare_ElapsedTime.py