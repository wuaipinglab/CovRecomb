#!/bin/bash
cd /home/soniali/Desktop/CovRecomb-Global-Version/Simulation_Test/Simulator_compare3SEQ_renew

echo "\n Start for simulation datasets generation and Comparison between the 3SEQ and CovRecomb method.\n"

# Ideal simulation, with less number of differential feature mutations; Figures S3A-B
do
	python3 Simulator_CovRecombTest.py -smp 6 -sr 1
	wait
	python3 3seq_run.py 
	wait
	python3 Compare_CoverageRate.py
	wait
	echo "\n End of test. \n"
done

echo `date`

# Ideal simulation, with more number of differential feature mutations; Figures S3C-D
do
	python3 Simulator_CovRecombTest.py -smp 8 -sr 1
	wait
	python3 3seq_run.py 
	wait
	python3 Compare_CoverageRate.py
	wait
	echo "\n End of test. \n"
done

echo `date`

# Analog simulation, with less number of differential feature mutations; Figures S3E-F
do
	python3 Simulator_CovRecombTest.py -smp 6 -sr 10
	wait
	python3 3seq_run.py 
	wait
	python3 Compare_CoverageRate.py
	wait
	echo "\n End of test. \n"
done

echo `date`

# Analog simulation, with more number of differential feature mutations; Figures S3G-H
do
	python3 Simulator_CovRecombTest.py -smp 8 -sr 10
	wait
	python3 3seq_run.py 
	wait
	python3 Compare_CoverageRate.py
	wait
	echo "\n End of test. \n"
done

echo `date`

# python3 Compare_ElapsedTime.py