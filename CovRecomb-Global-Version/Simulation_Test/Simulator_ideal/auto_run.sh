#!/bin/bash
###
 # @Author: Sonia-Ljy lijysunny@sina.com
 # @Date: 2022-07-07 14:48:45
 # @LastEditors: Sonia-Ljy lijysunny@sina.com
 # @LastEditTime: 2022-07-07 22:11:53
 # @FilePath: /Simulation_Test/Simulator_ideal/auto_run.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 

cd /home/soniali/Desktop/CovRecomb-Global-Version/Simulation_Test/Simulator_ideal


echo "\n Start for simulation datasets generation and Comparison between the 3SEQ and CovRecomb method.\n"

# Only one 'gene_candidate value' could be input at one time, and the recommended value is from 4 to 7 since the elapsed time for 3SEQ to analysis is too long if the sample size is more than 600
for gene_candidate in $*7
do
	python3 Simulator_CovRecombTest.py -ld [15,21,38,46] -gen [$gene_candidate]
	wait
	python3 3seq_run.py -ld [15,21,38,46] -gen [$gene_candidate]
	wait
	python3 Compare_CoverageRate.py -ld [15,21,38,46] -gen [$gene_candidate] -top 6
	wait
	echo "\n End of test. \n"
done

echo `date`

python3 Compare_ElapsedTime.py