#!/bin/bash

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/

echo -n "" > evalCoverage.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 5 5 30)
do
echo -n "" > input.txt
echo 1 757121 230 60 4 1 $i 0 50 10  >> input.txt
echo 15 15 35 35 >> input.txt
echo $i >> evalCoverage.txt
echo $i
echo "SimulateCoverage"
../build/simulator/mfSimulate < /cbcb/project-scratch/fdorri/Code/methylFlow/testing/input.txt > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt


echo "MethylFlowCoverage"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/ -l 0.05 -chr 0


echo "EvaluateCoverage"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/patterns.tsv /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalCoverage.txt 757121 757353 $i



#sed -e "s/$/$i/" eval.txt

done


echo -n "" > evalReadLength.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 40 5 120)
do
echo -n "" > input.txt
echo 1 757121 230 $i 4 1 20 0 50 10  >> input.txt
echo 15 15 35 35 >> input.txt
echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt
echo $i
echo "SimulateReadLength"
../build/simulator/mfSimulate < /cbcb/project-scratch/fdorri/Code/methylFlow/testing/input.txt > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt


echo "MethylFlowReadLength"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/ -l 0.001 -chr 0


echo "EvaluateReadLength"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/patterns.tsv /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalReadLength.txt 757121 757353 $i



#sed -e "s/$/$i/" eval.txt

done



echo -n "" > evalCpG.txt
echo -n "" > weight.txt
echo -n "" > match.txt

#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 70)
do
echo -n "" > input.txt
echo 1 757121 230 60 4 1 20 0 $i 10  >> input.txt
echo 15 15 35 35 >> input.txt
echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
echo "SimulateCpG"
../build/simulator/mfSimulate < /cbcb/project-scratch/fdorri/Code/methylFlow/testing/input.txt > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt


echo "MethylFlowCpG"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/ -l 0.05 -chr 0


echo "EvaluateCpG"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/patterns.tsv /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalCpG.txt 757121 757353 $i



#sed -e "s/$/$i/" eval.txt

done


