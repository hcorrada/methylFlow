#!/bin/bash

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
########  evaluation for different coverages ########
echo -n "" > evalAvgCoverage.txt
echo -n "" > mcfCoverage.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 5 3 100)
do
echo -n "" > evalCoverage.txt
echo -n "" > input.txt
echo 1 757121 230 60 10 1 $i 0 50 10 >> input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> input.txt
#echo $i >> evalCoverage.txt
#echo -n "   " >> evalCoverage.txt

echo $i
for j in {1..100}
do

echo "SimulateCoverage"
../build/simulator/mfSimulate < /cbcb/project-scratch/fdorri/Code/methylFlow/testing/input.txt > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt


echo "MethylFlowCoverage"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/ -l 0.05 -chr 0


echo "EvaluateCoverage"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/patterns.tsv /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalCoverage.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/mcfCoverage.txt 757121 757353 $i

done
echo "avgEval Start"
../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalCoverage.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalAvgCoverage.txt $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done

########  evaluation for different read lengthes ########

echo -n "" > evalAvgReadLength.txt
echo -n "" > mcfReadLength.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 10 5 150)
do
echo -n "" > evalReadLength.txt
echo -n "" > input.txt
echo 1 757121 230 $i 10 1 20 0 50 10 >> input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> input.txt
#echo $i >> evalReadLength.txt
#echo -n "   " >> evalReadLength.txt
echo $i
for j in {1..100}
do
echo "SimulateReadLength"
../build/simulator/mfSimulate < /cbcb/project-scratch/fdorri/Code/methylFlow/testing/input.txt > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt


echo "MethylFlowReadLength"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/ -l 0.001 -chr 0


echo "EvaluateReadLength"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/patterns.tsv /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalReadLength.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/mcfReadLength.txt 757121 757353 $i

done

../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalReadLength.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalAvgReadLength.txt $i

#sed -e "s/$/$i/" eval.txt

done


########  evaluation for different number of CpG sites ########

echo -n "" > evalAvgCpG.txt
echo -n "" > mcfCpG.txt
echo -n "" > weight.txt
echo -n "" > match.txt

#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
echo -n "" > evalCpG.txt
echo -n "" > input.txt
echo 1 757121 230 60 10 1 20 0 $i 10  >> input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/patterns.tsv



echo "SimulateCpG"

../build/simulator/mfSimulate < /cbcb/project-scratch/fdorri/Code/methylFlow/testing/input.txt > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt


echo "MethylFlowCpG"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/newOut.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/ -l 0.05 -chr 0


echo "EvaluateCpG"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/patterns.tsv /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalCpG.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/mcfCpG.txt 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"

../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalCpG.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/evalAvgCpG.txt $i
echo "avgCpGEval end"

done


