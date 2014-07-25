#!/bin/bash

### run :  sh cpg.sh par
### par = 0 > simple
### par = 1 > moderate
### par = 2 > Hard


## input file
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

# dataFlag = 0  >>> read CpG sites from file
# dataFlag > 0 >>>> dataFlag equals the number of cpg sites
# dataFlag < 0 >>> read the data from rest of the file

# freqFlag = 0 >>> randomly choose the frequency of each pattern
# freqFlag = 1 >>> read the frequency of patterns from rest of the file(second line)





########  evaluation for different number of CpG sites ########

if [ "$1" == 2 ]
then

echo "Hard Setting"

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/


echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/input.txt
echo 1 757121 230 60 10 1 20 0 $i 10  >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/match.txt


#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCpG"

../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard


echo "MethylFlowCpG"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/ -l 1 -chr 0


echo "EvaluateCpG"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"

../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard $i
echo "avgCpGEval end"

done


elif [ "$1" == 1 ]
then

echo "Moderate Setting"

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/


echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 30 2 120)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/input.txt
echo 1 757121 230 60 4 1 20 0 $i 10  >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/input.txt
echo 15 15 35 35 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/match.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCpG"

../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate


echo "MethylFlowCpG"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/ -l 1 -chr 0


echo "EvaluateCpG"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"

../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate $i
echo "avgCpGEval end"

done


elif [ "$1" == 0 ]
then

echo "Simple Setting"

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/


echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

#echo "var"  "u"   "v"    "weight "  >> weight.txt


for i in $(seq 3 2 120)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/input.txt
echo 1 757121 230 60 2 1 20 0 $i 10  >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/input.txt
echo 25 75 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/input.txt
#echo $i >> evalCpG.txt
#echo -n "   " >> evalCpG.txt
echo $i
for j in {1..100}
do

echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/match.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCpG"

../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple


echo "MethylFlowCpG"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/ -l 1 -chr 0


echo "EvaluateCpG"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple 757121 757353 $i


done
#sed -e "s/$/$i/" eval.txt
echo "avgEval Start"

../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple /cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple $i
echo "avgCpGEval end"

done


else

echo " your input should be 0, 1 or 2"

fi












