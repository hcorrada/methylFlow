#!/bin/bash

### run :  sh coverage.sh par1 par2

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard



## input file
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

# dataFlag = 0  >>> read CpG sites from file
# dataFlag > 0 >>>> dataFlag equals the number of cpg sites
# dataFlag < 0 >>> read the data from rest of the file

# freqFlag = 0 >>> randomly choose the frequency of each pattern
# freqFlag = 1 >>> read the frequency of patterns from rest of the file(second line)

if [ "$1" == 1 ];
then

if [ "$2" == 2 ];
then


echo "Hard Setting"

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/
########  evaluation for different coverages ########
echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 5 3 100)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/input.txt

#echo -n "" > eval.txt
#echo -n "" > input.txt
echo 1 757121 230 60 10 1 $i 0 50 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/input.txt
#echo $i >> evalCoverage.txt
#echo -n "   " >> evalCoverage.txt

echo $i
for j in {1..100}
do


echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/match.txt



#change directory to  /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCoverage"
../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard


echo "MethylFlowCoverage"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/  -l 1 -chr 0


echo "EvaluateCoverage"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard 757121 757353 $i

done
echo "avgEval Start"
../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done



elif [ "$2" == 1 ];
then

echo "Moderate Setting"


cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/
########  evaluation for different coverages ########
echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 5 3 100)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/input.txt

echo 1 757121 230 60 4 1 $i 0 50 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/input.txt
echo 15 15 35 35 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/input.txt
#echo $i >> evalCoverage.txt
#echo -n "   " >> evalCoverage.txt

echo $i
for j in {1..100}
do


echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/match.txt




#change directory to  /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCoverage"
../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate

echo "MethylFlowCoverage"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate  -l 1 -chr 0


echo "EvaluateCoverage"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate 757121 757353 $i

done
echo "avgEval Start"
../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done


elif [ "$2" == 0 ];
then


echo "Simple Setting"

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/
########  evaluation for different coverages ########
echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 5 3 100)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/input.txt
echo 1 757121 230 60 2 1 $i 0 50 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/input.txt
echo 25 75 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/input.txt
#echo $i >> evalCoverage.txt
#echo -n "   " >> evalCoverage.txt

echo $i
for j in {1..100}
do


echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/match.txt



#change directory to  /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCoverage"
../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple


echo "MethylFlowCoverage"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/  -l 1 -chr 0


echo "EvaluateCoverage"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple 757121 757353 $i

done
echo "avgEval Start"
../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done

else
echo " your input should be 0 , 1 or 2"
fi



elif [ "$1" == 0 ];
then

if [ "$2" == 2 ];
then


echo "Hard Setting"

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/
########  evaluation for different coverages ########
echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 5 3 100)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/input.txt

#echo -n "" > eval.txt
#echo -n "" > input.txt
echo 1 757121 230 60 10 1 $i 0 50 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/input.txt
echo 10 10 10 10 10 10 10 10 10 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/input.txt
#echo $i >> evalCoverage.txt
#echo -n "   " >> evalCoverage.txt

echo $i
for j in {1..100}
do


echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/match.txt



#change directory to  /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCoverage"
../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto


echo "MethylFlowCoverage"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/  -chr 0


echo "EvaluateCoverage"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto 757121 757353 $i

done
echo "avgEval Start"
../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done



elif [ "$2" == 1 ];
then

echo "Moderate Setting"


cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/
########  evaluation for different coverages ########
echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 5 3 100)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/input.txt

echo 1 757121 230 60 4 1 $i 0 50 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/input.txt
echo 15 15 35 35 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/input.txt
#echo $i >> evalCoverage.txt
#echo -n "   " >> evalCoverage.txt

echo $i
for j in {1..100}
do


echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/match.txt




#change directory to  /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCoverage"
../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto

echo "MethylFlowCoverage"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto  -chr 0


echo "EvaluateCoverage"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto 757121 757353 $i

done
echo "avgEval Start"
../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done


elif [ "$2" == 0 ];
then


echo "Simple Setting"

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/
########  evaluation for different coverages ########
echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt

for i in $(seq 5 3 100)
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/eval.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/input.txt
echo 1 757121 230 60 2 1 $i 0 50 10 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/input.txt
echo 25 75 >> /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/input.txt
#echo $i >> evalCoverage.txt
#echo -n "   " >> evalCoverage.txt

echo $i
for j in {1..100}
do


echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/match.txt



#change directory to  /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateCoverage"
../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto


echo "MethylFlowCoverage"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/  -chr 0


echo "EvaluateCoverage"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto 757121 757353 $i

done
echo "avgEval Start"
../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto /cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto $i

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
done

else
echo " your input should be 0 , 1 or 2"
fi


else
echo " your input should be 0 if it is automatic and 1 if lambda is hard coded"
fi