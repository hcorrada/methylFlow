#!/bin/bash

## input file
##>> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist ;

# dataFlag = 0  >>> read CpG sites from file
# dataFlag > 0 >>>> dataFlag equals the number of cpg sites
# dataFlag < 0 >>> read the data from rest of the file

# freqFlag = 0 >>> randomly choose the frequency of each pattern
# freqFlag = 1 >>> read the frequency of patterns from rest of the file(second line)


if [ "$1" = 1 ];
then

echo " different lambda start"

cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/
########  evaluation for different coverages ########
echo -n "" > evalAvg.txt
echo -n "" > mcf.txt
echo -n "" > weight.txt
echo -n "" > match.txt
echo -n "" > eval.txt
echo -n "" > num.txt

echo -n "" > input.txt
echo 1 757121 230 60 4 1 15 0 50 10 >> input.txt
echo 35 35 15 15 >> input.txt


for k in $(seq 0.001 0.1 20)
do
#k=`echo ""|awk 'END {print 2 ^ '$i' }'`
echo 'lambda'
echo -n $k >> num.txt
echo -n "   " >> num.txt

#echo $i >> evalCoverage.txt
#echo -n "   " >> evalCoverage.txt

echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/eval.txt

for j in {1..20}
do
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/shortRead.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/simPattern.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/patterns.tsv
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/weight.txt
echo -n "" > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/match.txt

#change directory to /cbcb/project-scratch/fdorri/Code/methylFlow/testing/
cd /cbcb/project-scratch/fdorri/Code/methylFlow/testing/


echo "SimulateLambda"
../build/simulator/mfSimulate /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/input.txt /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda



echo "MethylFlowLambda"
../build/methylFlow/methylFlow -i /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/shortRead.txt -o /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda -l $k -chr 0


echo "EvaluateLambda"
../build/evaluation/mfEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda 757121 757353 $k

wc -l  patterns.tsv >> num.txt


done
echo "avgEval Start"
../build/avgEval/avgEvaluation /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda /cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda $k
done

#sed -e "s/$/$i/" eval.txt
echo "avgEval end"
fi


