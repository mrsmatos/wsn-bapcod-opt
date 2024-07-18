#!/bin/bash
set -e
err=0
if [ $# -ne 7 ]
then
    echo "Usage: $0 outputSummary statistics test testDir parameters application saveLog"
    exit 1
fi

test=$3
testDir=$4
parameters=$5
application=$6

if [ $7 -ne 0 ]
then
    outputDest=$testDir/$1.log.out
    outputSummary=$1
    statistics=$2
else
    outputDest="/dev/null"
    outputSummary=$1_toBeDeleted
    statistics=$2_toBeDeleted
fi


# Test if the application exists.
if [ ! -f ./bin/$application ]
then
    echo "ERROR: ./bin/$application is not found."
    echo "       Did you build $application in release mode?"
    echo "       Are you in $application root directory?"
    exit 1
fi

# If the outputSummary file exist, we remove it.
if [ -f $testDir/$outputSummary ]
then
    rm $testDir/$outputSummary
fi

# Run the application
./bin/$application $parameters --DEFAULTPRINTLEVEL=0 --NameOfOutputSummaryFile=$testDir/$outputSummary --statistics=$testDir/$statistics >$outputDest || err=1

# State of the application: OK if stable, KO if not.
state="OK"

# If the ouputSummaryStable file is not found, then we copy outputSummary to outputSummaryStable
if [ ! -f $testDir/outputSummaryStable ]
then
    cp $testDir/$outputSummary $testDir/outputSummaryStable
fi

# Test if the result are the same as the last stable version.
# python -u ../../Scripts/python/TestOutPutSummary/testOutputSummary.py $testDir/$outputSummary $testDir/outputSummaryStable #commented and modified by Issam.
python -u ../../Scripts/python/TestOutputSummary/testOutputSummary.py $testDir/$outputSummary $testDir/outputSummaryStable

# If the previous test is not successful, then the run is not ok else it's stable
if [ $? -ne 0 ]
then
    echo "The current version is NOT stable."
    state="KO"
    echo "To reproduce:" ./bin/$application $parameters
    exit 1
#else
    #echo "The current version is stable."
fi

if [ $7 -ne 0 ]
then
    cd $testDir
    tar czf $1.log.tar.gz $1.log.out $outputSummary $statistics
    rm $1.log.out $outputSummary $statistics
    cd -
fi

exit $err

#cd $testDir

# Get the version from Git
#currentVersion=`git log --pretty="%H" HEAD^..HEAD`

# Add the test to the RUN_LIST file with the name, the version and the state
# echo "$outputSummary;$currentVersion;$state" >> RUN_LIST

# If state is OK then we add this version to STABLE_VERSION file
#if [ $state = "OK" ]
#then
#    todayFormated=`date +"%m-%d-%Y"`
#    echo "$currentVersion;$todayFormated" >> STABLE_VERSION
#fi

# Add the outputSummary and statistics to git
# git add $outputSummary $statistics RUN_LIST STABLE_VERSION
