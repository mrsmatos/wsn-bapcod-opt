#!/bin/sh

today=`date +"%m_%d_%Y"`

outputSummary=outputSummary-$today
statistics=statistics-$today

# You need to change the following variable according to the test
test=testSolveMIP
testDir=tests/nonRegressionTests/$test
parameters="-b $testDir/config/bcParameters.cfg -a $testDir/config/dsgapParameters.cfg -i tests/data/gapC-5-100.txt"
application=GeneralizedAssignmentDemo

sh tests/run_application_and_test_if_stable.sh "$outputSummary" "$statistics" "$test" "$testDir" "$parameters" "$application" "$1"
