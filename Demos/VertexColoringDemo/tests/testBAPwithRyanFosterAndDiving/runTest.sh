#!/bin/sh

today=`date +"%m_%d_%Y"`

outputSummary=outputSummary-$today
statistics=statistics-$today

# You need to change the following variable according to the test
test=testBAPwithRyanFosterAndDiving
testDir=tests/$test
parameters="-b $testDir/config/bc.cfg -a $testDir/config/app.cfg -i data/myciel3.col"
application=VertexColoringDemo

sh tests/run_application_and_test_if_stable.sh "$outputSummary" "$statistics" "$test" "$testDir" "$parameters" "$application" "$1"
