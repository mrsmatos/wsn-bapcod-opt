#!/bin/sh

if [ $# -ge 1 ]
then
    saveLog=$1
else
    saveLog="0"
fi

print_and_run()
{
    script=$1
	chmod +x $script
    echo $script
    sh $script $saveLog
}


print_and_run tests/nonRegressionTests/testBranchAndPrice/runTest.sh
print_and_run tests/nonRegressionTests/testColGenAutoSmoothing/runTest.sh

