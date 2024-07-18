#!/bin/sh

#This script must be run in BapcodProject folder (the root folder).

APPLICATIONS="
CapacitatedVehicleRouting
GeneralizedAssignment
VertexColoring
MultiCommodityFlow
BinPackingAndCuttingStock
"

DEMOS="
CuttingStockDemo  
GeneralizedAssignmentDemo  
VehicleRoutingWithTimeWindowsDemo  
VertexColoringDemo
"

cpu=`getconf _NPROCESSORS_ONLN`

set -e # Exit if any error

testedapps=""
errfile="err-file-app"
rm -f $errfile

# Configure the projects in Release mode
#cmake . -DCMAKE_BUILD_TYPE=Release


# Build app in parallel

for app in ${APPLICATIONS}
do
    # if the application folder exists
    if [ -d Applications/$app ]
    then
		# Go inside the application
	echo "=== Building $app ==="
	( ( cd Applications/$app && make clean && make -j$cpu ) > $app.log    2>&1 || (echo ; echo "=== ERROR on make $app" ; cat $app.log    ; touch $errfile ; exit 1) ) &
    fi
done

for demo in ${DEMOS}
do
    # if the application folder exists
    if [ -d Demos/$demo ]
    then
	# Go inside the application
	echo "=== Building $demo ==="
	( ( cd Demos/$demo && make clean && make -j$cpu ) > $demo.log    2>&1 || (echo ; echo "=== ERROR on make $demo" ; cat $demo.log    ; touch $errfile ; exit 1) ) &
    fi
done

wait

if [ -e $errfile ]; then
    echo "+++++ Error detected in applications compilation +++++"
    #cat *.log
    exit 1
fi

# Testing in parallel
err=0
for app in ${APPLICATIONS}
do
    # if the application folder exists
    if [ -d Applications/$app ]
    then
	echo "TESTING $app TESTING"
	testedapps="$testedapps $app"
	# Go inside the application
	time ( ( cd Applications/$app && sh tests/runTests.sh 1 ) > $app-test.log 2>&1 || (echo ; echo "=== ERROR on test $app" ; cat $app-test.log ; find Applications/$app/tests/ -name '*.out' -ls -exec cat {} \;  ; touch $errfile ; exit 1) ) || err=1 #&
	# Git Commit
	## cd tests
	
        #ONCE IN A WHILE THE SCRIPT SHOULD BE RUN WITH TWO LINES UNCOMMENTED
	#git commit -m "Tests are updated"
	#git push
	
        # Go back to the project root folder.
	##cd ../../..
    fi
done

for demo in ${DEMOS}
do
    # if the application folder exists
    if [ -d Demos/$demo ]
    then
	echo "TESTING $demo TESTING"
	testedapps="$testedapps $demo"
	# Go inside the application
	time ( ( cd Demos/$demo && sh tests/runTests.sh 1 ) > $demo-test.log    2>&1 || (echo ; echo "=== ERROR on test $demo" ; cat $demo-test.log ; find Demos/$demo/tests/ -name '*.out' -ls -exec cat {} \; ; touch $errfile ; exit 1) ) || err=1 #&
	# Git Commit
	## cd tests
	
        #ONCE IN A WHILE THE SCRIPT SHOULD BE RUN WITH TWO LINES UNCOMMENTED
	#git commit -m "Tests are updated"
	#git push
	
        # Go back to the project root folder.
	##cd ../../..
    fi
done

#wait

if [ -e $errfile ]; then
    echo "+++++ Error detected in applications/demos test +++++"
    #cat *.log
    exit 1
fi

touch no_app_is_no_error.log
echo "ALL IS FINE :" $testedapps
cat *.log
exit $err
