#! /usr/bin/env python

# -*-coding:UTF-8 -*

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

import sys, os, string

def usage():
    print("Usage: " + sys.argv[0] + " tolerance outputSummaryStable outputSummaryToBeTested")
    print("Example: " + sys.argv[0] + " 20 outputSummaryStable outputSummaryToBeTested")

def searchAndAppendValues(listToSearch, resultList):
    for line in listToSearch:
        if(line):
            values=line.split(" ")
            if (values[0] == "bcTimeMain:"):
                resultList.append(values[1:])

if (len(sys.argv) < 4):
    usage()
    sys.exit("Not enough parameters")


tolerance=round(float(sys.argv[1]))

staticOutputSummary=open(sys.argv[2])
textStaticOutputSummary=staticOutputSummary.read()
staticOutputSummary.close()

currentOutputSummary=open(sys.argv[3])
textCurrentOutputSummary=currentOutputSummary.read()
currentOutputSummary.close()

returnVal=0

textStaticOutputSummary=textStaticOutputSummary.split("\n")
textCurrentOutputSummary=textCurrentOutputSummary.split("\n")

if (len(textStaticOutputSummary) < 1):
    sys.exit("outputSummaryStable is empty")

if (len(textCurrentOutputSummary) < len(textStaticOutputSummary)):
    sys.exit("outputSummaryToBeTested must be at least the same size of the outputSummaryStable")

stableValues = list()
toBeTestedValues = list()

searchAndAppendValues(textStaticOutputSummary[1:], stableValues)
searchAndAppendValues(textCurrentOutputSummary[1:], toBeTestedValues)

if (len(stableValues) != len(toBeTestedValues)):
    sys.exit("Number of values to be tested are not the same: " + str(len(stableValues)) + " != " + str(len(toBeTestedValues)))

i=0
while(i < len(toBeTestedValues)) :
    stableVal = float(stableValues[i][0])
    testedVal = float(toBeTestedValues[i][0])

    if(stableVal*(1-2*tolerance/100.) >= testedVal):
        print bcolors.OKGREEN + "Better time: " + str(stableVal) + " >= " + str(testedVal) + bcolors.ENDC
    elif((stableVal*(1-tolerance/100.)) >= testedVal):
        print bcolors.OKBLUE + "Good time: " + str(stableVal) + " >= " + str(testedVal) + bcolors.ENDC
    elif((stableVal*(1+2*tolerance/100.)) <= testedVal):
        print bcolors.FAIL + "Worst time: " + str(stableVal) + " <= " + str(testedVal) + bcolors.ENDC
    elif((stableVal*(1+tolerance/100.)) <= testedVal):
        print bcolors.WARNING + "Bad time: " + str(stableVal) + " <= " + str(testedVal) + bcolors.ENDC
    else:
        print "Pratically Same time: " + str(stableVal) + " ~= " + str(testedVal)
    i += 1

os._exit(returnVal)
 



