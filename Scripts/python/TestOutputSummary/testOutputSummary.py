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
    print("Usage: " + sys.argv[0] + " outputSummaryStable outputSummaryToBeTested")

def searchAndAppendValues(listToSearch, resultList):
    for line in listToSearch:
        if(line):
            values=line.split(" ")
            if (len(values) > 2):
                if (values[0] == "CONSTANT"):
                    resultList.append(values[1:])

if (len(sys.argv) < 3):
    usage()
    sys.exit("Not enough parameters")


staticOutputSummary=open(sys.argv[1])
textStaticOutputSummary=staticOutputSummary.read()
staticOutputSummary.close()

currentOutputSummary=open(sys.argv[2])
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
    if(stableValues[i] != toBeTestedValues[i]):
        print bcolors.FAIL + "Different Values: Stable " + str(stableValues[i]) + " != Tested " + str(toBeTestedValues[i]) + bcolors.ENDC
        returnVal += 1
    i += 1

if(returnVal > 0):
    print bcolors.FAIL + "The current version is NOT stable." + bcolors.ENDC
else:
    print bcolors.OKGREEN + "The current version is stable." + bcolors.ENDC


os._exit(returnVal)
 



