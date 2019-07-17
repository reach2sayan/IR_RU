import numpy as np
import os

elems = []
sites = []
spec = {}
phases = []
allPhase = []
numOfparams = 0
valueList = ('energy','svib_ht')
phaseList = os.listdir("/users/ssamanta/atat/data/sqsdb/")
def readInput():
    sf = open("species.in","r")
    
    elems = (sf.read().strip("\n")).split(",")
    sf.close()
    files = os.listdir(".")
    tempPhases = [file for file in files if file in phaseList]

    for phase in tempPhases:
        fileSpec = "./{}/species.in".format(phase)
        pf = open(fileSpec,"r")
        rawsite = (pf.read().strip("\n")).split(" ")
        for s in range(0,len(rawsite)):
            asitespec = rawsite[s].split("=")
            
    #count how many parameters to be fit
    #for end-members G terms
    for phase in phases:
        for elem in elems:
            numOfparams = numberOfparams + 2
    #for L terms

    
