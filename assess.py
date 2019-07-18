import numpy as np
import os

elems = [] #holds all elements
sites = []
spec = {} #holds species hashed by site 
phases = [] #holds all phases calculated
terms = {}
allPhase = {} #holds all sites for all phase
numOfparams = 0 #records no. of params to fit with OC
valueList = ('energy','svib_ht')
phaseList = os.listdir("/users/ssamanta/atat/data/sqsdb/")

def readInput():
    sf = open("species.in","r")
    
    elems = (sf.read().strip("\n")).split(",")
    sf.close()
    files = os.listdir(".")
    tempPhases = [file for file in files if file in phaseList]

    for phase in tempPhases:
        sites = []
        spec = {}
        fileSpecPhase = "./{}/species.in".format(phase)
        pf = open(fileSpecPhase,"r")
        rawsite = (pf.read().strip("\n")).split(" ")
        for s in range(0,len(rawsite)):
            asitespec = rawsite[s].split("=")
            site = asitespec[0]
            sites.append(site)
            spec.update({site: asitespec[1].split(",")})
            
        allPhase.update({phase:spec})
        pf.close()
        
        atermraw = []
        termsTemp = []
        termPerPhase = []
        fileTerm = "./{}/terms.in".format(phase)
        tf = open(fileTerm,"r")
        for line in tf:
            line = line.strip("\n")
            atermraw = line.split(":")
            for item in atermraw:
                aterm = item.split(",")
                termsTemp.append(aterm)
            termPerPhase.append(termsTemp)
            termsTemp = []
        terms.update({phase:termPerPhase})

        numOfparams = 0
    for phase in terms:
        for term in terms[phase]:
            for siteTerm in term:
                numOfparams = numOfparams + (len(elems)//int(siteTerm[0]))*(int(siteTerm[1])+1)
    return [elems, allPhase, terms, numOfparams]

[elems, allPhase, terms, numOfparams] = readInput()
