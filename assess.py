import numpy as np
import os
import re

elems = [] #holds all elements
sites = []
spec = {} #holds species hashed by site 
phases = [] #holds all phases calculated
terms = {}
allPhase = {} #holds all sites for all phase
numOfparams = 0 #records no. of params to fit with OC
valueList = ('energy','svib_ht')
phaseList = os.listdir("/users/ssamanta/atat/data/sqsdb/")

#records information about phases, sublattice, species and levels of fit
def readInput():
    sf = open("species.in","r") #read general species.in 
    
    elems = (sf.read().strip("\n")).split(",") #remove the newline and extract elements separated by comma(,)
    sf.close()

    #read the folders which are named after the phases under consideration
    files = os.listdir(".")
    tempPhases = [file for file in files if file in phaseList]
    
    #iterate through the phases
    for phase in tempPhases:
        sites = [] #holds the different sites 
        spec = {} #hold the species hashed by sites
        fileSpecPhase = "./{}/species.in".format(phase) #opens the corresponding phase species.in file
        pf = open(fileSpecPhase,"r")
        rawsite = (pf.read().strip("\n")).split(" ") #remove newline and separate sites (by space)
        #iterate through sites
        for s in range(0,len(rawsite)):
            asitespec = rawsite[s].split("=") #get elements at each site
            site = asitespec[0]
            sites.append(site) #add to sites
            spec.update({site: asitespec[1].split(",")}) #update species at each site
            
        allPhase.update({phase:spec}) #update each phase with site and species info
        pf.close()
        
        atermraw = []
        termsTemp = []
        termPerPhase = [] #holds the terms to consider in the fit per phase
        fileTerm = "./{}/terms.in".format(phase) #read the terms.in file for each phase
        tf = open(fileTerm,"r")
        #loop through the line to calculate level and order of split, refer to compound energy formalism
        for line in tf:
            line = line.strip("\n") #remove newline
            atermraw = line.split(":") #remove order per sublattice
            for item in atermraw:
                aterm = item.split(",") #separate order and level for each sublattic
                termsTemp.append(aterm)
            termPerPhase.append(termsTemp)
            termsTemp = []
        terms.update({phase:termPerPhase}) #collect all terms for all levels for all phases

        numOfparams = 0
    
    #count the required no. of params to be fit by OC
    for phase in terms:
        for term in terms[phase]:
            for siteTerm in term:
                numOfparams = numOfparams + (len(elems)//int(siteTerm[0]))*(int(siteTerm[1])+1)
    return [elems, allPhase, terms, numOfparams]

#records the SGTE structure and physical constants for each element
def elemData(elems):
    elem_data = {}
    elemsf = open("/users/ssamanta/atat/data/sgte_elements.tdb","r")
    elemdata = elemsf.readlines()
    for line in elemdata:
        for element in elems:
            if "ELEM_"+element.upper() in line:
                tempstr = "enter element "+element.upper()+" "+element+" "+line.strip("ELEMENT ELEM_"+element.upper())
                tempstr = tempstr.strip("\n").strip("!")
                elem_data.update({element:tempstr})
    return elem_data
                
#records phase strcuture for each phases from the database
def phaseData(allPhase):
    phaseData = {} #stores phase structure in OC format per phase
    #iterate through all the phases for individual sites and species
    for name,val in allPhase.items():
        stringHead = "enter phase "+name #line 1: eg. enter phase FCC_A1
        stringDat = "CEF "+str(len(val))+" " #record model and total no. of sites, eg. CEF 2.0 (for HCP_A3)
        fphaseSpec = open("./"+name+"/species.in","r")
        fphaseMult = open("./"+name+"/mult.in","r")
        #records stoichometric multiplicity for each sublattice
        smult = fphaseMult.readlines()
        smult = smult[0] 
        smult = smult.strip('\n')
        smult = smult.split("\t") #sublattice multiplicity separated by tabs
        smult = smult[:-1] #remove stray spaces
        #records species in each subblattice for the given phase
        sspec = fphaseSpec.readlines()
        sspec = sspec[0]
        sspec = sspec.strip('\n')
        sspec = sspec.split(" ")
        strdet = ""
        #iterates over sublattice to record the element with its multiplicity in OC format
        for st,val in zip(sspec,smult):
            ssublat = st.split('=')
            smultlat = val.split('=')
            strdet = strdet + str(float(smultlat[1])) + " "
            for s in ssublat[1].split(","):
                strdet = strdet+ s.upper()+","
            strdet = strdet[:-1]  #remove trailing comma
            strdet = strdet + "; " #adds semi-colon separating sublattices
        stringDat = stringDat + strdet
        stringDat = stringDat[:-1] #removes the trailing semi-colon
        phaseData.update({name:[stringHead,stringDat]})
    return phaseData
        
def tpfunc(elems,allphase,terms, numOfparams):
    sgtefunc = {}
    abinfunc = {}
    reffunc = {}

    

[elems, allPhase, terms, numOfparams] = readInput()
elem_data = elemData(elems)
phase_data = phaseData(allPhase)
