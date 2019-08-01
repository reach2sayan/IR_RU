import numpy as np
import os
import re
import shutil
from itertools import product,chain,combinations,permutations

elems = [] #holds all elements
sites = []
spec = {} #holds species hashed by site 
phases = [] #holds all phases calculated
terms = {}
allPhase = {} #holds all sites for all phase
numOfparams = 0 #records no. of params to fit with OC
valueList = ('energy','svib_ht')
phaseList = os.listdir("/users/ssamanta/atat/data/sqsdb/")

#===================================================================
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
        pf.close()
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
        termdata = tf.readlines()
        tf.close()
        #loop through the line to calculate level and order of split, refer to compound energy formalism
        for line in termdata:
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
#====================================================================

#records the SGTE structure and physical constants for each element
def elemData(elems):
    elem_data = {}
    elemsf = open("/users/ssamanta/atat/data/sgte_elements.tdb","r")
    elemdata = elemsf.readlines()
    elemsf.close()
    for line in elemdata:
        for element in elems:
            if "ELEM_"+element.upper() in line:
                tempstr = "ENTER ELEMENT "+element.upper()+" "+element+" "+line.strip("ELEMENT ELEM_"+element.upper())
                tempstr = tempstr.strip("\n").strip("!")
                tempstr = tempstr[:-1]
                elem_data.update({element:tempstr})
    return elem_data
                
#records phase structure for each phases from the database
def phaseData(allPhase):
    phaseData = {} #stores phase structure in OC format per phase
    #iterate through all the phases for individual sites and species
    for name,val in allPhase.items():
        stringHead = "ENTER PHASE "+name #line 1: eg. enter phase FCC_A1
        stringDat = "CEF "+str(len(val))+" " #record model and total no. of sites, eg. CEF 2.0 (for HCP_A3)
        fphaseSpec = open("./"+name+"/species.in","r")
        fphaseMult = open("./"+name+"/mult.in","r")
        #records stoichometric multiplicity for each sublattice
        smult = fphaseMult.readlines()
        fphaseMult.close()
        smult = smult[0] 
        smult = smult.strip('\n')
        smult = smult.split("\t") #sublattice multiplicity separated by tabs
        smult = smult[:-1] #remove stray spaces
        #records species in each subblattice for the given phase
        sspec = fphaseSpec.readlines()
        fphaseSpec.close()
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
#==============================================================

#records the SGTE, ABIN and REFERENCE functions        
def tpfunc(elems,allPhase,terms,numOfparams):
    sgtefunc = [] #contains the SGTE functions
    abinfunc = [] #contains the ab-initio functions (to be fit)
    reffunc = [] #contains the reference functions
    
    funcHead = "ENTER TP " #every function begins like this
    fsgte = open("/users/ssamanta/atat/data/sgte_freee.tdb","r")
    sgte_freee = fsgte.readlines()
    fsgte.close()
    for line in sgte_freee:
        for element in elems:
            for phase in allPhase:
                if "SGTE_"+phase+"_ELEM_"+element.upper() in line:
                    #sindex = line.find(">")
                    line = re.sub(r"<NL>","\n",line)
                    line = line[:-1]+ " REFDUM"
                    line = re.sub(r"FUNCTION ","",line)
                    line = re.sub(r"_ELEM","",line)
                    line = re.sub(r"   298.15", "FUN 298.15",line)
                    line = funcHead + line
                    line = re.sub(r"! ","",line)
                    sgtefunc.append(line)
    paramCounter = 0
    for phase in allPhase:
        phaseDict = allPhase[phase]
        for v in product(*phaseDict.values()):
            line = ""
            elem = "".join(v)
            elem = elem.upper()
            line = line + funcHead + "ABIN_" + phase +"_"+ elem + " FUN 298.15 A" + str(paramCounter) + "+A" + str(paramCounter+1) + "*T; 10000 N REFDUM"
            paramCounter = paramCounter + 2
            abinfunc.append(line)

    for phase in allPhase:
        for element in elems:
            line = ""
            line = line + funcHead + "REF_" + phase + "_" + element.upper() + " FUN 298.15 ABIN_" + availPhaseElem(phase,element.upper())[0] + "_" + element.upper() + " - SGTE_" + availPhaseElem(phase,element.upper())[0] + "_" + element.upper() + "; 10000 N REFDUM"
            reffunc.append(line)
    return [sgtefunc,abinfunc,reffunc,paramCounter]
#==================================================================

#checks if the phase exists in SGTE, if not then gets the stable element
def availPhaseElem(phase,elem):
    #checks in the SGTE database, if found, then cool
    fstableElem = open("/users/ssamanta/atat/data/sgte_freee.tdb","r")
    stableElem = fstableElem.readlines()
    fstableElem.close()
    for line in stableElem:
        if phase+"_ELEM_"+elem in line:
            return [phase,elem]
    
    #if not found then replace the phase with the phase that is stable as per the element
    stableElem = []
    fstableElem = open("/users/ssamanta/atat/data/sgte_elements.tdb","r")
    stableElem = fstableElem.readlines()
    for line in stableElem:
        if "ELEM_"+elem.upper() in line:
            line = line.split(" ")
            return [line[2],elem]
#===================================================================

def paramEndmem(allPhase,elems):
    
    param_data = []
    paramHead = "ENTER PARAM "
    for phase in allPhase:
        phaseDict = allPhase[phase]
        for v in product(*phaseDict.values()):
            paramline = "G("+phase+","+":".join(v).upper()+";0) 298.15 "
            fphaseMult = open("./"+phase+"/mult.in","r")
            #records stoichometric multiplicity for each sublattice
            smult = fphaseMult.readlines()
            fphaseMult.close()
            smult = smult[0] 
            smult = smult.strip('\n')
            smult = smult.split("\t") #sublattice multiplicity separated by tabs
            smult = smult[:-1]
            mult = []
            for item in smult:
                tempmult = item.split("=")
                mult.append(tempmult[1])

            #following the strategy from sqs2tdb
            paramline = paramline + "ABIN_"+phase+"_"+"".join(v).upper() + " - "

            for item in zip(v,mult):
                paramline = paramline + item[1]+"*REF_"+phase.upper()+"_"+item[0].upper()+" - "
                
            paramline = paramline[:-3] #removing trailing -ve sign and space
            paramline = paramline+"; 10000 N"
            param_data.append(paramHead+paramline)
            
    return param_data
#=================================================================

def paramMixing(terms,elems,allPhase,paramCounter):
    
    param_mixing_data = []
    params = {}
    for phase,term in terms.items():
        for oneterm in term:
            elemComb = list(chain.from_iterable(combinations(elems,r) for r in range(len(term)+1)))[1:]
            elemdum = [",".join(e).upper() for e in elemComb]
            elemComb = elemdum
            #print(phase+" "+str(elemComb))
            lev = 0
            for sublats in oneterm:
                #print(phase+" "+str(sublats))
                if int(sublats[1]) >= lev:
                    lev = int(sublats[1])
                perm = list(permutations(elemComb,len(oneterm)))
                permdum = [":".join(e) for e in perm]
                perm = permdum
        params.update({phase:[perm, lev]})
    for phase,key in params.items():
        for item in key[0]:
            oneitem = item.split(":")
            #print(phase+" "+str(oneitem))
            for pieces in oneitem:
                if "," in pieces:
                    for i in range(key[1]+1):
                        line = "L("+phase+","+":".join(oneitem)+";"+str(i)+")"
                        line = line + " A"+str(paramCounter) + "+A"+str(paramCounter+1)+"*T"+"; 10000 N REFDUM"
                        paramCounter = paramCounter+2
                        param_mixing_data.append("ENTER PARAM "+line)

    return [param_mixing_data, paramCounter]
#================================================================

#records the different experimental enthalpy values
def addExpEnthalpy(allPhase,expCounter):
    
    experiments_E = {} #record all enthalpies hashed by phase
    for phase in allPhase:
        exp = []
        filePhaseConc = open("./{}/conc_energy.out".format(phase),"r")
        filePhaseEnergy = open("./{}/y_energy.in".format(phase),"r")
        conc = filePhaseConc.readlines()
        energy = filePhaseEnergy.readlines()
        filePhaseConc.close()
        filePhaseEnergy.close()
        for c,e in zip(conc,energy):
            c = c.strip("\n")
            c = c.replace("\t"," ")
            e = e.strip("\n")
            st = "EX"+str(expCounter)+" "+str(c)+" "+e
            expCounter = expCounter+1
            exp.append(st)
        if len(exp) == 0:
            continue
        experiments_E.update({phase:exp})
    
    return [experiments_E, expCounter]
#=================================================================

#records the different experimental entropy values
def addExpEntropy(allPhase,expCounter):
    experiments_S = {} #record all entropy values
    for phase in allPhase:
        exp = []
        filePhaseConc = open("./{}/conc_svib_ht.out".format(phase),"r")
        filePhaseEnergy = open("./{}/y_svib_ht.in".format(phase),"r")
        conc = filePhaseConc.readlines()
        energy = filePhaseEnergy.readlines()
        filePhaseConc.close()
        filePhaseEnergy.close()
        for c,e in zip(conc,energy):
            c = c.strip("\n")
            c = c.replace("\t"," ")
            e = e.strip("\n")
            st = "EX"+str(expCounter)+" "+str(c)+" "+e
            expCounter = expCounter+1
            exp.append(st)
        if len(exp) == 0:
            continue
        experiments_S.update({phase:exp})
    
    return [experiments_S, expCounter]
#==================================================================

#writes to the asses.OCM file
def writeToFile(elem_data, phase_data, param_data, param_mixing_data, sgtefunc, abinfunc, reffunc, experiments_E, experiments_S,elems,paramCounter,expCounter):
    fasses = open("asses.OCM","w+")
    fasses.write("NEW Y\n\n")
    fasses.write("ENTER OPT_COEF {}\n\n".format(paramCounter))

    #writes physical parameter of the elements
    for line in elem_data.values():
        fasses.write(line+"\n")

    #writes the phase, structures and model (CEF used here)
    fasses.write("\n")
    for line in phase_data.values():
        fasses.write(line[0]+"\n")
        fasses.write(line[1]+"\n")

    #writes the SGTE functions from database
    fasses.write("\n")
    for line in sgtefunc:
        fasses.write(line+"\n")
    
    #writes the Ab-initio functions (to be fit)
    fasses.write("\n")
    for line in abinfunc:
        fasses.write(line+"\n")

    #writes the reference functions (in sqs2tdb strategy)
    fasses.write("\n")
    for line in reffunc:
        fasses.write(line+"\n")

    #writes end-member parameters
    fasses.write("\n")
    for line in param_data:
        fasses.write(line+"\n")
        
    #writes mixing (L) parameters
    fasses.write("\n")
    for line in param_mixing_data:
        fasses.write(line+"\n")

    #writes enthalpy experiments
    fasses.write("\n")
    for phase,exp in experiments_E.items():
        fasses.write("\n")
        fasses.write("ENTER MANY_EQUILIBRIA")
        compcond = []
        colcount = 1
        for element in elems:
            compcond.append("X("+element+")=@"+str(colcount))
            colcount = colcount+1
        fasses.write("\n")
        fasses.write("CONDITION P=1e5 T=1000 N=1 "+" ".join(compcond)+"\n")
        fasses.write("EXPERIMENT H("+phase+")=@"+str(colcount)+":0.1"+"\n")
        fasses.write("TABLE_START\n")
        for line in exp:
            fasses.write(line+"\n")
        fasses.write("TABLE_END\n")

    #writes entropy experiments
    fasses.write("\n")
    for phase,exp in experiments_S.items():
        fasses.write("\n")
        fasses.write("ENTER MANY_EQUILIBRIA")
        compcond = []
        colcount = 1
        for element in elems:
            compcond.append("X("+element+")=@"+str(colcount))
            colcount = colcount+1
        fasses.write("\n")
        fasses.write("CONDITION P=1e5 T=1000 N=1 "+" ".join(compcond)+"\n")
        fasses.write("EXPERIMENT S("+phase+")=@"+str(colcount)+":0.1"+"\n")
        fasses.write("TABLE_START\n")
        for line in exp:
            fasses.write(line+"\n")
        fasses.write("TABLE_END\n")

    #write range of experiments to consider
    fasses.write("\n")
    fasses.write("SET RANGE_EXP 2 {}".format(expCounter))

    #write initial guess for params
    fasses.write("\n")
    for i in range(0,paramCounter):
        fasses.write("SET VAR {} 100\n".format(i))

    fasses.write("\nOPT 1000\nSET INTER")
    fasses.close()
#===================================================================

#shortens function names for OC compatibility
def shortenFun():
    
    #we first record each function name and assign a unique short name
    fileMacro = open("asses.OCM","r")
    macro = fileMacro.readlines()
    functions = {}
    funccount = 0
    for line in macro:
        if "ENTER TP" in line:
            lines = line.split(" ")
            functions.update({lines[2]:"FUNC{}".format(funccount)})
            funccount = funccount + 1
    macro.clear()

    #now re-read the file, and replace all old functions with new names
    shortFileMacro = open("asses.OCM.tmp","w+")
    fileMacro.seek(0,0)
    macro = fileMacro.readlines()
    for line in macro:
        for oldfunc,newfunc in functions.items():
            if oldfunc in line:
                line = line.replace(oldfunc,newfunc)
        shortFileMacro.write(line)

    shortFileMacro.close()
    fileMacro.close()
#===================================================================

[elems, allPhase, terms, numOfparams] = readInput()
elem_data = elemData(elems)
phase_data = phaseData(allPhase)
[sgtefunc, abinfunc, reffunc, paramCounter] = tpfunc(elems,allPhase,terms,numOfparams)

param_data = paramEndmem(allPhase,elems)
[param_mixing_data, paramCounter] = paramMixing(terms,elems,allPhase,paramCounter)
expCounter = 0
[experiments_E,expCounter] = addExpEnthalpy(allPhase,expCounter)
[experiments_S,expCounter] = addExpEntropy(allPhase,expCounter)
if os.path.exists("asses.OCM"):
    os.remove("asses.OCM")
writeToFile(elem_data, phase_data, param_data, param_mixing_data, sgtefunc, abinfunc, reffunc,experiments_E,experiments_S,elems,paramCounter,expCounter)

if os.path.exists("asses.OCM.tmp"):
    os.remove("asses.OCM.tmp")

shortenFun()
os.remove("asses.OCM")
shutil.move(os.path.join(".","asses.OCM.tmp"),os.path.join(".","asses.OCM"))

