import numpy as np
import matplotlib.pyplot as plt
import math

def liquid_analysis():
    #read the files
    f0 = open('E0','r')
    fK = open('EK','r')
    s0 = f0.readlines()
    s1 = fK.readlines()

    #remove linebreak and convert to float
    s0 = [s.strip('\n') for s in s0]
    s0 = [float(i) for i in s0]
    s1 = [s.strip('\n') for s in s1]
    s1 = [float(i) for i in s1]

    s = np.add(s0,s1)
    #print(s)
    #each timestep in MD is 1.5 fs

    #print(time)
    time = [i*1.5 for i in range(0,len(s))]
    plt.plot(time,s)
    #plt.show()
    
    #a shall store the values in time series from equlibriation
    a = s[500:]/54
    #print(a)
    t = len(a) #stores the total no. of points in analysis

    print('loading completed')

    tBs = list()
    p = list()
    stdevs = list()
    rho = list()

    plt.plot(a)
    #plt.show()
    
    for i in range(1,int(t/10)):
        
        tB = i #try out different block sizes for block average
        nB = int(t/tB) #extract integer part for no. of blocks

        A_b = np.zeros((nB-1,1))

        for j in range(0,nB-1):
            A_b[j,0] = np.mean(a[j*tB:j*tB+1])

        #print(A_b)
        tBs.append(tB)
        stdevs.append(np.std(a[0:tB:t]))
        b = np.array(a[0:t:tB])
        #print(a[0:t:tB])
        #input()
        #print(len(stdevs))
        covar = np.cov(b[0:len(b)-1],b[1:len(b)])
        #print((covar))
        rho.append(covar[0,1]/covar[0,0])
        p.append(tB*np.var(A_b)/np.var(a))
        print(type(tBs))

    plt.plot(np.array(tBs), np.array(p))
    plt.xlabel('block time')
    plt.ylabel('decay time')
    plt.legend()
    plt.show()
        
liquid_analysis()
