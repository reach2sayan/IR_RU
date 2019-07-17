import numpy as np

def liquid_analysis():
    f0 = open('E0','r')
    fK = open('EK','r')
    s0 = f0.readlines()
    s1 = fK.readlines()
    print(s0)

liquid_analysis()
