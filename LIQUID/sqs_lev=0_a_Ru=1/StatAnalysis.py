import matplotlib.pyplot as plt
import numpy as np

def findBlockAverage():
    
    press = []
    fpress = open('press','r')
    p = fpress.readlines()
    fpress.close()
    for item in p:
        item.strip('\n')
        press.append(float(item))
    
    time = [i*3 for i in range(len(press))]
    press = press[1000:]
    t= len(press)
    tBs = []
    p = []
    stds = []

    for i in range(1,int(t/5)):
        tB = i
        nB = int(t/tB)

        tBs.append(tB)
        stds.append(np.std(press[0:t:tB]))
    
    #print(tBs)
    #print(stds)
    plt.plot(tBs,stds)
    plt.show()
    
    tB = 25 
    sterr = np.std(press[0:t:tB])/np.sqrt(len(press[0:t:tB]))
    val = np.mean(press[0:t:tB])

    print(val)
    print(sterr)
    
    fresult = open("avg_press.out","w+")
    fresult.write("val="+str(val)+" "+"dev="+str(sterr)+"\n")
    fresult.close()
press = findBlockAverage()

