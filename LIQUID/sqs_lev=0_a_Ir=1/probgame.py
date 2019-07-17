import numpy as np
import matplotlib.pyplot as plt

def singlegame():
    count = 0
    add = 0
    while add < 1:
        add = add + np.random.randint(0,1000)/1000
        count = count + 1
    return count

X = list()
for i in range(1,100000):
    X.append(singlegame())

A = [np.mean(X), np.var(X)]
print(A)
plt.hist(X,bins='auto')
plt.show()
