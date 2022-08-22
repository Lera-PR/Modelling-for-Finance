import random
import math
import matplotlib.pyplot as plt

class Wiener_Process:
    def __init__(self,T,m):
        W=[0]
        S=0
        delta_t=T/m
        for i in range(0,m):
            inc=random.gauss(0,math.sqrt(delta_t))
            S=S+inc
            W.append(S)
        self.W=W
        
    def plot(self):
        fig = plt.figure()
        plt.plot(self.W)
        plt.grid()


T=1
m=200
My_Wiener=Wiener_Process(T,m)
My_Wiener.plot()

Int=[0] #here I calculate integral W(s)ds on (0,T)
delta_t=T/m
for i in range(1,m):
    temp=Int[i-1]+My_Wiener.W[i]*delta_t
    Int.append(temp)

x=np.linspace(0,T,m)
plt.plot(x,Int)


Int=[0] #here I calculate integral W(s)dW(s) on (0,T)
for i in range(1,m-1):
    temp=0.5*My_Wiener.W[i]**2 - 0.5*i*delta_t
    Int.append(temp)
Int.append(temp)

plt.plot(x,Int)
plt.legend(['My Wiener', 'First integral','Second integral'])
plt.show()
