import random
import math
import matplotlib.pyplot as plt
import numpy as np

class Poisson_Process:
    def __init__(self,T,m,psi):
        P=[0]
        S=0
        delta_t=T/m
        for i in range(m-1):
            inc=np.random.poisson(psi*delta_t)
            S=S+inc
            P.append(S)
        self.P=P
        self.T=T
        self.m=m
        self.psi=psi
        
    def plot(self):
        x=np.linspace(0,self.T,self.m)
        plt.plot(x,self.P,color='b')
        plt.grid()

        
class Compensated_Poisson_Process:
    def __init__(self,My_PP):
        self.T=My_PP.T
        self.m=My_PP.m
        self.psi=My_PP.psi
        delta_t=self.T/self.m
        P=[0]
        for i in range(1,m):
            temp=My_PP.P[i]-self.psi*delta_t*i
            P.append(temp)
        self.P=P
    
    def plot(self):
        x=np.linspace(0,self.T,self.m)
        plt.plot(x,self.P)
        plt.grid()


T=30
m=1000
psi=1
N=20 # will plot 20 Poisson process paths
figure, axis = plt.subplots(2)
x=np.linspace(0,T,m)
for i in range(N):
    My_Poisson=Poisson_Process(T,m,psi)
    axis[0].plot(x, My_Poisson.P, color='b')
    My_Compensated_PP=Compensated_Poisson_Process(My_Poisson)
    axis[1].plot(x,My_Compensated_PP.P, color='r')
axis[0].grid()
axis[1].grid()
plt.show()
