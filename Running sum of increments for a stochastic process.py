#!/usr/bin/env python
# coding: utf-8

# In[21]:


import random
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class Wiener_Process:
    def __init__(self,T,m):
        W=[0]
        S=0
        self.T=T
        self.m=m
        delta_t=T/m
        for i in range(1,m):
            inc=random.gauss(0,math.sqrt(delta_t))
            S=S+inc
            W.append(S)
        self.W=W
        
    def plot(self):
        x=np.linspace(0,self.T,self.m)
        fig = plt.figure()
        plt.plot(x,self.W)
        plt.grid()
        
class Geometric_Wiener_Process:
    def __init__(self,My_WP,mu,sigma,S0):
        T=My_WP.T
        m=My_WP.m
        S=np.ones([m,1])
        S[0]=S0
        delta_t=T/m
        for i in range(1,m):
            S[i]=S[i-1]*math.exp((mu-sigma*sigma/2)*delta_t)*math.exp(sigma*(My_WP.W[i]-My_WP.W[i-1]))
        self.S=S
        self.mu=mu
        self.sigma=sigma
        self.T=T
        self.m=m
    
    def plot(self):
        x=np.linspace(0,self.T,self.m)
        plt.plot(x,self.S,'b')
        plt.grid()

        
def running_sum_of_squared_increments(S):
    res=np.zeros([len(S)-1,1])
    res[0]=np.power((S[1]-S[0]),2)
    for i in range(1,len(S)-1):
        res[i]=res[i-1]+np.power((S[i+1]-S[i]),2)
    return res


# In[29]:


T=3
S0=0.7
delta_t=0.01
m=int(T/delta_t)
sigma=random.uniform(0.1,0.75)
mu=random.uniform(0.01,0.1)
print("Parameters: mu=",mu," sigma=",sigma)

for i in range (10):
    My_Wiener = Wiener_Process(T,m)
    My_GBM=Geometric_Wiener_Process(My_Wiener,mu,sigma,S0)
    My_GBM.plot()
    running_sum=running_sum_of_squared_increments(My_GBM.S)
    x=np.linspace(0,T,m-1)
    plt.plot(x,running_sum,'r')
    
plt.grid()
plt.show()


# In[32]:


df = pd.read_csv('HO.PA.csv')
S = df.Close #here we have closing prices of HO.PA stock
plt.plot(S)
plt.grid()

delta_t=1
running_sum=running_sum_of_squared_increments(S)
plt.plot(running_sum)


# In[ ]:




