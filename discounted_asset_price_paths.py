#!/usr/bin/env python
# coding: utf-8

# In[214]:


import random
import math
import matplotlib.pyplot as plt
import numpy as np


# In[215]:


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
    def __init__(self,My_WP,mu,sigma):
        T=My_WP.T
        m=My_WP.m
        S=np.ones([m,1])
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
        
class Saving_Account:
    def __init__(self,r,T,m):
        self.T=T
        self.r=r
        self.m=m
        delta_t=T/m
        M=np.ones([m,1])
        for i in range(1,m):
            M[i]=M[i-1]*math.exp(r*delta_t)
        self.M=M
        
    def plot(self):
        x=np.linspace(0,self.T,self.m)
        plt.plot(x,self.M,'r')
        plt.grid()
        


# In[240]:


T=10
m=1000
mu=0.15
sigma=0.1

Average_GM=np.zeros([m,1])

n=1000 #In this loop we generate 1000 paths of GBM and plot 10 of them. We do this to find average
for i in range(0,n):
    My_Wiener=Wiener_Process(T,m)
    My_Geometric_Wiener=Geometric_Wiener_Process(My_Wiener,mu,sigma)
    Average_GM=np.add(Average_GM,My_Geometric_Wiener.S)
    if(i<10):
        My_Geometric_Wiener.plot()
    
for i in range(0,m):
    Average_GM[i]=Average_GM[i]/n
x=np.linspace(0,T,m)
plt.plot(x,Average_GM,'r')
plt.grid()

plt.show()


# In[244]:


r=0.05
My_Saving_Account=Saving_Account(r,T,m)

Average_S_over_M=np.zeros([m,1])

n=100 #In this loop we generate 1000 paths of GBM and plot 10 of them. We do this to find average
for i in range(0,n):
    My_Wiener=Wiener_Process(T,m)
    My_Geometric_Wiener=Geometric_Wiener_Process(My_Wiener,mu,sigma)
    My_S_over_M=np.divide(My_Geometric_Wiener.S,My_Saving_Account.M)
    Average_S_over_M=np.add(Average_S_over_M,My_S_over_M)
    if(i<10):
        plt.plot(x, My_S_over_M,'b')
        
for i in range(0,m):
    Average_S_over_M[i]=Average_S_over_M[i]/n
    
plt.plot(x,Average_S_over_M,'r')
    
plt.grid()
plt.show()


# In[243]:


#now let's try risk-neutral measure. Instead of the drift mu we use interest rate r.
Average_S_over_M=np.zeros([m,1])

n=100 #In this loop we generate 1000 paths of GBM and plot 10 of them. We do this to find average
for i in range(0,n):
    My_Wiener=Wiener_Process(T,m)
    My_Geometric_Wiener=Geometric_Wiener_Process(My_Wiener,r,sigma)
    My_S_over_M=np.divide(My_Geometric_Wiener.S,My_Saving_Account.M)
    Average_S_over_M=np.add(Average_S_over_M,My_S_over_M)
    if(i<10):
        plt.plot(x, My_S_over_M,'b')
        
for i in range(0,m):
    Average_S_over_M[i]=Average_S_over_M[i]/n
    
plt.plot(x,Average_S_over_M,'r')
    
plt.grid()
plt.show()


# In[ ]:




