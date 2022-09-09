#!/usr/bin/env python
# coding: utf-8

# In[17]:


import random
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mpl_toolkits import mplot3d


# In[18]:


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
        
        


# In[19]:


def F(d):
    return norm.cdf(d) #this is a cumulitive function of a standard normal distribution

def Delta_of_European_call_option(t,S,K,m,r,sigma,T):
    d1=(np.log(S/K)+(T-t)*(r+0.5*sigma**2))/(sigma*np.sqrt(T-t))
    return F(d1)

def Gamma_of_European_call_option(t,S,K,m,r,sigma,T):
    d1=(np.log(S/K)+(T-t)*(r+0.5*sigma**2))/(sigma*np.sqrt(T-t))
    Gamma=np.exp(-d1**2/2)/(S*sigma*np.sqrt(T-t)*np.sqrt(2*math.pi))
    return Gamma

def Value_of_European_call_option(t,S,K,m,r,sigma,T):
    d1=(np.log(S/K)+(T-t)*(r+0.5*sigma**2))/(sigma*np.sqrt(T-t))
    d2=d1-sigma*np.sqrt(T-t)
    V=S*F(d1)-K*np.exp(-r*(T-t))*F(d2)
    return V


# In[23]:


T=1
m=1000
d=T/m #step
sigma=0.2
K=0.95 #strike price
S0=1
r=0.1
My_W=Wiener_Process(T,m)
My_GBM=Geometric_Wiener_Process(My_W,r,sigma,S0)
My_GBM.plot()

V=np.zeros([m,1])
t=np.linspace(0,T,m)
for i in range(0,m):
    V[i]=Value_of_European_call_option(i*d,My_GBM.S[i],K,m,r,sigma,T)
plt.plot(t,V,'m')

Delta=np.zeros([m,1])
for i in range(0,m):
    Delta[i]=Delta_of_European_call_option(i*d,My_GBM.S[i],K,m,r,sigma,T)
plt.plot(t,Delta,'g')

PnL=np.zeros([m,1])
PnL[0]=V[0]-Delta[0]*My_GBM.S[0]
for i in range(1,m-1):
    PnL[i]=PnL[i-1]*np.exp(r*d)-My_GBM.S[i]*(Delta[i]-Delta[i-1])
PnL[m-1]=PnL[m-2]*np.exp(r*d)-max(My_GBM.S[m-1]-K,0)+Delta[m-2]*My_GBM.S[m-1]
plt.plot(t,PnL,'r')


# In[ ]:





# In[ ]:




