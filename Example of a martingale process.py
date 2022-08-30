#!/usr/bin/env python
# coding: utf-8

# In[30]:


import numpy as np
import random
import matplotlib.pyplot as plt
import math


# In[48]:


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
        plt.plot(x,self.S)
        plt.grid()


# In[58]:


#We need two GBM process
mu_x=0.04
sigma_x=0.12 # this can vary

mu_y=0.0625-0.15*sigma_x
if(mu_y<0):
    print("beta is negative")
sigma_y=0.15

T=2
m=200
My_Wiener=Wiener_Process(T,m)
X=Geometric_Wiener_Process(My_Wiener,mu_x,sigma_x)
Y=Geometric_Wiener_Process(My_Wiener,mu_y,sigma_y)
X.plot()
Y.plot()

lambd=0

dt=T/m
Z=np.zeros([m,1])
for i in range(m):
    Z[i]=2*X.S[i]/Y.S[i] - lambd*dt*i
x=np.linspace(0,T,m)
plt.plot(x,Z)
plt.grid()
plt.show()


# In[ ]:




