#!/usr/bin/env python
# coding: utf-8

# In[52]:


import random
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure


# In[53]:


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
        plt.plot(x,self.W)
        plt.grid()


# In[ ]:





# In[55]:


delta_t=T/m
r=0.05
sigma=0.2
sigma_J=0.5
mu_J=0
psi=1
N=15
T=5
m=150
x=np.linspace(0,T,m)

fig = plt.figure()
plt.rcParams["figure.figsize"] = (8,8)
fig, axis = plt.subplots(2)

for n in range(N):
    X=[np.log(100)] # will store X here
    S=[100]
    My_WP=Wiener_Process(T,m)
    for i in range(m-1):
        temp=X[-1]+(r-(sigma**2)/2 - psi*(np.exp((sigma_J**2)/2)-1))*delta_t
        J=random.gauss(mu_J,sigma_J)
        swith=np.random.poisson(psi*delta_t)
        temp=temp+sigma*(My_WP.W[i+1]-My_WP.W[i])+J*swith
        X.append(temp)
        temp=S[-1]*(1+(r - psi*(np.exp((sigma_J**2)/2)-1))*delta_t + sigma*(My_WP.W[i+1]-My_WP.W[i])+(np.exp(J)-1)*swith)
        S.append(temp)
    axis[0].plot(x, X, color='b')
    axis[1].plot(x, S, color='b')

axis[0].grid()
axis[1].grid()
plt.show()


# In[62]:


#Ok, I will try to add density 
t=np.linspace()
mu_bar=r-0.5*sigma**2 - psi*(np.exp((sigma_J**2)/2)-1)
mu=X[0]+t*mu_bar+mu_J*psi*t
var=t*sigma**2+psi*t*sigma_J**2
y=np.linspace(0,10,10000)
f_y=np.exp(-((y-mu)**2)/(2*var))/np.sqrt(2*math.pi*var)
plt.plot(y,f_y)
plt.grid()


# In[ ]:




