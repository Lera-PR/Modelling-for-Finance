#!/usr/bin/env python
# coding: utf-8

# In[137]:


import random
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mpl_toolkits import mplot3d


# In[170]:


class function:
    def __init__(self,domain,rang):
        if(len(domain)!=len(rang)):
            print("Error: function is not bijective")
        else:
            self.domain=domain
            self.range=rang
            
    def plot(self):
        fig = plt.figure()
        plt.plot(self.domain,self.range)
        plt.grid()
        
    def derivative(self):
        der=[]
        temp=(self.range[1]-self.range[0])/(self.domain[1]-self.domain[0])
        der.append(temp)
        der.append(temp)
        for i in range(1,len(self.range)-1):
            temp=(self.range[i+1]-self.range[i])/(self.domain[i+1]-self.domain[i])
            der.append(temp)
        #print(self.range[1:len(self.range)],der)
        derivative=function(self.domain,der)
        return derivative
    
    def value_at_point(self,x):
        # if x is in the domain array, we are good. If not, we find two closest points, and take average
        if x<self.domain[0] or x>self.domain[len(self.domain)-1]:
            print("x is not in the domain")
            return -11111
        else:
            r=self.domain[0]
            i=0
            if(r==x):
                return self.range[0]
            else:
                while (r<x):
                    i=i+1
                    r=self.domain[i]
                return 0.5*(self.range[i]+self.range[i-1])
        
    def integral(self,l,r): #will calculate the integral of my function at the given interval (which has to be in the domain)
        if (l<self.domain[0]) or (l>self.domain[-1]):
            print("Left bound of the interval is outside of the domain")
            return -11111
        elif (r<self.domain[0]) or (r>self.domain[-1]):
            print("Right bound of the interval is outside of the domain")
            return -11111
        else: # the interval is inside the domain
            S=0
            i=0
            while (self.domain[i]<l):
                i=i+1
            #Found the beginning of the interval in the domain
            while (self.domain[i]<r):
                S=S+0.5*(self.range[i]+self.range[i+1])*(self.domain[i+1]-self.domain[i])
                i=i+1
            return S


# In[171]:


def sigma_imp(K):
    if(K>3):
        return 0.510-0.591*3+0.376*9-0.105*27+0.011*81
    else:
        return 0.510-0.591*K+0.376*K**2-0.105*K**3+0.011*K**4

def d_1(S0,r,K,T,t0):
    sigma=sigma_imp(K)
    res=(np.log(S0/K)+(r+0.5*sigma**2)*(T-t0))/(sigma*np.sqrt(T-t0))
    return res

def d_2(S0,r,K,T,t0):
    sigma=sigma_imp(K)
    res=(np.log(S0/K)+(r+0.5*sigma**2)*(T-t0))/(sigma*np.sqrt(T-t0))-sigma*np.sqrt(T-t0)
    return res

def F(d):
    return norm.cdf(d)

def V_c(S0,r,K,T,t0):
    d1=d_1(S0,r,K,T,t0)
    d2=d_2(S0,r,K,T,t0)
    V_c=S0*F(d1)-np.exp(-r*(T-t0))*K*F(d2)
    return V_c

def V_p(S0,r,K,T,t0):
    d1=d_1(S0,r,K,T,t0)
    d2=d_2(S0,r,K,T,t0)
    V_p=np.exp(-r*(T-t0))*K*F(-d2)-S0*F(-d1)
    return V_p
#Ok, so far so good; For different values K we can calculate the prices for simple call and put options,


# In[172]:


#To make life easier (and to practice) I will make a multiplication operation for two functions.
def multiplication(Func_1,Func_2):
    if(Func_1.domain[0]!=Func_2.domain[0]) or (Func_1.domain[len(Func_1.domain)-1]!=Func_2.domain[len(Func_2.domain)-1]):
        print("Cannot multiply these functions on the given domains")
        return -1
    else:
        new_domain=np.concatenate((Func_1.domain,Func_2.domain),axis=0)
        new_domain.sort()
        new_range=[]
        for x in new_domain:
            temp=Func_1.value_at_point(x)*Func_2.value_at_point(x)
            new_range.append(temp)
        Result=function(new_domain,new_range)
        return Result


# In[193]:


K=np.linspace(0,10,1000) #all possible values for K
S0=1
r=0
T=4
t0=0
#I will also need second derivatives of the payoff functions.
#(a) H(T,S(T))=max(S^2(T)-1.2S(T),0)
S1=np.linspace(0,S0-0.001,100) # possible prices of the stock at time T
S2=np.linspace(S0,100,10000)
S=np.concatenate((S1,S2),axis=0)
H=[]
for s in S1:
    temp=max(4-s**3,0)+max(s-2,0)
    H.append(temp)

for s in S2:
    temp=max(4-s**3,0)+max(s-2,0)
    H.append(temp)

Payoff_func=function(S,H)

First_der=Payoff_func.derivative()
Second_der=First_der.derivative() #OK, we calculated second derivative; Now we need the integral
I1=[]
I2=[]

for s in S1:
    temp=Second_der.value_at_point(s)*V_p(S0,r,s,T,t0)
    I1.append(temp)
    
for s in S2:
    temp=Second_der.value_at_point(s)*V_c(S0,r,s,T,t0)
    I2.append(temp)
    
Func_1=function(S1,I1)
Func_2=function(S2,I2)

res=Func_2.integral(S2[0],S2[-1])+Func_1.integral(S1[0],S1[-1])+Payoff_func.value_at_point(S0)
print(res)


# In[ ]:





# In[ ]:





# In[ ]:




