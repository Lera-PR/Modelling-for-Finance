#!/usr/bin/env python
# coding: utf-8

# In[177]:


import random
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import norm
from mpl_toolkits import mplot3d


# In[178]:


class function:
    def __init__(self,domain,rang):
        if(len(domain)!=len(rang)):
            print("Error: function is not bijective")
        else:
            self.domain=domain
            self.range=rang
            
    def plot(self):
        fig = plt.figure()
        plt.rcParams["figure.figsize"] = (8,5.5)
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


# In[209]:


def d_1(S0,r,K,T,t0,sigma):
    res=(np.log(S0/K)+(r+0.5*sigma**2)*(T-t0))/(sigma*np.sqrt(T-t0))
    return res

def d_2(S0,r,K,T,t0,sigma):
    res=(np.log(S0/K)+(r+0.5*sigma**2)*(T-t0))/(sigma*np.sqrt(T-t0))-sigma*np.sqrt(T-t0)
    return res

def F(d):
    return norm.cdf(d)

def V_c(S0,r,K,T,t0,sigma):
    d1=d_1(S0,r,K,T,t0,sigma)
    d2=d_2(S0,r,K,T,t0,sigma)
    V_c=S0*F(d1)-np.exp(-r*(T-t0))*K*F(d2)
    return V_c

def V_p(S0,r,K,T,t0,sigma):
    d1=d_1(S0,r,K,T,t0,sigma)
    d2=d_2(S0,r,K,T,t0,sigma)
    V_p=np.exp(-r*(T-t0))*K*F(-d2)-S0*F(-d1)
    return V_p


# In[221]:


def find_r(k,K_array):
    i=0
    r=K_array[0]
    while (r<k):
        i=i+1
        r=K_array[i]
    return i

def linear_interpolation(K_array,imp_vol_array):
    K=np.linspace(3,22,300)
    Y=[]
    for k in K:
        if(k<K_array[0]):
            Y.append(imp_vol_array[0])
        elif(k>K_array[-1]):
            Y.append(imp_vol_array[-1])
        else:
            r=find_r(k,K_array) #found indexes
            l=r-1
            
            a=(imp_vol_array[r]-imp_vol_array[l])/(K_array[r]-K_array[l]) #found coefficients for linear function
            b=imp_vol_array[r]-a*K_array[r]
            
            temp=a*k+b
            Y.append(temp)
    Lin_fun=function(K,Y)
    return Lin_fun

def NN_interpolation(K_array,imp_vol_array):
    K=np.linspace(3,22,100)
    Y=[]
    for k in K:
        if(k<K_array[0]):
            Y.append(imp_vol_array[0])
        elif(k>K_array[-1]):
            Y.append(imp_vol_array[-1])
        else:
            r=find_r(k,K_array) #found indexes
            l=r-1
            diff_r=K_array[r]-k
            diff_l=k-K_array[l]
            
            if(diff_r>diff_l):
                Y.append(imp_vol_array[l])
            else:
                Y.append(imp_vol_array[r])
    fun=function(K,Y)
    return fun

def Cubic_spline_interpolation(K_array,imp_vol_array):
    a=imp_vol_array
    h=[]
    for i in range(len(a)-1):
        h.append(K_array[i+1]-K_array[i])
        
    M=np.zeros([len(a),len(a)])
    for i in range(len(a)):
        if i==0:
            M[i][i]=1
        elif i==len(a)-1:
            M[i][i]=1
        else:
            M[i][i]=2*(h[i-1]+h[i])
            M[i][i+1]=h[i]
            M[i][i-1]=h[i-1]
    r=[]
    r.append(0)
    for i in range(1,len(a)-1):
        temp=3*(((a[i+1])-a[i])/h[i] - (a[i]-a[i-1])/h[i-1])
        r.append(temp)
    r.append(0)
    M=np.linalg.inv(M)
    
    c=M.dot(r)
    c=c.tolist()
    
    b=[]
    for i in range(1,len(c)):
        temp=(a[i]-a[i-1])/h[i-1] +(2*c[i]+c[i-1])*h[i-1]/3
        b.append(temp)
    
    d=[]
    for i in range(1,len(c)):
        temp=(c[i]-c[i-1])/(3*h[i-1])
        d.append(temp)
    
    #OK, we found all the coefficients for the cubic splines
    K=np.linspace(3,22,300)
    Y=[]
    for k in K:
        if(k<K_array[0]):
            Y.append(a[0])
        elif(k>K_array[-1]):
            Y.append(a[-1])
        else:
            i=find_r(k,K_array) #found indexes
            Y.append(a[i]+b[i-1]*(k-K_array[i])+c[i]*(k-K_array[i])**2+d[i-1]*(k-K_array[i])**3)
    fun=function(K,Y)
    return fun


# In[222]:


def check_call_butterfly_arbitrage(V_c_function):
    check=0
    for f in V_c_function.range:
        if f<-0.001:
            print(f)
            return 1
    return 0


# In[223]:


K_array=[3.28, 5.46, 8.2, 10.93, 13.66, 16.39, 19.12, 21.86]
imp_vol_array=[0.3137, 0.2249, 0.1491, 0.0909, 0.0685, 0.0809, 0.0945, 0.1063]
Lin_fun_vol=linear_interpolation(K_array,imp_vol_array)

NN_fun_vol=NN_interpolation(K_array,imp_vol_array)

Cub_fun_vol=Cubic_spline_interpolation(K_array,imp_vol_array)


T=1
t0=0
S0=10.5
r=0.04
K=np.linspace(3,22,100)
Y=[]
for k in K:
    imp_sigma=NN_fun_vol.value_at_point(k)
    temp=V_c(S0,r,k,T,t0,imp_sigma)
    Y.append(temp)

New_func=function(K,Y)
First_der=New_func.derivative()
density=First_der.derivative()
for i in range(0,len(K)):
    density.range[i]=density.range[i]*np.exp(r*(T-t0))
density.plot()

check=check_call_butterfly_arbitrage(density)
if(check==0):
    print("No butterfly arbitrage")
else:
    print("There is butterfly arbitrage")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




