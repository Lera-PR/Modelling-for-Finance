import random
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mpl_toolkits import mplot3d


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
        
        
def F(d):
    return norm.cdf(d) #this is a cumulitive function of a standard normal distribution

#next two functions simply calculate value of a cash-or-nothing call option and its delta
def Delta_of_cash_or_nothing_call_option(t,S,K,m,r,sigma,T,A):
    d2=(np.log(S/K)+(T-t)*(r-0.5*sigma**2))/(sigma*np.sqrt(T-t))
    delta=(A/S)*(np.exp(-r*(T-t)))/(sigma*np.sqrt(T-t))*np.exp(-d2**2/2)
    return delta

def Value_of_cash_or_nothing_call_option(t,S,K,m,r,sigma,T,A):
    d2=(np.log(S/K)+(r-0.5*sigma**2)*(T-t))/(sigma*np.sqrt(T-t))
    V=A*np.exp(-r*(T-t))*F(d2)
    return V


from matplotlib.pyplot import figure
figure(figsize=(10, 8), dpi=80)

T=1
m=1000
d=T/m #step
sigma=0.2
K=0.95 #strike price
A=1.2 #payoff in case the option is in the money
S0=1
r=0.1
My_W=Wiener_Process(T,m)
My_GBM=Geometric_Wiener_Process(My_W,r,sigma,S0)
My_GBM.plot()

V=np.zeros([m,1])
t=np.linspace(0,T,m)
for i in range(0,m):
    V[i]=Value_of_cash_or_nothing_call_option(i*d,My_GBM.S[i],K,m,r,sigma,T,A)
plt.plot(t,V,'m')

Delta=np.zeros([m,1])
for i in range(0,m):
    Delta[i]=Delta_of_cash_or_nothing_call_option(i*d,My_GBM.S[i],K,m,r,sigma,T,A)
plt.plot(t,Delta)
plt.show()
