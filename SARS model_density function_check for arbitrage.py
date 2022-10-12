import random
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mpl_toolkits import mplot3d

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

def a_hat(K,alpha,S0,beta,r,T,t0):
    SF=S0*np.exp(r*(T-t0))
    A=(SF*K)**((1-beta)/2)
    B=1+((1-beta)**2)*((np.log(SF/K))**2)/24+((1-beta)**4)*((np.log(SF/K))**4)/1920
    res=alpha/(A*B)
    return res

def c_hat(K,gamma,alpha,S0,r,T,t0,beta):
    SF=S0*np.exp(r*(T-t0))
    A=(SF*K)**((1-beta)/2)
    B=np.log(SF/K)*gamma/alpha
    res=A*B
    return res

def g_hat(K,gamma,alpha,S0,r,T,t0,beta,rho):
    x=c_hat(K,gamma,alpha,S0,r,T,t0,beta)
    ar=(np.sqrt(1-2*rho*x +x**2) + x - rho)/(1-rho)
    res=np.log(ar)
    return res
    
    
def imp_volatility(K,gamma,alpha,S0,r,T,t0,beta,rho):
    a=a_hat(K,alpha,S0,beta,r,T,t0)
    c=c_hat(K,gamma,alpha,S0,r,T,t0,beta)
    g=g_hat(K,gamma,alpha,S0,r,T,t0,beta,rho)
    SF=S0*np.exp(r*(T-t0))
    A=1+T*((((1-beta)**2)*alpha**2)/(24*(SF*K)**(1-beta))+rho*beta*gamma*alpha/(4*(SF*K)**((1-beta)/2))+(2-3*rho**2)*gamma**2/24)
    res=a*c*A/g
    return res


def d_1(S0,r,K,T,t0,alpha,beta,gamma,rho):
    sigma=imp_volatility(K,gamma,alpha,S0,r,T,t0,beta,rho)
    res=(np.log(S0/K)+(r+0.5*sigma**2)*(T-t0))/(sigma*np.sqrt(T-t0))
    return res

def d_2(S0,r,K,T,t0,alpha,beta,gamma,rho):
    sigma=imp_volatility(K,gamma,alpha,S0,r,T,t0,beta,rho)
    res=(np.log(S0/K)+(r+0.5*sigma**2)*(T-t0))/(sigma*np.sqrt(T-t0))-sigma*np.sqrt(T-t0)
    return res

def F(d):
    return norm.cdf(d)

def V_c(S0,r,K,T,t0,alpha,beta,gamma,rho):
    d1=d_1(S0,r,K,T,t0,alpha,beta,gamma,rho)
    d2=d_2(S0,r,K,T,t0,alpha,beta,gamma,rho)
    V_c=S0*F(d1)-np.exp(-r*(T-t0))*K*F(d2)
    return V_c

def check_call_spread_arbitrage(V_c_function):
    check=0
    for f in V_c_function.range:
        if f>0:
            V_c_function.plot()
            return 1
    return 0

def check_call_calendar_arbitrage(V_c_function):
    check=0
    for f in V_c_function.range:
        if f<0:
            V_c_function.plot()
            return 1
    return 0

def check_call_butterfly_arbitrage(V_c_function):
    check=0
    for f in V_c_function.range:
        if f<0:
            V_c_function.plot()
            return 1
    return 0    

K=np.linspace(0.001,10,1000)
Y=[]
T=2.5
t0=0
beta=0.5
S0=5.6
r=0.015
alpha=0.2
rho=-0.7
gamma=0.35
for k in K:
    temp=V_c(S0,r,k,T,t0,alpha,beta,gamma,rho)
    Y.append(temp)
V_c_function=function(K,Y)
First_der=V_c_function.derivative()
density=First_der.derivative()
for i in range(0,len(K)):
    density.range[i]=density.range[i]*np.exp(r*(T-t0))

check=check_call_spread_arbitrage(First_der)
if(check==0):
    print("No spread arbitrage")
else:
    print("There is spread arbitrage")
check=check_call_butterfly_arbitrage(density)
if(check==0):
    print("No butterfly arbitrage")
else:
    print("There is butterfly arbitrage")
    
T=np.linspace(0.001,5,1000)
for k in K:
    V_c_k=[]
    for t in T:
        temp=V_c(S0,r,k,t,t0,alpha,beta,gamma,rho)
        V_c_k.append(temp)
    V_c_function=function(T,V_c_k)
    First_der=V_c_function.derivative()
    check=check_call_calendar_arbitrage(First_der)
    if(check==1):
        print("There is calendar arbitrage")
        break
if check==0:
    print("No calendar arbitrage")
