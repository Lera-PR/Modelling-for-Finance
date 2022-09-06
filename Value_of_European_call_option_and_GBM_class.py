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
        

        

T=1 #these are parameters for European call option
m=1000
d=T/m # time step
sigma=0.4
K=10 #strike price
S0=10
r=0.05

My_W=Wiener_Process(T,m) #create a Wiener process 

My_GBM=Geometric_Wiener_Process(My_W,r,sigma,S0) #and create GBM based on the process and parameters
My_GBM.plot()


def F(d):
    return norm.cdf(d) #this is a cumulitive function of a standard normal distribution

def Delta_of_European_call_option(t,S,K,m,r,sigma,T):            # calculate Delta of European_call_option given current time
    d1=(np.log(S/K)+(T-t)*(r+0.5*sigma**2))/(sigma*np.sqrt(T-t)) #t, current price of stock S, and specified parameters
    return F(d1)                                                 #of the option

def Gamma_of_European_call_option(t,S,K,m,r,sigma,T):                 # calculate Gamma of European_call_option given current time
    d1=(np.log(S/K)+(T-t)*(r+0.5*sigma**2))/(sigma*np.sqrt(T-t))      #t, current price of stock S, and specified parameters
    Gamma=np.exp(-d1**2/2)/(S*sigma*np.sqrt(T-t)*np.sqrt(2*math.pi))  #of the option
    return Gamma

def Value_of_European_call_option(t,S,K,m,r,sigma,T):             # calculate value of European_call_option given current time
    d1=(np.log(S/K)+(T-t)*(r+0.5*sigma**2))/(sigma*np.sqrt(T-t))  #t, current price of stock S, and specified parameters
    d2=d1-sigma*np.sqrt(T-t)                                      #of the option
    V=S*F(d1)-K*np.exp(-r*(T-t))*F(d2)
    return V


fig = plt.figure()                                                #here we calculate and plot values of European call 
ax = plt.axes(projection='3d')                                    #option given GBM for stock prices
t=np.linspace(d,T-d,m)
V=np.zeros([m,1])
for i in range(0,m):
    V[i]=Value_of_European_call_option(t[i],My_GBM.S[i],K,m,r,sigma,T)
ax.scatter(t, My_GBM.S, V,color='red')                                

S_max=22                                                         #here we calculate and plot values of European call
t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
S=np.linspace(0.01,S_max,m)                                      #price
X,Y = np.meshgrid(t, S)
V=np.zeros([m,m])
for i in range(0,m):
    for j in range(0,m):
        V[j][i]=Value_of_European_call_option(t[i],S[j],K,m,r,sigma,T)

ax.plot_surface(X, Y, V)
ax.set_xlabel('t')
ax.set_ylabel('S(t)')
ax.set_zlabel('V(t,S(t))');

ax.view_init(20, -120)



fig = plt.figure()                                              #here we calculate and plot Deltas of European call
ax = plt.axes(projection='3d')                                  #option given GBM for stock prices
t=np.linspace(d,T-d,m)
D=np.zeros([m,1])
for i in range(0,m):
    D[i]=Delta_of_European_call_option(t[i],My_GBM.S[i],K,m,r,sigma,T)
ax.scatter(t, My_GBM.S, D,color='red')

S_max=22                                                        #here we calculate and plot Deltas of European call
t=np.linspace(d,T-d,m)                                          #option for all pairs (t,S) where t<T and S is a stock
S=np.linspace(0.01,S_max,m)                                     #price
X,Y = np.meshgrid(t, S)
D=np.zeros([m,m])
for i in range(0,m):
    for j in range(0,m):
        D[j][i]=Delta_of_European_call_option(t[i],S[j],K,m,r,sigma,T)

ax.plot_surface(X, Y, D)
ax.set_xlabel('t')
ax.set_ylabel('S(t)')
ax.set_zlabel('Delta');

ax.view_init(20, -120)


fig = plt.figure()                                               #here we calculate and plot Gammas of European call
ax = plt.axes(projection='3d')                                   #option given GBM for stock prices
t=np.linspace(d,T-d,m)
Gamma=np.zeros([m,1])
for i in range(0,m):
    Gamma[i]=Gamma_of_European_call_option(t[i],My_GBM.S[i],K,m,r,sigma,T)
ax.scatter(t, My_GBM.S, Gamma,color='red')

S_max=22                                                         #here we calculate and plot Gammas of European call
t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
S=np.linspace(0.01,S_max,m)                                      #price
X,Y = np.meshgrid(t, S)
Gamma=np.zeros([m,m])
for i in range(0,m):
    for j in range(0,m):
        Gamma[j][i]=Gamma_of_European_call_option(t[i],S[j],K,m,r,sigma,T)

ax.plot_surface(X, Y, Gamma)
ax.set_xlabel('t')
ax.set_ylabel('S(t)')
ax.set_zlabel('Gamma');

ax.view_init(20, -120)





