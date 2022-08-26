import random
import math
import matplotlib.pyplot as plt
import numpy as np

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
        
class Saving_Account: #simple deterministic process 
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
        

T=10
m=1000
mu=0.15 #parameters for GMB
sigma=0.1

r=0.05 #interest rate for my saving account
My_Saving_Account=Saving_Account(r,T,m)

Average_S_over_M=np.zeros([m,1]) #here I will keep my average discounted asset price path, calculated with Monte Carlo

n=100 #In this loop we generate n discounted GBM paths and plot 10 of them
for i in range(0,n):
    My_Wiener=Wiener_Process(T,m) #create Wiener process first
    My_Geometric_Wiener=Geometric_Wiener_Process(My_Wiener,mu,sigma) #use Wiener process to genearate GBM path
    My_S_over_M=np.divide(My_Geometric_Wiener.S,My_Saving_Account.M) #discount it with saving account
    Average_S_over_M=np.add(Average_S_over_M,My_S_over_M) 
    if(i<10):
        plt.plot(x, My_S_over_M,'b') #plot some dicounted process
        
for i in range(0,m):
    Average_S_over_M[i]=Average_S_over_M[i]/n #calculate average discounted process and plot it
    
plt.plot(x,Average_S_over_M,'r')
    
plt.grid()
plt.show()

#If drift mu>r, we have sub-martingale process. If mu=r, then martingale and if mu<r (why to buy this asset then?), it is a super-martingale
