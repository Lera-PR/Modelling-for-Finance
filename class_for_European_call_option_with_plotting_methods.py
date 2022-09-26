import random
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mpl_toolkits import mplot3d

def F(d):
    return norm.cdf(d) #this is a cumulitive function of a standard normal distribution

class European_call_option:
    def __init__(self,t,S,K,m,r,sigma,T): #just want all the parameters for some European call option
        self.t0=t
        self.S0=S
        self.K=K
        self.m=m
        self.r=r
        self.sigma=sigma
        self.T=T
        
    def Value_of_European_call_option(self):
        d1=(np.log(self.S0/self.K)+(self.T-self.t0)*(self.r+0.5*self.sigma**2))/(self.sigma*np.sqrt(self.T-self.t0))
        d2=d1-self.sigma*np.sqrt(self.T-self.t0)
        V=self.S0*F(d1)-self.K*np.exp(-self.r*(self.T-self.t0))*F(d2)
        return V
    
    def Delta_of_European_call_option(self):
        d1=(np.log(self.S0/self.K)+(self.T-self.t0)*(self.r+0.5*self.sigma**2))/(self.sigma*np.sqrt(self.T-self.t0)) #t, current price of stock S, and specified parameters
        return F(d1)
    
    def Gamma_of_European_call_option(self):                 # calculate Gamma of European_call_option given current time
        d1=(np.log(self.S0/self.K)+(self.T-self.t0)*(self.r+0.5*self.sigma**2))/(self.sigma*np.sqrt(self.T-self.t0))      #t, current price of stock S, and specified parameters
        Gamma=np.exp(-d1**2/2)/(self.S0*self.sigma*np.sqrt(self.T-self.t0)*np.sqrt(2*math.pi))  #of the option
        return Gamma
    
    def dV_dK_of_European_call_option(self):
        d1=(np.log(self.S0/self.K)+(self.T-self.t0)*(self.r+0.5*self.sigma**2))/(self.sigma*np.sqrt(self.T-self.t0))
        d2=d1-self.sigma*np.sqrt(self.T-self.t0)
        
        dV_dK=-self.S0*np.exp(-d1**2/2)/(np.sqrt(2*math.pi)*self.K*self.sigma*np.sqrt(self.T-self.t0))
        dV_dK=dV_dK-F(d2)*np.exp(-self.r*(self.T-self.t0))
        dV_dK=dV_dK+np.exp(-self.r*(self.T-self.t0))*np.exp(-d2**2/2)/(np.sqrt(2*math.pi)*self.sigma*np.sqrt(self.T-self.t0))
        return dV_dK
    
    def d2V_dK2_of_European_call_option(self):
        d1=(np.log(self.S0/self.K)+(self.T-self.t0)*(self.r+0.5*self.sigma**2))/(self.sigma*np.sqrt(self.T-self.t0))
        d2=d1-self.sigma*np.sqrt(self.T-self.t0)
        
        res=-self.S0*np.exp(-d1**2/2)*(d1/(self.sigma*np.sqrt(self.T-self.t0))-1)/(np.sqrt(2*math.pi)*np.sqrt(self.T-self.t0)*self.sigma*K**2)
        res=res+np.exp(-self.r*(self.T-self.t0))*np.exp(-d2**2/2)/(np.sqrt(2*math.pi)*np.sqrt(self.T-self.t0)*K*self.sigma)
        res=res+np.exp(-self.r*(self.T-self.t0))*np.exp(-d2**2/2)/(np.sqrt(2*math.pi)*(self.T-self.t0)*self.sigma**2*K)
        return res
    
    
class Dynamics_of_European_call_option:
    def __init__(self,t,S,K,m,r,sigma,T):
        self.t0=t
        self.S0=S
        self.K=K
        self.m=m
        self.r=r
        self.sigma=sigma
        self.T=T
        
    def plot_dynamics_of_value_of_ECO_t_S(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')                                    #option given GBM for stock prices                             

        S_max=22                                                         #here we calculate and plot values of European call
        t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
        S=np.linspace(0.01,S_max,m)                                      #price
        X,Y = np.meshgrid(t, S)
        V=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],S[j],self.K,self.m,self.r,self.sigma,self.T)
                V[j][i]=EC.Value_of_European_call_option()

        ax.plot_surface(X, Y, V)
        ax.set_xlabel('t')
        ax.set_ylabel('S(t)')
        ax.set_zlabel('V(t,S(t))');

        ax.view_init(20, -120)
        
    def plot_dynamics_of_value_of_ECO_t_K(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')                                    #option given GBM for stock prices                             

        K_max=15                                                         #here we calculate and plot values of European call
        t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
        K=np.linspace(0.01,K_max,m)                                      #price
        X,Y = np.meshgrid(t, K)
        V=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],self.S0,K[j],self.m,self.r,self.sigma,self.T)
                V[j][i]=EC.Value_of_European_call_option()

        ax.plot_surface(X, Y, V)
        ax.set_xlabel('t')
        ax.set_ylabel('K')
        ax.set_zlabel('V(t,S(t))');

        ax.view_init(20, -120)
        
    def plot_dynamics_of_delta_of_ECO_t_S(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')
        
        S_max=22                                                        #here we calculate and plot Deltas of European call
        t=np.linspace(d,T-d,m)                                          #option for all pairs (t,S) where t<T and S is a stock
        S=np.linspace(0.01,S_max,m)                                     #price
        X,Y = np.meshgrid(t, S)
        D=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],S[j],self.K,self.m,self.r,self.sigma,self.T)
                D[j][i]=EC.Delta_of_European_call_option()

        ax.plot_surface(X, Y, D)
        ax.set_xlabel('t')
        ax.set_ylabel('S(t)')
        ax.set_zlabel('Delta');

        ax.view_init(20, -120)
        
    def plot_dynamics_of_delta_of_ECO_t_K(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')
        
        K_max=15                                                        #here we calculate and plot Deltas of European call
        t=np.linspace(d,T-d,m)                                          #option for all pairs (t,S) where t<T and S is a stock
        K=np.linspace(0.01,K_max,m)                                     #price
        X,Y = np.meshgrid(t, K)
        D=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],self.S0,K[j],self.m,self.r,self.sigma,self.T)
                D[j][i]=EC.Delta_of_European_call_option()

        ax.plot_surface(X, Y, D)
        ax.set_xlabel('t')
        ax.set_ylabel('K')
        ax.set_zlabel('Delta');

        ax.view_init(20, -120)
        
    def plot_dynamics_of_gamma_of_ECO_t_S(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')
        
        S_max=22                                                         #here we calculate and plot Gammas of European call
        t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
        S=np.linspace(0.01,S_max,m)                                      #price
        X,Y = np.meshgrid(t, S)
        Gamma=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],S[j],self.K,self.m,self.r,self.sigma,self.T)
                Gamma[j][i]=EC.Gamma_of_European_call_option()

        ax.plot_surface(X, Y, Gamma)
        ax.set_xlabel('t')
        ax.set_ylabel('S(t)')
        ax.set_zlabel('Gamma');

        ax.view_init(20, -120)
        
        
    def plot_dynamics_of_gamma_of_ECO_t_K(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')
        
        K_max=15                                                         #here we calculate and plot Gammas of European call
        t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
        K=np.linspace(0.01,K_max,m)                                      #price
        X,Y = np.meshgrid(t, K)
        Gamma=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],self.S0,K[j],self.m,self.r,self.sigma,self.T)
                Gamma[j][i]=EC.Gamma_of_European_call_option()

        ax.plot_surface(X, Y, Gamma)
        ax.set_xlabel('t')
        ax.set_ylabel('K')
        ax.set_zlabel('Gamma');

        ax.view_init(20, -120)
        
        
    def plot_dynamics_of_dV_dK_of_ECO_t_S(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')
        
        S_max=22                                                         #here we calculate and plot Gammas of European call
        t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
        S=np.linspace(0.01,S_max,m)                                      #price
        X,Y = np.meshgrid(t, S)
        dV_dK=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],S[j],self.K,self.m,self.r,self.sigma,self.T)
                dV_dK[j][i]=EC.dV_dK_of_European_call_option()

        ax.plot_surface(X, Y, dV_dK)
        ax.set_xlabel('t')
        ax.set_ylabel('S(t)')
        ax.set_zlabel('dV_dK');

        ax.view_init(20, -120)
        
    def plot_dynamics_of_dV_dK_of_ECO_t_K(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')
        
        K_max=15                                                         #here we calculate and plot Gammas of European call
        t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
        K=np.linspace(0.01,K_max,m)                                      #price
        X,Y = np.meshgrid(t, K)
        dV_dK=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],self.S0,K[j],self.m,self.r,self.sigma,self.T)
                dV_dK[j][i]=EC.dV_dK_of_European_call_option()

        ax.plot_surface(X, Y, dV_dK)
        ax.set_xlabel('t')
        ax.set_ylabel('K')
        ax.set_zlabel('dV_dK');

        ax.view_init(20, -120)
        
    def plot_dynamics_of_d2V_dK2_of_ECO_t_S(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')
        
        S_max=22                                                         #here we calculate and plot Gammas of European call
        t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
        S=np.linspace(0.01,S_max,m)                                      #price
        X,Y = np.meshgrid(t, S)
        d2V_dK2=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],S[j],self.K,self.m,self.r,self.sigma,self.T)
                d2V_dK2[j][i]=EC.d2V_dK2_of_European_call_option()

        ax.plot_surface(X, Y, d2V_dK2)
        ax.set_xlabel('t')
        ax.set_ylabel('S(t)')
        ax.set_zlabel('d2V_dK2');

        ax.view_init(20, -120)
        
    def plot_dynamics_of_d2V_dK2_of_ECO_t_K(self):
        fig = plt.figure()                                                #here we calculate and plot values of European call 
        ax = plt.axes(projection='3d')
        
        K_max=15                                                         #here we calculate and plot Gammas of European call
        t=np.linspace(d,T-d,m)                                           #option for all pairs (t,S) where t<T and S is a stock
        K=np.linspace(0.01,K_max,m)                                      #price
        X,Y = np.meshgrid(t, K)
        d2V_dK2=np.zeros([m,m])
        for i in range(0,m):
            for j in range(0,m):
                EC=European_call_option(t[i],self.S0,K[j],self.m,self.r,self.sigma,self.T)
                d2V_dK2[j][i]=EC.d2V_dK2_of_European_call_option()

        ax.plot_surface(X, Y, d2V_dK2)
        ax.set_xlabel('t')
        ax.set_ylabel('K')
        ax.set_zlabel('d2V_dK2');

        ax.view_init(20, -120)


T=1 #these are parameters for European call option
m=100
d=T/m # time step
sigma=0.4
K=10 #strike price
S0=10
r=0.05


My_European_Call_option_Dynamics=Dynamics_of_European_call_option(0,S0,K,m,r,sigma,T)
My_European_Call_option_Dynamics.plot_dynamics_of_value_of_ECO_t_S()
My_European_Call_option_Dynamics.plot_dynamics_of_value_of_ECO_t_K()

My_European_Call_option_Dynamics.plot_dynamics_of_dV_dK_of_ECO_t_S()
My_European_Call_option_Dynamics.plot_dynamics_of_dV_dK_of_ECO_t_K()

My_European_Call_option_Dynamics.plot_dynamics_of_d2V_dK2_of_ECO_t_S()
My_European_Call_option_Dynamics.plot_dynamics_of_d2V_dK2_of_ECO_t_K()
