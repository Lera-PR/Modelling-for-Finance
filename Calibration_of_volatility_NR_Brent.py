import random
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mpl_toolkits import mplot3d


def F(d):
    return norm.cdf(d) #this is a cumulitive function of a standard normal distribution

def Value_of_European_call_option(t,S,K,r,sigma,T): # value of a European call option
    d1=(np.log(S/K)+(T-t)*(r+0.5*sigma**2))/(sigma*np.sqrt(T-t))
    d2=d1-sigma*np.sqrt(T-t)
    V=S*F(d1)-K*np.exp(-r*(T-t))*F(d2)
    return V



#I calibrate parameter for volatility in the next two function using Combined Root Finding algorithm
#and using Brent's method. I am looking for the root of the function g(), and need its derivative for N-R method. 
#Value function V() for European call option is derived with BS model

def g(t0,S0,K,T,r,V_mkt,sigma): #this is function I optimise
    g=V_mkt - Value_of_European_call_option(t0,S0,K,r,sigma,T)
    return g

def dg(t0,S0,K,T,r,V_mkt,sigma): #this is the derivative of the function I optimise
    d1=(np.log(S0/K)+(T-t0)*(r+0.5*sigma**2))/(sigma*np.sqrt(T-t0))
    d2=d1+sigma*np.sqrt(T-t0)
    dg=-np.sqrt(T-t0)*K*np.exp(-r*(T-t0))*np.exp(-d2**2/2)
    return dg

    
def Newton_Raphson_upd(t0,S0,K,T,r,V_mkt):
    sigma_l=0
    sigma_r=1 #this would be huge volatility
    
    g_l=g(t0,S0,K,T,r,V_mkt,sigma_l)
    g_r=g(t0,S0,K,T,r,V_mkt,sigma_r) #calculate the values of g() at the ends of the interval
    
    if g_l*g_r>0: #if the signs are the same, there is no root on this interval (use monotonicity of g())
        print("No root at this interval")
    else:
        sigma_k=(sigma_l+sigma_r)/2 #first 'guess' for the calibrated sigma
        
        g_k=g(t0,S0,K,T,r,V_mkt,sigma_k) 
        dg_k=dg(t0,S0,K,T,r,V_mkt,sigma_k)
        delta=-g_k/dg_k #calculate everything for N-R update
        
        while np.abs(delta/sigma_k)>0.00001: #while we haven't detected the root we update the sigma
            new_sigma=sigma_k+delta
            
            if new_sigma <sigma_l or new_sigma>sigma_r: #if new sigma is outside of the initial interval, we update the interval
                g_new=g(t0,S0,K,T,r,V_mkt,new_sigma)
                
                if g_new*g_l>0:
                    sigma_l=sigma_k
                if g_new*g_l<0:
                    sigma_r=sigma_k
    
                new_sigma=(sigma_l+sigma_r)/2 #and calculate new sigma
        
            g_new=g(t0,S0,K,T,r,V_mkt,new_sigma)
            dg_new=dg(t0,S0,K,T,r,V_mkt,new_sigma)
            delta=-g_new/dg_new
            sigma_k=new_sigma #make N-R update
    return sigma_k

def Brent_method(t0,S0,K,T,r,V_mkt):
    sigma_l=0
    sigma_r=1 #this would be huge volatility
    g_l=g(t0,S0,K,T,r,V_mkt,sigma_l)
    g_r=g(t0,S0,K,T,r,V_mkt,sigma_r) #calculate the values of g() at the ends of the interval
    
    if g_l*g_r>0: #if the signs are the same, there is no root on this interval (use monotonicity of g())
        print("No root at this interval")
    else:
        sigma_k=(sigma_l+sigma_r)/2 #first 'guess' for the calibrated sigma
        sigma_k_1=(sigma_l+sigma_r)/4
        sigma_k_2=3*(sigma_l+sigma_r)/4
        
        g_k=g(t0,S0,K,T,r,V_mkt,sigma_k)
        g_k_1=g(t0,S0,K,T,r,V_mkt,sigma_k_1)
        g_k_2=g(t0,S0,K,T,r,V_mkt,sigma_k_2)
        
        while np.abs(g_k)>0.00001: #while we haven't detected the root we update the sigma
            new_sigma=sigma_k*g_k_1*g_k_2/((g_k-g_k_1)*(g_k-g_k_2))
            new_sigma=new_sigma+sigma_k_1*g_k*g_k_2/((g_k_1-g_k_2)*(g_k_1-g_k))
            new_sigma=new_sigma+sigma_k_2*g_k*g_k_1/((g_k_2-g_k)*(g_k_2-g_k_1))
            
            if new_sigma <sigma_l or new_sigma>sigma_r: #if new sigma is outside of the initial interval, we update the interval
                g_new=g(t0,S0,K,T,r,V_mkt,new_sigma)
                
                if g_new*g_l>0:
                    sigma_l=sigma_k
                if g_new*g_l<0:
                    sigma_r=sigma_k
    
                new_sigma=(sigma_l+sigma_r)/2 #and calculate new sigma
            
            sigma_k_2=sigma_k_1
            sigma_k_1=sigma_k
            sigma_k=new_sigma
            g_k=g(t0,S0,K,T,r,V_mkt,new_sigma)
            g_k_1=g(t0,S0,K,T,r,V_mkt,sigma_k_1)
            g_k_2=g(t0,S0,K,T,r,V_mkt,sigma_k_2)
    return sigma_k
        
T=1
K=120 #strike price
S0=100
r=0.05
V_mkt=2
t0=0

sigma=Newton_Raphson_upd(t0,S0,K,T,r,V_mkt)
print(sigma)

sigma=Brent_method(t0,S0,K,T,r,V_mkt)
print(sigma)
