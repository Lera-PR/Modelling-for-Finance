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
        for i in range(1,len(self.range)-1):
            temp=(self.range[i+1]-self.range[i])/(self.domain[i+1]-self.domain[i])
            der.append(temp)
        #print(self.range[1:len(self.range)],der)
        derivative=function(self.domain[1:len(self.range)],der)
        return derivative
    
    def value_at_point(self,x):
        i=0
        for d in self.domain:
            temp=d-x
            if abs(temp)<0.001:
                res=self.range[i]
                return res
            i=i+1
        print("No such argument in the domain")
        return -111111

def root_NR(My_function,l,r):
    #will try finding all the roots with MC
    v_l=My_function.value_at_point(l)
    v_r=My_function.value_at_point(r)
    m=-11111
    if (v_l*v_r<0): #there at least one root at this random interval
        m=(l+r)/2 # at least one of the roots is present at this interval
        v_m=My_function.value_at_point(m)
            
        der=My_function.derivative()
        dv_m=der.value_at_point(m) # calculate value of the derivative at the point m
            
        delta=-v_m/dv_m #calculate everything for N-R update
            
        while np.abs(delta)>0.001: #while we haven't detected the root we update the candidate
            root_to_be=m+delta
                
            if root_to_be<l or root_to_be>r: #if the candidate for the root is outside of the initial interval, we update the interval
                v_new=My_function.value_at_point(root_to_be)
                if v_new*v_l>0:
                    l=m
                if v_new*v_l<0:
                    r=m
                root_to_be=(l+r)/2 #and calculate new sigma
            
            v_new=My_function.value_at_point(root_to_be)
            dv_new=der.value_at_point(root_to_be)
            delta=-v_new/dv_new
            m=root_to_be #make N-R update
    return m

def root_BM(My_function,l,r):
    v_l=My_function.value_at_point(l)
    v_r=My_function.value_at_point(r)
    k=-11111
    if v_l*v_r<0: #if the signs are the same, there is no root on this interval (use monotonicity of g())
        k=(l+r)/2 #first 'guess' for the calibrated sigma
        k_1=(l+r)/4
        k_2=3*(l+r)/4
        
        v_k=My_function.value_at_point(k)
        v_k_1=My_function.value_at_point(k_1)
        v_k_2=My_function.value_at_point(k_2)
        
        while np.abs(v_k)>0.001: #while we haven't detected the root we update the sigma
            root_to_be=k*v_k_1*v_k_2/((v_k-v_k_1)*(v_k-v_k_2))
            root_to_be=root_to_be+k_1*v_k*v_k_2/((v_k_1-v_k_2)*(v_k_1-v_k))
            root_to_be=root_to_be+k_2*v_k*v_k_1/((v_k_2-v_k)*(v_k_2-v_k_1))
            
            if root_to_be<l or root_to_be>r: #if new sigma is outside of the initial interval, we update the interval
                if v_new*v_l>0:
                    l=k
                if v_new*v_l<0:
                    r=k
    
                root_to_be=(l+r)/2 #and calculate new sigma
            
            k_2=k_1
            k_1=k
            k=root_to_be
            v_k=My_function.value_at_point(k)
            v_k_1=My_function.value_at_point(k_1)
            v_k_2=My_function.value_at_point(k_2)
    return k


#Newton-Raphson method returns one root of the function on a given interval, if it exists.
#Brent method does the same. They both return -11111 if there is no root at the interval (or if there is an even number of points)

l=-5
r=5
x=np.linspace(l,r,10000)
y=np.zeros([len(x),1])

for i in range(len(x)):
    y[i]=(np.exp(x[i])+np.exp(-x[i]))/2 - 2*x[i]

My_func=function(x,y)
My_func.plot()

m = np.random.choice(My_func.domain, size=1)
v_m=My_func.value_at_point(m)
while(v_m>0):
    m = np.random.choice(My_func.domain, size=1)
    v_m=My_func.value_at_point(m)
    
root=root_NR(My_func,l,m)
print(root)
root=root_NR(My_func,m,r)
print(root)

root=root_BM(My_func,l,m)
print(root)
root=root_BM(My_func,m,r)
print(root)

#here I find two roots of the given function. However, I have to use additional information - known number of roots.
