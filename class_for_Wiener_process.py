import random
import math
import matplotlib.pyplot as plt

class Wiener_Process:
    def __init__(self,T,m):
        W=[0]
        S=0
        delta_t=T/m
        for i in range(0,m):
            inc=random.gauss(0,delta_t)
            S=S+inc
            W.append(S)
        self.W=W
    
    def print(self):
        l=len(self.W)
        for i in range(0,l):
            print(self.W[i])
        
    def plot(self):
        fig = plt.figure()
        plt.plot(self.W)
        plt.grid()


T=2
m=1000
My_Wiener=Wiener_Process(T,m)
My_Wiener.plot()
