#Here I will use homemade MLE to predict future stock prices for TESLA
import numpy as np
import random
import math
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('TSLA.csv')
closing_prices = df.Close #here we have closing prices of TESLA stock
plt.plot(closing_prices)
plt.grid()


#apply log-transformation to get to a simple Brownian motion, because we assume that closing prices follow GBM
l_closing_prices=np.zeros([len(closing_prices),1])
for i in range(len(closing_prices)):
    l_closing_prices[i]=math.log(closing_prices[i])


m=len(closing_prices)
delta_t=1
mu=(l_closing_prices[m-1] - l_closing_prices[0])/(m*delta_t) #MLE of historical drift mu

sqsigma=0
for i in range(0,m-1):
    sqsigma=sqsigma+np.power(l_closing_prices[i+1]-l_closing_prices[i]-mu*delta_t,2)

sqsigma=sqsigma/(m*delta_t) #MLE of historical volatility sigma^2
sigma=math.sqrt(sqsigma)
print(mu,sigma)

#now we will predict future prices
n=70 #number of periods we predict prices for
N=7 #number of Monte Carlo runs
X=np.zeros([n,1])
plt.plot(closing_prices[1700:]) #I plot only part of the closing prices to make the picture nicer
for i in range(0,N): #N Monte Carlo runs
    X[0]=random.gauss(l_closing_prices[-1]+mu*delta_t,sqsigma*delta_t)
    for j in range(1,n):
        X[j]=random.gauss(X[j-1]+mu*delta_t,sqsigma*delta_t) #we know conditional distribution of X_t+1 on X_t
    X=np.exp(X) #X will store predicted prices
    plt.plot(range(m,m+n),X,'r')

plt.grid()
plt.show()
