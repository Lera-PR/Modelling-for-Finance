import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
from mpl_toolkits import mplot3d



nflx_option = pd.read_csv("netflix_options.csv")
nflx_option["Expiration Date"] = pd.to_datetime(nflx_option["Expiration Date"])

dt.date(2019, 4, 8) #let's say it is a today date

nflx_option["timeToMaturity"] = (nflx_option["Expiration Date"].dt.date - dt.date(2019, 4, 8))/ dt.timedelta(days=1)
#here we added an additional column, calculating time to maturity

x = nflx_option["timeToMaturity"]
y = nflx_option["Strike"]
z = nflx_option["IV"]
fig= plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('time to maturity')
ax.set_ylabel('strike')
ax.set_zlabel('Implied vol')
ax.scatter(x, y, z, c='red', marker='o')
