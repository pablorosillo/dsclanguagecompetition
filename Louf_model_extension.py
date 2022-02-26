import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib as mpl

#--------------------------------------------------------------------------------------------
#----COEXISTENCE OF LANGUAGES THORUGH PARAMETER MODIFICATION. LOUF ET AL. MODEL--------------
#--------------------------------------------------------------------------------------------

#Master equations of the Louf et al. model in terms of the r parameter
def func(var,t,r,q,s):
    x, y=var
    dydt=[ r*(1 - x - y)*s*(x + q*(1 - x - y)) - x*(1 - s)*(y + (1 - q)*(1 - x - y))  , r*(1 - x - y)*(1 - s)*(y + (1 - q)*(1 - x - y)) - y*s*(x + q*(1 - x - y))]
    return dydt


mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rc('font',**{'family':'serif','serif':['Times']})



t=np.linspace(2020, 2200, 1000) #TIme range in which the extension of the evolution is going to be done
y0=[0.100411738340564, 0.753073207228281] #Initial condition (2020 data from data fitting) of valenciá
#Parameters modified to reach coexistence
q=0.5
s=0.5
mu=1/50
c=1/20
r=mu/(c*(1-mu))


sol=odeint(func, y0, t, args=(r, q, s)) #Solve the master equation


#Code for plotting the solution of the master equation
plt.plot(t, sol[:,0], color=(191/255, 151/255, 198/255), linewidth=5, label=r'$p_{A}$')
plt.plot(t, sol[:,1], color=(195/255, 219/255, 164/255), linewidth=5, label=r'$p_{B}$')
plt.plot(t, 1-sol[:,0]-sol[:,1], color=(166/255, 223/255, 247/255), linewidth=5, label=r'$p_{AB}$')
plt.xlabel('Year', fontsize=15)
plt.legend(prop={'size': 13})
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.show()

