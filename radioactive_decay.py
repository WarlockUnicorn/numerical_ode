# =============================================================================
#  ***** radioactive_decay.py *****
#  Python script to solve first order ODE for radioactive decay.
#
#  A = -dN/dt = lambda*N
#
#  Author:     Ryan Clement
#  Created:    July 2021
#
#  Change Log:
#  Who:
#  Date:       MM/DD/YYY
#  What:
#
#  Who:
#  Date:       MM/DD/YYYY
#  What:
# =============================================================================

import numpy as np
from matplotlib import pyplot as plt


### Functions
def anal(n0,l,t):
    return n0*np.exp(-l*t)

## Integrators    
def euler(n0,l,scale):
    nOld = n0 
    nNew = n0
    dt = 1.0/(l*scale)
    tArr = [0.0]
    nArr = [n0]
    t = 0.0
    while (nNew >= nStop):
        nNew = nOld*(1.0 - l*dt)
        nOld = nNew
        nArr.append(nOld)
        t += dt
        tArr.append(t)
    print(f"Forward Euler stop time, dt: {t}, {dt}")
    return tArr, nArr

def eulerStepsB(n0,l,scale,nSteps):
    nOld = n0 
    nNew = n0
    dt = 1.0/(l*scale)
    tArr = [0.0]
    nArr = [n0]
    t = 0.0
    for i in range(nSteps):
        nNew = nOld/(1.0 + l*dt)
        nOld = nNew
        nArr.append(nOld)
        t += dt
        tArr.append(t)
    print(f"Backward Euler stop time, dt: {t}, {dt}")
    return tArr, nArr

def eulerSteps(n0,l,scale,nSteps):
    nOld = n0 
    nNew = n0
    dt = 1.0/(l*scale)
    tArr = [0.0]
    nArr = [n0]
    t = 0.0
    for i in range(nSteps):
        nNew = nOld*(1.0 - l*dt)
        nOld = nNew
        nArr.append(nOld)
        t += dt
        tArr.append(t)
    print(f"Forward Euler stop time, dt: {t}, {dt}")
    return tArr, nArr

# Backward Euler
def eulerB(n0,l,scale):
    nOld = n0 
    nNew = n0
    dt = 1.0/(l*scale)
    tArr = [0.0]
    nArr = [n0]
    t = 0.0
    while (nNew >= nStop):
        nNew = nOld/(1.0 + l*dt)
        nOld = nNew
        nArr.append(nOld)
        t += dt
        tArr.append(t)
    print(f"Backward Euler stop time, dt: {t}, {dt}")
    return tArr, nArr

def trap(n0,l,scale):
    nOld = n0
    nNew = n0
    dt = 1.0/(l*scale)
    d = 1.0/(2.0*scale)
    rat = (1.0 - d)/(1.0 + d)
    tArr = [0.0]
    nArr = [n0]
    t = 0.0
    while (nNew >= nStop):
        nNew = nOld*rat
        nOld = nNew
        t += dt
        tArr.append(t)
        nArr.append(nOld)
    print(f"Trapezoid Rule stop time, dt: {t}, {dt}")
    return tArr, nArr

def rk2(n0,l,scale):
    nOld = n0
    nNew = n0
    dt = 1.0/(l*scale)
    dt2 = dt**2
    l2 = l**2
    tArr = [0.0]    # t0 = 0.0
    nArr = [n0]
    t = 0.0
    while (nNew >= nStop):
        nNew = nOld*(1.0 - l*dt + l2*dt2/2.0)
        nOld = nNew
        t += dt
        tArr.append(t)
        nArr.append(nOld)
    print(f"RK2 stop time, dt: {t}, {dt}")
    return tArr, nArr

# Plotters
def plotFE(n0,l):
    tA = np.linspace(0,tAnal,100)
    nA = anal(n0,l,tA)
    t2, n2 = euler(n0,l,2.0)
    t4, n4 = euler(n0,l,4.0)
    t20, n20 = euler(n0,l,20.0)
    figE, axE = plt.subplots()
    axE.set_title('Radioactive Decay')
    axE.set_xlabel('t')
    axE.set_ylabel('Decay Material')
    axE.plot(t20, n20,'*',color='red',label='Euler 20')
    axE.plot(t4, n4,'*',color='blue',label='Euler 4')
    axE.plot(t2, n2,'x',color='purple',label='Euler 2')
    axE.plot(tA,nA,color='black',label='Analytic')
    axE.legend()
    # figE.savefig('radioactive_FE.png')

def plotUnstabFE(n0,l):
    tp3, np3 = eulerSteps(1.0/2.1,10)
    ts = 10.0/(l*(1.0/2.1))
    tA = np.linspace(0,ts,100)
    nA = anal(n0,l,tA)
    figEU, axEU = plt.subplots()
    axEU.set_title('Radioactive Decay')
    axEU.set_xlabel('t')
    axEU.set_ylabel('Decay Material')
    axEU.plot(tp3, np3,'-*',color='red',label='Euler 0.3')
    axEU.plot(tA,nA,color='black',label='Analytic')
    axEU.legend()
    # figEU.savefig('radioactive_unstable_FE.png')

def plotEulerB(n0,l):
    tA = np.linspace(0,tAnal,100)
    nA = anal(n0,l,tA)
    tb2, nb2 = eulerB(n0,l,2.0)
    tb4, nb4 = eulerB(n0,l,4.0)
    tb20, nb20 = eulerB(n0,l,20.0)
    figBE, axBE = plt.subplots()
    axBE.set_title('Radioactive Decay')
    axBE.set_xlabel('t')
    axBE.set_ylabel('Decay Material')
    axBE.plot(tb20, nb20,'*',color='red',label='Backward Euler 20')
    axBE.plot(tb4, nb4,'*',color='blue',label='Backward Euler 4')
    axBE.plot(tb2, nb2,'x',color='purple',label='Backward Euler 2')
    axBE.plot(tA,nA,color='black',label='Analytic')
    axBE.legend()
    # figBE.savefig('radioactive_BE.png')

# Unstable Case BE comparo
def plotBE_Comparo(n0,l):
    tp3, np3 = eulerSteps(n0,l,1.0/2.1,10)
    tp3B, np3B = eulerStepsB(n0,l,1.0/2.1,10)
    ts = 10.0/(l*(1.0/2.1))
    tA = np.linspace(0,ts,100)
    nA = anal(n0,l,tA)
    figBU, axBU = plt.subplots()
    axBU.set_title('Radioactive Decay')
    axBU.set_xlabel('t')
    axBU.set_ylabel('Decay Material')
    axBU.plot(tp3, np3,'-*',color='red',label='Forward Euler 0.3')
    axBU.plot(tp3B, np3B,'-*',color='blue',label='Backward Euler 0.3')
    axBU.plot(tA,nA,color='black',label='Analytic')
    axBU.legend()
    # figBU.savefig('radioactive_unstable_BE.png')

def plotFEvBE(n0,l):
    tA = np.linspace(0,tAnal,100)
    nA = anal(n0,l,tA)
    t2, n2 = euler(n0,l,2.0)
    tb2, nb2 = eulerB(n0,l,2.0)
    figC, axC = plt.subplots()
    axC.set_title('Radioactive Decay')
    axC.set_xlabel('t')
    axC.set_ylabel('Decay Material')
    axC.plot(t2, n2,'-*',color='red',label='Forward Euler 2')
    axC.plot(tb2, nb2,'-*',color='blue',label='Backward Euler 2')
    axC.plot(tA,nA,color='black',label='Analytic')
    axC.legend()
    # figC.savefig('radioactive_FE_BE.png')

def plotFEvBEvT(n0,l):
    tA = np.linspace(0,tAnal,100)
    nA = anal(l,n0,tA)
    t2, n2 = euler(n0,l,2.0)
    tb2, nb2 = eulerB(n0,l,2.0)
    tTrap, nTrap = trap(n0,l,2.0)
    figTrap, axTrap = plt.subplots()
    axTrap.set_title('Radioactive Decay')
    axTrap.set_xlabel('t')
    axTrap.set_ylabel('Decay Material')
    axTrap.plot(t2, n2,'-o',color='red',label='Forward Euler 2')
    axTrap.plot(tb2, nb2,'-o',color='blue',label='Backward Euler 2')
    axTrap.plot(tTrap, nTrap,'-o',color='purple',label='Trapezoid 2')
    axTrap.plot(tA,nA,color='black',label='Analytic')
    axTrap.legend()
    # figTrap.savefig('radioactive_FE_BE_Trap.png')


if __name__ == '__main__':
    n0 = 1.0
    l = 1.0                         # Decay constant (lambda)
    nStop = 0.01                    # Stop amount
    tAnal = -np.log(nStop)          # Analytic stop time
    print(f"Analytic Stop Time: {tAnal}")
    plotFEvBEvT(n0, l)
    plotBE_Comparo(n0, l)
     
