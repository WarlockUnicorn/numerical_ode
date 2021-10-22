# =============================================================================
#  ***** earth_orbit_sun_ch3.py *****
#  Python script to solve and plot numerical solution for Newtonian
#  orbit of Earth around Sun using RK4.
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


# Earth Orbit Parameters
aEarth = 1.0                                                # semi-major in AU
eEarth = 0.0167086                                          # eccentricity
bEarth = aEarth*np.sqrt(1 - eEarth**2)
mu     = 4.0*np.pi**2                                       # mu = 4Pi**2 (units!)
twoPi  = 2.0*np.pi                      
vEarth = twoPi*np.sqrt((1 + eEarth)/(1.0 - eEarth))     # VisViva @ Perihelion

# Earth Orbit as Ellipse
def rOrbit(a,e,t):
    return a*(1.0 - e**2)/(1.0 + e*np.cos(t))
def xOrbit(a,e,t):
    return rOrbit(a,e,t)*np.cos(t)
def yOrbit(a,e,t):
    return rOrbit(a,e,t)*np.sin(t)
t = np.linspace(0,twoPi,100)
xE = xOrbit(aEarth,eEarth,t)
yE = yOrbit(aEarth,eEarth,t)

# Earth Orbit Initial Conditions
x0  = 1.0 - eEarth
y0  = 0.0
vx0 = 0.0
vy0 = vEarth
nSteps = 12    # to complete one orbit
dt  = 1.0/nSteps
dt6 = dt/6.0

# Newton's Equations
def calcR3(x,y):
    return (x**2 + y**2)**(3/2)
# xD = vx
# yD = vy
# vxD = -mu*x/r3
# vyD = -mu*y/r3

# Earth Orbit via RK4
def update(x,y,vx,vy):
    # k1
    r3 = calcR3(x, y)
    k1x = vx
    k1y = vy
    k1vx = -mu*x/r3
    k1vy = -mu*y/r3
    xu = x + k1x*dt/2.0
    yu = y + k1y*dt/2.0
    # k2
    r3 = calcR3(xu, yu)
    k2x = vx + k1vx*dt/2.0
    k2y = vy + k1vy*dt/2.0
    k2vx = -mu*xu/r3
    k2vy = -mu*yu/r3
    # k3
    xu = x + k2x*dt/2.0
    yu = y + k2y*dt/2.0
    r3 = calcR3(xu, yu)
    k3x = vx + k2vx*dt/2.0
    k3y = vy + k2vy*dt/2.0
    k3vx = -mu*xu/r3
    k3vy = -mu*yu/r3
    # k4
    xu = x + k3x*dt
    yu = y + k3y*dt
    r3 = calcR3(xu, yu)
    k4x = vx + k3vx*dt
    k4y = vy + k3vy*dt
    k4vx = -mu*xu/r3
    k4vy = -mu*yu/r3
    # Weighted Sum
    kx = k1x + 2.0*k2x + 2.0*k3x + k4x
    ky = k1y + 2.0*k2y + 2.0*k3y + k4y
    kvx = k1vx + 2.0*k2vx + 2.0*k3vx + k4vx
    kvy = k1vy + 2.0*k2vy + 2.0*k3vy + k4vy
    # Updates
    xN = x + dt6*kx
    yN = y + dt6*ky
    vxN = vx + dt6*kvx
    vyN = vy + dt6*kvy
    return xN,yN,vxN,vyN

xN = x0
yN = y0
vxN = vx0
vyN = vy0
xArr = [x0]
yArr = [y0]
vxArr = [vx0]
vyArr = [vy0]
tArr = [0]
for i in range(nSteps):
    xN, yN, vxN, vyN = update(xN,yN,vxN,vyN)
    xArr.append(xN)
    yArr.append(yN)
    vxArr.append(vxN)
    vyArr.append(vyN)
    tArr.append(dt*(i+1))

print(f'xAnal: {x0}, xComp: {xArr[-1]}')
print('X Percent Difference: ', 100*(x0-xArr[-1])/x0)
print(f'yAnal: {y0}, yComp: {yArr[-1]}')
# print('Y Percent Difference: ', 100*(y0-yArr[-1])/y0)
print(f'vxAnal: {vx0}')
print(f'vyAnal: {vy0}')

    
# Orbit plot
figE, axE = plt.subplots()
axE.set_title('Earth Orbit')
axE.set_xlabel('x (AU)')
axE.set_ylabel('y (AU)')
axE.plot(xArr,yArr,'+',color='red',label='RK4')
axE.plot(xE, yE, color='blue',label='Analytic')
axE.plot(aEarth*eEarth,0,'o',color='orange',label='Sun')
axE.grid('true')
axE.axis('equal')
axE.set(xlim=(-1.1,1.1),ylim=(-1.1,1.1))
axE.legend()
# figE.savefig('earth_orbit.png')

# Angular momentum plot
xA = np.array(xArr)
yA = np.array(yArr)
vxA = np.array(vxArr)
vyA = np.array(vyArr)
rxvAnal = (x0*vy0 - y0*vx0)
rxvAnalArr = [rxvAnal]*len(tArr)
rxvAA = np.array(rxvAnalArr)
rxvCompArr = xA*vyA - yA*vxA
perDiffArr = 100*(rxvAA - rxvCompArr)/rxvAA
figL, axL = plt.subplots()
axL.set_title('Angular Momentum')
axL.set_xlabel('t (periods)')
axL.set_ylabel(r'r $x$ v')
axL.plot(tArr,rxvAA,color='blue',label='Exact')
axL.plot(tArr,rxvCompArr,color='red',label='Computed')
axL.legend()
# figL.savefig('earth_orbit_Lcon.png')

# Transformed Polar Example
eEarth = 0.0167086                                          # eccentricity
mu     = 4.0*np.pi**2                                       # mu = 4Pi**2 (units!)
aS     = mu*(1.0 - eEarth**2)
ma     = 1.0/(1.0 - eEarth**2)
twoPi  = 2.0*np.pi
print(f'Analytical Minimum Radius: {1.0-eEarth} au')
print(f'Analytical Maximum Radius: {1.0+eEarth} au')                    
def q2r(q):
    return 1.0/q
# IC
q0  = 1.0/(1.0 - eEarth)
z0  = 0.0
nSteps = 12
dTh = twoPi/nSteps
def newQZ_FE(dTh, q, z):
    c = 1.0/(1.0 - eEarth**2)
    qN  = q + z*dTh
    zN  = z + (c-q)*dTh
    return qN, zN
def newQZ_RK(dTh, q, z):
    c = 1.0/(1.0 - eEarth**2)
    dTh6 = dTh/6.0
    aq0 = q
    az0 = z
    kq1 = az0
    kz1 = c - aq0
    aq1 = aq0 + dTh*kq1/2.0
    az1 = az0 + dTh*kz1/2.0
    kq2 = az1
    kz2 = c - aq1
    aq2 = aq0 + dTh*kq2/2.0
    az2 = az0 + dTh*kz2/2.0
    kq3 = az2
    kz3 = c - aq2
    aq3 = aq0 + dTh*kq3
    az3 = az0 + dTh*kz3
    kq4 = az3
    kz4 = c - aq3
    kq  = kq1 + 2.0*kq2 + 2.0*kq3 + kq4
    kz  = kz1 + 2.0*kz2 + 2.0*kz3 + kz4
    qN  = aq0 + dTh6*kq
    zN  = az0 + dTh6*kz
    return qN, zN
rArr = [q2r(q0)]
thArr = [0.0]
qO = q0
zO = z0
for i in range(1,nSteps+1):
    qN, zN = newQZ_RK(dTh,qO,zO)
    rArr.append(q2r(qN))
    thArr.append(i*dTh)
    qO = qN
    zO = zN
figP, axP = plt.subplots(subplot_kw={'projection': 'polar'})
axP.plot(thArr,rArr,'.',color='blue',label='Earth')
axP.plot(0,0,'o',color='orange',label='Sun')
axP.set_title('Earth Orbit')
axP.legend(loc='best', bbox_to_anchor=(1.1, 1.2))
# figP.savefig('earth_orbit_PolarTrans.png')

    
print(f'Min computed radius: {min(rArr)} au')
print(f'Max computed radius: {max(rArr)} au')

