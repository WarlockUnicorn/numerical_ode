# =============================================================================
#   ***** marshak_waves.py *****
#   Python script to solve and plot Marshak Waves (n=0 & n=3).
#
#   Author:     Ryan Clement
#   Created:    September 2021
#
#   Change Log:
#   Who:
#   Date:       MM/DD/YYY
#   What:
#
#   Who:
#   Date:       MM/DD/YYYY
#   What:
# =============================================================================

# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
# from scipy.integrate import solve_bvp
from scipy.optimize import fsolve

# ***************************
#    n=0 Marshak Wave       *
# ***************************
tMaxN0 = 1.2311729703 # LANL
# MS0 = lambda t,s: np.dot(np.array([[0,1], [0,-(t + 12.0*(s[0]**2)*s[1])/(4.0*s[0]**3)]]),s)
# G0 = lambda tm,t: np.power((3.0*tm*(tm-t)/4.0), (1.0/3.0))

def ms0(t,s):
    m11 = -(t/(4.0*(s[0]**3)) + 3.0*s[1]/s[0])
    # m = np.array([[0,1], [0,m11]])
    # return np.dot(m,s)
    return np.array([s[1],s[1]*m11])

def gEnd0(tm,t):
    return np.power(3.0*tm*(tm-t)/4.0, 1.0/3.0)

def gpEnd0(tm,t):
    tmp = np.power(3.0*tm*(tm-t)/4.0, -2.0/3.0)
    return -tm*tmp/4.0

def shootM0(tm):
    t0 = tm - dt0
    gE = gEnd0(tm,t0)
    vE = gpEnd0(tm,t0)
    y0 = np.array([gE,vE]).ravel()
    sol = solve_ivp(ms0,(t0,0),y0,t_eval=t_eval,method='Radau',rtol=1e-9) 
    g0 = sol.y[0]
    return g0[-1] - 1.0

def fin0(tm):
    tm, = fsolve(shootM0,tm)
    print(f'n=0 tm = {tm:.8f}')
    t0 = tm-dt0
    t_eval = np.linspace(t0,0,1000)
    gE = gEnd0(tm,t0)
    vE = gpEnd0(tm,t0)
    y0 = np.array([gE,vE]).ravel()
    sol = solve_ivp(ms0,(t0,0),y0,t_eval=t_eval,method='Radau',rtol=1e-9)
    tSol = sol.t
    tSol = np.insert(tSol,0,tm)
    gSol = sol.y[0]
    gSol = np.insert(gSol,0,0)
    return tSol, gSol


# ***************************
#    n=3 Marshak Wave       *
# ***************************
tMaxN3 = 1.1199349391 # LANL

def ms3(t,s):
    m11 = -(t/(7.0*(s[0]**6)) + 6.0*s[1]/s[0])
    m = np.array([[0,1], [0,m11]])
    return np.dot(m,s)

def gEnd3(tm,t):
    return np.power(6.0*tm*(tm-t)/7.0, 1.0/6.0)

def gpEnd3(tm,t):
    tmp = np.power(6.0*tm*(tm-t)/7.0, -5.0/6.0)
    return -tm*tmp/7.0

def shootM3(tm):
    t0 = tm - dt0
    gE = gEnd3(tm,t0)
    vE = gpEnd3(tm,t0)
    y0 = np.array([gE,vE]).ravel()
    sol = solve_ivp(ms3,(t0,0),y0,t_eval=t_eval,method='Radau',rtol=1e-9) 
    g0 = sol.y[0]
    return g0[-1] - 1.0

def fin3(tm):
    tm, = fsolve(shootM3,tm)
    print(f'n=3 tm = {tm:.8f}')
    t0 = tm-dt0
    t_eval = np.linspace(t0,0,1000)
    gE = gEnd3(tm,t0)
    vE = gpEnd3(tm,t0)
    y0 = np.array([gE,vE]).ravel()
    sol = solve_ivp(ms3,(t0,0),y0,t_eval=t_eval,method='Radau',rtol=1e-9)
    tSol = sol.t
    tSol = np.insert(tSol,0,tm)
    gSol = sol.y[0]
    gSol = np.insert(gSol,0,0)
    return tSol, gSol


# Higher Order Forward Stuff (n=0)
# def forward0():
#     dt = 1e-6
#     g0 = 1.0
#     gp0 = -1.0/4.0
#     t_eval = np.linspace(0,1.3,1000)
#     sol = solve_ivp(ms0,(0,1.3),[g0,gp0],t_eval=t_eval,method='Radau',rtol=1e-9)
#     tSol = sol.t
#     gSol = sol.y[0]
#     return tSol, gSol

# Second Order Finite Difference Stuff (n=0)
# def ryan0(i,dx,gi,gim1):
#     ci = 2.0*(gim1**4) - 4.0*(gi**4) - i*(dx**2)*gim1
#     F = lambda g: i*dx*dx*g + 2.0*(g**4) + ci
#     g, = fsolve(F,gi)
#     return g
# gim1 = 1.0
# gim12 = 1.0
# dx = 1e-5
# gi = gim1 - dx/4.0
# gi2 = gim12 - dx/2.0
# g = 1.0
# g2 = 1.0
# i = 1
# tSR0 = [0,dx]
# gSR0 = [gim1,gi]
# gSR02 = [gim12,gi2]
# while (g > 0):
#     g = ryan0(i,dx,gi,gim1)
#     g2 = ryan0(i,dx,gi2,gim12)
#     i+=1
#     gim1 = gi
#     gi = g
#     tSR0.append(i*dx)
#     gSR0.append(g)
#     gim12 = gi2
#     gi2 = g2
#     gSR02.append(g2)


# BVP Method
# Nice try but no...
# initF = lambda t, tm: (1.0 - (t/tm)**(6))**(1./6.)
# tm = 1.23
# dt = 1e-12
# tEnd = tm - dt
# t_bvp = np.linspace(0, tEnd, 999999)
# g_guess0 = initF(t_bvp,tm)
# g_guess = np.vstack((g_guess0,np.zeros(999999)))
# # print('g_guess: ',g_guess)
# def bc(bcL,bcR):
#     return np.array([bcL[0]-1.0,bcR[0]])
# g_sol = solve_bvp(ms0,bc,t_bvp,g_guess,max_nodes=1000000,tol=1e-6)
# g_y = g_sol.y[0]
# g_x = g_sol.x
# # print(g_y)
# plt.plot(t_bvp,g_guess0,color='blue')
# plt.plot(g_x,g_y,color='red')
# plt.show()




# ***************************
#    Plot Results           *
# ***************************
plotIt = 1
if plotIt:
    t_eval = np.linspace(1.0,0,1000)
    dt0 = 1e-6
    tS0, gS0 = fin0(1.2)

    t_eval = np.linspace(1.0,0,1000)
    dt0 = 1e-6
    tS3, gS3 = fin3(1.1)
    # tF, gF = forward0()

    plt.figure()
    plt.xlim(0,1.5)
    plt.ylim(0,1.1)
    plt.plot(tS0,gS0,label='n=0 B',color='blue')
    # plt.plot(tF,gF,label='n=0 F',color='red')
    # plt.plot(tSR0,gSR0,label='n=0 F',color='black')
    # plt.plot(tSR0,gSR02,label='n=0 F2',color='red')
    plt.plot(tS3,gS3,label='n=3 B',color='red')
    plt.xlabel(r'Normalized Distance: $\mathbf{\xi}$')
    plt.ylabel(r'Normalized Temperature: $\mathbf{g}$')
    plt.title('Marshak Waves')
    plt.legend()
    plt.show()




















# Test problem with known solution (v0=34.5 m/s)
# M = lambda s: np.array( [ [0,1], [0,-9.8/s[1]] ] )
# F = lambda t,s: np.dot(M(s),s)
# F = lambda t,s: np.dot(np.array( [ [0,1], [0,-9.8/s[1]] ] ),s)
# y0 = 0
# v0 = 34.5 # First guess
# t_eval = np.linspace(0,5,10)
# def shoot(v0):
#     sol = solve_ivp(F,[0,5],[y0,v0],t_eval=t_eval)
#     y = sol.y[0]
#     return y[-1] - 50

# # v0, = fsolve(shoot,v0)
# print(f'v0 = {v0}')
# sol = solve_ivp(F,[0,5],[y0,v0],t_eval=t_eval)
# yp = sol.y[1]
# print(f'v5 = {yp[-1]}')

# plt.figure()
# plt.plot(sol.t,sol.y[0])
# plt.plot(0,0,'ro')
# plt.plot(5,50,'ro')
# plt.xlabel('Time (s)')
# plt.ylabel('Altitude (m)')
# plt.title(f'Solution: v={v0:.7f} m/s')
# plt.show()

# F = lambda t,s: np.dot(np.array( [ [0,1], [0,-9.8/s[1]] ] ),s)
# y5 = 50
# v5 = -14.5 # First guess
# t_eval = np.linspace(5,0,10)
# sol = solve_ivp(F,[5,0],[y5,v5],t_eval=t_eval)
# plt.figure()
# plt.plot(sol.t,sol.y[0])
# plt.plot(0,0,'ro')
# plt.plot(5,50,'ro')
# plt.xlabel('Time (s)')
# plt.ylabel('Altitude (m)')
# plt.title(f'Solution: v={v5:.7f} m/s')
# plt.show()
# print(sol.t)
# print(sol.y[0])



