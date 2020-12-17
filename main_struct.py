import math as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
import logging
logging.basicConfig(level=logging.DEBUG, format=" %(asctime)s - %(levelname)s - %(message)s")

# %%
#from BEM_Functions import*
logging.debug('start of program')
headers = ['r', 'pitch_s', 'mass_d', 'EIy', 'EIz', 'twist']
df_bladestruc = pd.read_csv('bladestruc_mod.txt', sep=',',names=headers, skiprows=1)   

p = df_bladestruc.pitch_s.to_numpy()
beta = df_bladestruc.twist.to_numpy()
EIy = df_bladestruc.EIy.to_numpy()
EIz = df_bladestruc.EIz.to_numpy()

V = 11
#data = np.genfromtxt('loads' + str(V) + '.txt', delimiter=',')
data=pd.read_csv('loads11_mod.txt', sep=',', skiprows=(1))

print('stop')

x = data.iloc[:, 0]
pz = data.iloc[:, 1] * 10**3  # pN
py = data.iloc[:, 2] * 10**3  # pT

fig0=plt.figure()
plt.plot(x,pz,label='Pz [kN/m]')
plt.xlabel('radial position')
plt.ylabel('[kN/m]')
plt.grid()


plt.plot(x,py,label='Py [kN/m]')
plt.xlabel('radial position')
plt.ylabel('[kN/m]')
plt.grid()
plt.legend()
plt.show()


Ty = np.zeros(len(x))
Tz = np.zeros(len(x))
My = np.zeros(len(x))
Mz = np.zeros(len(x))



for j in range(len(x) - 1):
    i = len(x) - 1 - j
    Ty[i - 1] = Ty[i] + 0.5 * (py[i - 1] + py[i]) * (x[i] - x[i - 1])
    Tz[i - 1] = Tz[i] + 0.5 * (pz[i - 1] + pz[i]) * (x[i] - x[i - 1])

    My[i - 1] = My[i] - Tz[i] * (x[i] - x[i - 1]) - (1 / 6 * pz[i - 1] + 1 / 3 * pz[i]) * (x[i] - x[i - 1])**2
    Mz[i - 1] = Mz[i] + Ty[i] * (x[i] - x[i - 1]) + (1 / 6 * py[i - 1] + 1 / 3 * py[i]) * (x[i] - x[i - 1])**2
# %%
for hide_plot in range(1):
    fig1=plt.figure()
    plt.title('Ty & Tz')
    plt.plot(x, Ty / 10**6, label='Ty [kN]', color='red', linestyle='--')
    plt.plot(x, Tz / 10**6, label='Tz [kN]',color='orange')
    plt.xlabel('radial position')
    plt.ylabel('[kN]')
    plt.grid()
    plt.legend()
    
    fig2=plt.figure()
    plt.title('My & Mz')
    plt.plot(x, My / 10**6, label='My [kN/M]',color='green',linestyle='--')
    plt.plot(x, Mz / 10**6, label='Mz [kN/M]', color='lime')
    plt.legend()
    plt.xlabel('radial position')
    plt.ylabel('[kN/M]')
    plt.grid()
    plt.show()

# %%

M1 = np.zeros(len(x))
M2 = np.zeros(len(x))
K1 = np.zeros(len(x))
K2 = np.zeros(len(x))
Ky = np.zeros(len(x))
Kz = np.zeros(len(x))
EI1 = EIy
EI2 = EIz
for i in range(len(x)):
    M1[i] = My[i] * m.cos(np.deg2rad(beta[i] + p[i])) - Mz[i] * m.sin(np.deg2rad(beta[i] + p[i]))
    M2[i] = My[i] * m.sin(np.deg2rad(beta[i] + p[i])) + Mz[i] * m.cos(np.deg2rad(beta[i] + p[i]))

    K1[i] = M1[i] / EI1[i]
    K2[i] = M2[i] / EI2[i]

    Kz[i] = -K1[i] * m.sin(np.deg2rad(beta[i] + p[i])) + K2[i] * m.cos(np.deg2rad(beta[i] + p[i]))
    Ky[i] = K1[i] * m.cos(np.deg2rad(beta[i] + p[i])) + K2[i] * m.sin(np.deg2rad(beta[i] + p[i]))


thetay = np.zeros(len(x))
thetaz = np.zeros(len(x))
uy = np.zeros(len(x))
uz = np.zeros(len(x))
for i in range(len(x) - 1):
    thetay[i + 1] = thetay[i] + 0.5 * (Ky[i + 1] + Ky[i]) * (x[i + 1] - x[i])
    thetaz[i + 1] = thetaz[i] + 0.5 * (Kz[i + 1] + Kz[i]) * (x[i + 1] - x[i])

    uy[i + 1] = uy[i] + thetaz[i] * (x[i + 1] - x[i]) + (1 / 6 * Kz[i + 1] + 1 / 3 * Kz[i]) * (x[i + 1] - x[i])**2
    uz[i + 1] = uz[i] - thetay[i] * (x[i + 1] - x[i]) - (1 / 6 * Ky[i + 1] + 1 / 3 * Ky[i]) * (x[i + 1] - x[i])**2

plt.plot(x, uy, label='Uy')
plt.plot(x, uz, label='Uz')
plt.xlabel('r [m]')
plt.ylabel('Deflection [m]')
plt.yticks(np.arange(0, 10, 1))
plt.grid()
plt.legend()
plt.show()
