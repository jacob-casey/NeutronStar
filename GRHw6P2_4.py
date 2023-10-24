#
from scipy.integrate import odeint as bob
import numpy as np
import matplotlib.pyplot as plt

# e = float(input('Enter e: '))
# l = float(input('Enter lz: '))
# r0 = float(input('Enter r0: '))
global e, l, r0, T, u0, t, r, P, u
e = 0.9622504
l = 3.6756
r0 = 12

n = 10001
T = np.linspace(0, 1000, n)

u0 = np.sqrt((e ** 2) - ((1 - 2 / r0) * (((l ** 2) / (r0 ** 2)) + 1)))

t = np.empty_like(T)
r = np.empty_like(T)
P = np.empty_like(T)
u = np.empty_like(T)

# initial conditions
z0 = [0, r0, 0, u0]
t[0] = 0
r[0] = r0
P[0] = 0
u[0] = u0

def pend(t, r, P, u, T):
    t= z0[0]
    r= z0[1]
    P= z0[2]
    u = z0[3]
    dtdT = e * ((1 - (2 / r)) ** (-1))*0
    drdT = u
    dPdT = l / (r ** 2)
    dudT = -(1 / r ** 2) + ((l ** 2) / (r ** 3)) - 3 * ((l ** 2) / (r ** 4))

    dzdT = [dtdT, drdT, dPdT, dudT]
    return dzdT


for i in range(1, n):

    z = bob(pend, z0, T, args=(e, l, r0))
    t[i] = z[0][0]
    r[i] = z[0][1]
    P[i] = z[0][2]
    u[i] = z[0][3]
    z0 = z[1]

fig1, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
ax1.plot(T, t, 'b-', label='t(Tau)')
ax1.set_ylabel('Time t')


ax2.plot(T, r, 'r--', label='R (Tau)')
ax2.set_ylabel('R')


ax3.plot(T, P, 'b--', label='P(Tau)')
ax3.set_ylabel('$\phi')


ax4.plot(T, u, 'g--', label='u(Tau)')
ax4.set_ylabel('4-Velocity Ur')
plt.xlabel('Proper Time Tau')

fig2, (ax5) = plt.subplots(1, subplot_kw=dict(projection = 'polar'),)
ax5.set_rlim([0,20])
ax5.plot(P,r)

plt.show()

