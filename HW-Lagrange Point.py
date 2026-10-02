import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const

G = const.G.value          # Newtons Gravitational  constant
M = const.M_earth.value   # Earth mass
m = 7.348e22         # Moon mass
w = 2.662e-6         # Angular velocity
R = 3.844e8          # distance from Earth to Moon


def f(r):
    return (G*M)/(r**2)-(G*m)/((R-r)**2)-(w**2)*r


def df(r):
    return -2*(G*M)/(r**3)-2*(G*m)/((R-r)**3)-(w**2)


def r_roots(r_guess, tolerance=1e4, N=100):
    r= r_guess
    for i in range(N):
        fr= f(r)
        dfr= df(r)
        r_next= r-(fr/dfr)
        if abs(r_next-r)< tolerance:
            return r_next
        r = r_next
    return r


r_i = R*0.8
L1_point = r_roots(r_i)

print(f"The distance r form Earth to the L1 point is: {L1_point:,.4e}")

r= np.linspace(R*0.05, R*0.95, 1000)
f_vals = f(r)
plt.figure(figsize=(10, 6))
plt.plot(r, f_vals, color='purple', linestyle='solid', label='Net Force $f(r)$')
plt.plot(L1_point, f(L1_point), 'bo', markersize=5, label=f'L1 Point ({L1_point:.4e} m)')
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.title('Net Gravitational Force and L1 Lagrange Point', fontsize=10)
plt.xlabel('Distance from Earth Center $r$ (meters)', fontsize=10)
plt.ylabel('Net Force Metric $f(r)$', fontsize=10)
plt.legend(fontsize=11)
plt.show()

