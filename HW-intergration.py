import matplotlib.pyplot as plt
import numpy as np

def f(t):
    return np.exp(-t**2)

# limits of integration and step size/increments
N=30
a=0.0
b= 3.0
dx=0.1

# Calculation of integral using trapezoid method
k= np.arange(1, N)
s= 0.5*f(a)+ 0.5*f(b)+ np.sum(f(a+k*dx))
print(s*dx)

# Array of x-values from 0.0 to 3.0
x_v= np.arange(a, 3.1, 0.1)
f_v=f(x_v)

# The cumulative integral of E(x)
cut_area=0.5*dx*(f_v[:-1]+f_v[1:])  # individual area cuts
E_x= np.concatenate(([0.0], np.cumsum(cut_area)))


# graph of E(x) with highlighted area under curve and approximate integral value
plt.figure(figsize=(8, 5))
plt.plot(x_v, E_x, color='purple', linestyle="solid")
plt.fill_between(x_v, E_x, color='lightblue', alpha=0.5, label='Integrated Area')
plt.title(r'Numerical Evalution of E(x)=$\int_0^x (e^{-t^2})dt$ via Trapezoid Method')
plt.scatter(x_v[-1], E_x[-1], color='red', s=50, zorder=5)
plt.annotate(f'Approximate plateau: {float(E_x[-1]):.6f}',
             xy=(x_v[-1], E_x[-1]),
             xytext=(x_v[-1] - 0.9, E_x[-1] - 0.1),
             arrowprops=dict(arrowstyle="->", color='black'))
plt.xlabel('x', fontsize=13)
plt.ylabel('E(x)', fontsize=13)
plt.show()
