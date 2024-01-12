import numpy as np
from numpy import sin, exp, pi
from sympy import var, plot_implicit
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------------------------------------------------

# # N = 2
# # dx = 1/(2**N-1)
# # k = np.arange(2**N)
# # y = []
# # for kval in k:
# #     sum = 0
# #     for j in range(2**N):
# #         sum += sin(dx * j * pi) * exp(dx*1j*j*k)
# #     y.append(sum)
    
# # print(y)

# # Define the initial condition function f(x)
# def f(x):
#     return (1/np.pi**2)*np.sin(np.pi * x)  # Just an example; you can change this function

# # Define the x range for plotting
# x = np.linspace(0, 1, 400)

# t = np.pi-0.1

# # Calculate u(x, t=1) using the implicit solution u(x, t) = f(x - u*t)
# def u_at_t_one(x, t=t):
#     return f(x - t * f(x))

# # Calculate u(x, t=1) values
# u_values = u_at_t_one(x)

# # Plotting
# plt.figure(figsize=(10, 5))
# plt.plot(x, u_values, label=r'$u(x, t=1)$', color='blue')
# plt.title(r'Plot of $u(x, t=1)$ for given initial condition $f(x)$')
# plt.xlabel(r'$x$')
# plt.ylabel(r'$u(x, t=1)$')
# plt.grid(True)
# plt.legend()
# plt.show()

# ------------------------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------------------------


