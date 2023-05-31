import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd

#  Define the objective function with two FM terms. Add more FM terms as
# required to get a good fit. (H +- Hc) => try all combination of +- and
# chose the combination with lowest resnorm and without any discontinuity
def objective_func(x, H):
    FM1 = (2 * x[2] / np.pi) * np.arctan((H + x[0]) / x[0] * np.tan(np.pi * x[4] / 2))
    FM2 = (2 * x[3] / np.pi) * np.arctan((H + x[1]) / x[1] * np.tan(np.pi * x[5] / 2))
    AFM = x[6] * H
    return FM1 + FM2 + AFM

# Provide the initial guess for the parameters. Change as per base material
# to get a good convergence
initial_guess = [7.1, 6.1, 0.3, 0.4, 0.01, 0.02, 2.3e-5]

# Load the experimental data using pandas
data = pd.read_csv('experimental_data.csv')
x = data['H'].values
y = data['M'].values

# Define the sum of squares error function
def error_func(params):
    return np.sum((objective_func(params, x) - y) ** 2)

# Perform curve fitting
result = minimize(error_func, initial_guess, method='Nelder-Mead')

# Extract the optimized parameter values
fit_params = result.x
resnorm = result.fun

Hc1_optimized, Hc2_optimized, Ms1_optimized, Ms2_optimized, S1_optimized, S2_optimized, chi_optimized = fit_params

# Display the optimized parameter values
print('Optimized Parameter Values:')
print('Hc1:', Hc1_optimized)
print('Hc2:', Hc2_optimized)
print('Ms1:', Ms1_optimized)
print('Ms2:', Ms2_optimized)
print('S1:', S1_optimized)
print('S2:', S2_optimized)
print('chi:', chi_optimized)
print('resnorm:', resnorm)

# Separation of the desired fitted model variables 
FM1 = [(2 * Ms1_optimized / np.pi) * np.arctan((val + Hc1_optimized) / Hc1_optimized * np.tan(np.pi * S1_optimized / 2)) for val in x]
FM2 = [(2 * Ms2_optimized / np.pi) * np.arctan((val + Hc2_optimized) / Hc2_optimized * np.tan(np.pi * S2_optimized / 2)) for val in x]
AFM = chi_optimized * x
Total_M = [FM1[i] + FM2[i] + AFM[i] for i in range(len(x))]

# Plotting (optional)
plt.plot(x, Total_M)
plt.scatter(x, y, s=10, edgecolors='black', facecolors='blue', linewidth=0.5)
plt.show()
