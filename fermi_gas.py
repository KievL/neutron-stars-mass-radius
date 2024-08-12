import numpy as np
from scipy.constants import pi, hbar, c

def k_to_x(k, m):
  return k/(m*c)

def e0(m):
  return (m**4)*(c**5)/((pi**2)*(hbar**3))

def pressure_density_x(x, m):
  return e0(m)*(((2*(x**3) - 3*x)*((1+ x**2)**(0.5)) + 3*np.arcsinh(x))/24)

def energy_density_x(x, m):
  return e0(m)*(((2*(x**3) + x)*(1+ x**2)**(1/2)) - (np.arcsinh(x)))/8