from scipy.constants import atomic_mass, e, c
from astropy.constants import M_sun

#constantes da seção no SI
n0 = 0.16/(1e-45) #nucleons /m³
E0 = 22*1e6*e #joule
sigma = 2.112 #admimensional
A = -118.2*1e6*e #joule
B = 65.39*1e6*e #joule
S0 = 30.0*1e6*e #joule
sun_mass = M_sun

def p_n(n):
  u = n/n0

  termo1 = 2*E0*(u**(5/3))/3
  termo2 = A*(u**2)/2
  termo3 = B*sigma*(u**(sigma+1))/(sigma+1)

  return n0*(termo1+termo2+termo3)

def p_n_assymetric(n, alpha):
  u = n/n0

  termo1 = p_n(n)
  termo2 = n0*alpha**2*( (((2**(2/3)) - 1)*E0*( (2*(u**(5/3))/3)-(u**2) ) ) + S0*(u**2))

  return termo1+termo2

def e_n_assymetric(n, alpha):
  u = n/n0

  termo1=atomic_mass*c**2                 
  termo2=(2**(2/3))*E0*(u**(2/3))
  termo3 = A*u/2
  termo4 = B*(u**sigma)/(1+sigma)
  termo5 = (S0 - (E0*((2**(2/3))-1)))*u

  return (termo1+termo2+termo3+termo4+termo5)*n