from scipy.constants import pi, G, c
import numpy as np
import math

def dp_dr_tov(**kwargs):  
  p_r_ = kwargs['p_r_'] # p(r)
  
  r_ = kwargs['r_'] # radius
  
  m_r_ = kwargs['m_r_'] # m_bar(r) : m(r) / sun_mass
  
  spline_ = kwargs['spline'] # cubic spline : energy_density(p(r))

  e_r_ = spline_(p_r_)
  t1 = (1 + (p_r_/e_r_))
  t2 = 1+ (4*pi*(r_**3)*p_r_)/(m_r_*(c**2))
  t3 = 1 - (2*G*m_r_)/(r_*(c**2))

  return -G*spline_(p_r_)*m_r_*t1*t2/(((r_*c)**2)*t3)

def dm_dr_tov(**kwargs):  
  e_r_ = kwargs['e_r_'] # energy_density(r)
  
  r_ = kwargs['r_'] # radius

  return (4*pi*(r_**2)*e_r_)/(c**2)

# Runge-Kutta 4
def edo_runge_kutta(**kwargs):  
  gdm = kwargs['gdm'] # dm_dr
  gdp = kwargs['gdp'] # dp_dr
  m0_ = kwargs['m0_'] # initial mass
  p0_ = kwargs['p0_'] # initial pressure density
  e0_ = kwargs['e0_'] # initial energy density
  r_ini = kwargs['r_ini'] # initial radius, close to 0
  r_fim = kwargs['r_fim'] # final radius, bigger than the expected max radius
  num_steps = kwargs['num_steps'] # steps
  spline_ = kwargs['spline_'] # spline energy_density(p)

  h = (r_fim - r_ini)/num_steps # step-size
  
  r_values = np.linspace(r_ini, r_fim, num_steps) # radius array
  m_values = np.zeros(num_steps)# mass array
  p_values = np.zeros(num_steps) # p array
  e_values = np.zeros(num_steps) # energy density array

  # setting initial values
  m_values[0] = m0_
  p_values[0] = p0_
  e_values[0] = e0_

  # computing values
  for i in range(1, num_steps):
    k1m = gdm(e_r_=e_values[i-1], r_=r_values[i-1], p_r_ = p_values[i-1], spline=spline_)
    k1p = gdp(e_r_=e_values[i-1],r_=r_values[i-1], p_r_ = p_values[i-1], m_r_=m_values[i-1], spline=spline_)

    k2m = gdm(e_r_=spline_(p_values[i-1] + k1p/2), r_=r_values[i-1] + h/2, p_r_ = p_values[i-1] + k1p/2, spline=spline_)
    k2p = gdp(e_r_=spline_(p_values[i-1] + k1p/2),r_=r_values[i-1]+ h/2, p_r_ = p_values[i-1] + k1p/2, m_r_=m_values[i-1] + k1m/2, spline=spline_)

    k3m = gdm(e_r_=spline_(p_values[i-1] + k2p/2), r_=r_values[i-1] + h/2, p_r_ = p_values[i-1] + k2p/2, spline=spline_)
    k3p = gdp(e_r_=spline_(p_values[i-1] + k2p/2),r_=r_values[i-1]+ h/2, p_r_ = p_values[i-1] + k2p/2, m_r_=m_values[i-1] + k2m/2, spline=spline_)

    k4m = gdm(e_r_=spline_(p_values[i-1] + k3p), r_=r_values[i-1] + h, p_r_ = p_values[i-1] + k3p, spline=spline_)
    k4p = gdp(e_r_=spline_(p_values[i-1] + k3p),r_=r_values[i-1]+ h, p_r_ = p_values[i-1] + k3p, m_r_=m_values[i-1] + k3m, spline=spline_)


    m_values[i] = m_values[i-1] + h * (k1m + (k2m*2) + (k3m*2) + k4m)/6
    p_values[i] = p_values[i-1] + h * (k1p + (k2p*2) + (k3p*2) + k4p)/6
    e_values[i] = spline_(p_values[i])

    # end condition (if pressure density is less or equal to 0)
    if(p_values[i]<=0):
      r_values = r_values[:i+1]
      p_values = p_values[:i+1]
      e_values = e_values[:i+1]
      m_values = m_values[:i+1]

      break

  return m_values, r_values, p_values

# mass-radius relation
def mass_radius(**kwargs):
  gdm_ = kwargs['gdm'] # dm_dr
  gdp_ = kwargs['gdp'] # dp_dr
  m0_1 = kwargs['m0_'] # initial mass
  spline_1 = kwargs['spline'] # cubic spline energy_density(p)
  r_fim_ = kwargs['r_fim_'] # final radius
  ps_ini = kwargs['ps_ini'] # initial pressure
  ps_fim = kwargs['ps_fim'] # final pressure
  modelo = kwargs['modelo'] # model name

  num_steps_=5000 # steps

  results_model = [] # store the results

  ps_log = np.logspace(np.log10(ps_ini), np.log10(ps_fim), 100) # target pressure density values

  
  # progress prompt
  len_ps = len(ps_log)
  progress_printed = False
  print("Computing model: ", modelo)
  processos_rk = []
  for pos,i in enumerate(ps_log):

    # printing progress
    progress = math.floor((pos+1)*100/len_ps)
    if(progress%10==0):
      if(not progress_printed):
        print(progress,"% finished")
        progress_printed=True
    else:
      progress_printed=False

      ep_0 = spline_1(i)

    # resolving with rk4
    mas, r_c , p_trash = edo_runge_kutta(gdm = gdm_, gdp=gdp_ ,m0_ = m0_1, p0_= i, e0_= spline_1(i), r_ini = 1e-5, r_fim=r_fim_,num_steps =  num_steps_, spline_ = spline_1)

    # appending results in results array
    results_model.append([mas[-1], r_c[-1], i])

    results_polytrope_array = np.array(results_model)

  print("Finished for ", modelo)

  # separating each list
  masses = results_polytrope_array[:,0].tolist()
  radii = results_polytrope_array[:,1].tolist()
  pressures = results_polytrope_array[:,2].tolist()

  # cleaning invalid values
  index = 0
  len_masses = len(masses)

  while(index<len_masses):
    mas = masses[index]
    if(math.isnan(mas) or math.isinf(mas)):
      masses.pop(index)
      pressures.pop(index)
      radii.pop(index)
      index-=1
      len_masses-=1
    index+=1

  return pressures, masses, radii