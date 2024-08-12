from fermi_gas import pressure_density_x, energy_density_x, k_to_x
from scipy.constants import m_n, m_p, m_e, c

def total_energy_density(kn,kp):
    xn = k_to_x(kn, m_n)
    xp = k_to_x(kp, m_p)
    xe = k_to_x(kp, m_e) # kp=ke pelo equilíbrio eletrostático

    e_proton=energy_density_x(xp, m_p)
    e_electron=energy_density_x(xe, m_e)
    e_neutron=energy_density_x(xn, m_n)

    return e_proton+e_electron+e_neutron

def total_pressure_density(kn, kp):
    xn = k_to_x(kn, m_n)
    xp = k_to_x(kp, m_p)
    xe = k_to_x(kp, m_e) # kp=ke pelo equilíbrio eletrostático

    p_proton=pressure_density_x(xp, m_p)
    p_electron=pressure_density_x(xe, m_e)
    p_neutron=pressure_density_x(xn, m_n)

    return p_proton+p_electron+p_neutron

def kp_kn(kn_):
    kn=kn_
    part1 = (((kn*c)**2)+((m_n**2)*(c**4))-((m_e**2)*(c**4)))**2
    part2 = 2*(m_p**2)*(c**4)*(((kn*c)**2)+((m_n**2)*(c**4))+((m_e**2)*(c**4)))
    part3 = ((m_p**4)*(c**8))
    part4 = 2*c*((((kn*c)**2)+((m_n**2)*(c**4)))**(1/2))

    return ((part1-part2+part3)**(1/2))/part4

def p_root(kn, kp, p):
    return total_pressure_density(kn, kp)-p

def root_kn(p2):
    a = 0
    b = 1.0e+3
    err = 1.0e-6

    fa = p_root(a,kp_kn(a),p2)
    fb = p_root(a,kp_kn(b),p2)

    if (fa*fb) > 0:
                print(" No convergence garanteed:" , p2)
    if (fa*fb) < 0:
        while (abs(b-a)/abs(b) > err):
            xm = (a + b)/2
            fm = p_root(xm,kp_kn(xm),p2)
            fa = p_root(a,kp_kn(a),p2)
            fb = p_root(b,kp_kn(b),p2)

            if fm == 0:
                break
            else:
                if (fa*fm) < 0:
                    b = xm
                else:
                    a = xm
        return xm