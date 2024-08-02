import numpy as np

def phi(x):
        return 0.1818 * np.exp(-3.2*x) + 0.5099*np.exp(-0.9423*x) + \
        0.2802 * np.exp(-0.4029*x) + 0.02817 * np.exp(-0.2016*x)

#z1, z2: atomic charges
# epsilon: softening radius
def make_ZBL_potential(z1, z2, epsilon = 0):
    a = 0.8854 * 0.529 / (z1**0.23+z2**0.23) # 0.529 = Bohr radius

    prefactor = 14.3983 # 1/(4pi*e0) * e^2 in units of eV * Angstrom 
    def f(r):
        return prefactor * z1 * z2 / np.sqrt(r*r+epsilon**2) * phi(r / a)
    return f
