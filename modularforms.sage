# This code is ONLY for half-integer indices, as of now
# TODO: 
#   - make it more modular, 
#   - integer indices 

# Constructing the modular form $\Delta_{k/2}(Z)$

# The modular form Delta_{k/2,1/2} of weight k/2 and index 1/2 can be obtained as the additive list of the Jacobi form Psi_{k/2,1/2} of \Gamma_1(N) of weight k/2 and index 1/2

q = var('q')
r = var('r')

num_truncate = 10 # Max order of stuff

# Generate the expansion for \theta_1
def theta1_qexpansion(num_truncate=10):
    i = var('i') #Dummy index
    theta1_term = (-1)^i*(r^(i+1/2))*q^(((i+1/2)^2)/2)
    theta1_series = sum(theta1_term, i, -num_truncate, num_truncate)
    return theta1_series

theta1_exp = theta1_qexpansion()

q1 = qexp_eta(ZZ[['q']], num_truncate).truncate(num_truncate)
eta_exp = q1*q^(1/24)

# Weight of psi is k+1/2
# Jacobi form Psi_{k/2,1/2} of \Gamma_1(N) of weight k/2 and index 1/2
def psi_compute(k=10, N=1):
    psi_exp = theta1_exp*eta_exp^((k-4)/2)*(eta_exp(q=q^N))^((k+2)/2)
    return psi_exp.expand()

def psi_coeff(q_power,r_power):
    return psi_exp.coeff(q^(q_power/2)).coeff(r^(r_power/2))

# vim: ft=python
