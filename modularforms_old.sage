# Constructing the modular form $\Delta_{k/2}(Z)$

# The modular form Delta_{k/2,1/2} of weight k/2 and index 1/2 can be obtained as the additive lift of the Jacobi form Psi_{k/2,1/2} of \Gamma_1(N) of weight k/2 and index 1/2

q = var('q')
r = var('r')
s = var('s')

# Remove the fractional exponents
qh = var('qh')
rh = var('rh')
sh = var('sh')

num_truncate = 5 # Max order of eta_expansion

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

########################################
#  Constructing the N=6 Modular Forms  #
########################################


#####################
# The seed is a modular form of the group \Gamma(N*q_num, N)
N=6
q_num=2


# Generate the expansion for \Psi_{1,1/2}
#eta_product_exp = eta_exp.subs(q=q^2)*eta_exp.subs(q=q^3)*eta_exp.subs(q=q^6)/(eta_exp^2)

# Need to divide the eta functions, PowerSeries works better than SymbolicExpression
q1 = qexp_eta(QQ[['q']], num_truncate)
q1_prec = q1.prec()

qt = qexp_eta(QQ[['q']], num_truncate).truncate(num_truncate)

S = PowerSeriesRing(QQ, q)
#q2 = var('q2') 
#q2 = S(q2)# q2 := q^2

q2 = qt.subs(q=q^2); q2 = S(q2).O(2*q1_prec)
q3 = qt.subs(q=q^3); q3 = S(q3).O(3*q1_prec)
q6 = qt.subs(q=q^6); q6 = S(q6).O(6*q1_prec)

qp = q2*q3*q6/q1^2
qp = qp.truncate(num_truncate)
qp = q^(3/8)*qp

eta_product_exp = qp.expand()
###

psi_one_half = theta1_exp*eta_product_exp

psi_one_half = psi_one_half.subs(q=qh^2, r=rh^2)
psi_one_half = psi_one_half.simplify_exp()
psi_one_half = psi_one_half.expand()

# Information about \psi_{1,1/2}
weight = 1
index = 1/2

# HECKE OPERATOR
# Conditions on the sum
#   ad=m; (a,Nq=1), b = 1,...,d-1

# This is what we need. Do not change this.
mu=1

# Need a arbit large number to truncate the series
some_high_number = 5

# Since we have
#   m = mu mod q_num
# mu=1, and q_num=2 just means that m runs over all odd numbers

m_list = range(mu, some_high_number, q_num)

# CHARACTER 
# chi = (-1)^(c/N + b)
# We need chi(sigma_a) and for that b,c=0 mod N*q_num. 
# chi(sigma_a) = (-1)^(q_num*c' + b'*N*q_num)
# Since N,q_num are both even, the character is trivial.

# Compute the additive lift for mu
lift_exp = 0
for m in m_list:

    # Compute the Hecke Operator T_m

    # ad=m  means that a,d are all divisors of m
    m_divisors = divisors(m)
    a_list = [a for a in m_divisors if gcd(a,N*q_num)==1]

    print a_list

    lift_exp_m = 0
    for a in a_list:
        d = m/a
        b_list = range(0,d)

        print a,d,b_list

        for b in b_list:
            lift_exp_m_term = psi_one_half.subs(qh=qh^(a/d)*e^(pi*I*b*q_num/d), rh=rh^a)
            lift_exp_m += lift_exp_m_term

        if not lift_exp_m==0:
            #print(lift_exp_m)
            lift_exp_m = lift_exp_m.simplify_exp()
        
        lift_exp_m *= d^(-weight)

    lift_exp += m^(weight-1)*lift_exp_m*sh^(2*m*index)

lift_exp = lift_exp.simplify_exp()
lift_exp = lift_exp.expand()


#####################
# SQUARED, with k=2

eta_product_exp = (qp^2).expand()
####

phi_two_one = (theta1_exp^2)*eta_product_exp
# Remove the fractional exponents
qh = var('qh')
rh = var('rh')
sh = var('sh')

phi_two_one = phi_two_one.subs(q=qh^2, r=rh^2)
phi_two_one = phi_two_one.simplify_exp()
phi_two_one = phi_two_one.expand()

# vim: ft=python
