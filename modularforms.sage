########################################
#  Constructing the N=6 Modular Forms  #
########################################

num_truncate=20
 #This is the where the Hecke Operators T_m go upto
m_max= 10

q = var('q')
r = var('r')
s = var('s')

# To remove the fractional exponents
qh = var('q_h')
rh = var('r_h')
sh = var('s_h')

# Generate the expansion for \theta_1
def theta1_qexpansion(num_truncate=10):
    i = var('i') #Dummy index
    theta1_term = (-1)^i*(r^(i+1/2))*q^(((i+1/2)^2)/2)
    theta1_series = sum(theta1_term, i, -num_truncate, num_truncate)
    return theta1_series

#####################
# The seed is a modular form of the group \Gamma(N*q_num, N)
N=6
q_num=2

# Need to divide the eta functions, PowerSeries works better than SymbolicExpression as it is faster and also keeps track of the order of the series.
q1 = qexp_eta(QQ[['q']], num_truncate)
q1_prec = q1.prec()
qt = qexp_eta(QQ[['q']], num_truncate).truncate(num_truncate)

# Power Series Ring over the field of Rational Numbers QQ in the variable q
S = PowerSeriesRing(QQ, q)

q2 = qt.subs(q=q^2); q2 = S(q2).O(2*q1_prec)
q3 = qt.subs(q=q^3); q3 = S(q3).O(3*q1_prec)
q6 = qt.subs(q=q^6); q6 = S(q6).O(6*q1_prec)

qpp = q2*q3*q6/q1^2
qp = qpp.truncate(num_truncate)
qp = q^(3/8)*qp
#qo = var('qo')
#qop = qp.subs(q=qo^8).simplify_exp()


eta_product_exp = qp.expand()
eta_product_prec = qpp.prec()
###

# Generate the theta series only uptil the order of the eta product expansion. The higher terms are useless, and it seems to difficult to truncate a symbolic expression later.
theta1_num_truncate = ceil((sqrt(3*eta_product_prec)-1)/2)
theta1_exp = theta1_qexpansion(num_truncate=theta1_num_truncate)
theta1_exp = theta1_exp.simplify_exp().expand()

#theta1_nfexp = theta1_exp.subs(q=qo^8, r=rh^2).simplify_exp().expand() #No fractions expansion


# This is valid only upto order Min(Order(theta1_exp),Order(eta_product_exp))^2. 
# Need to find a way to truncate it, to get rid of the higher order junk
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

 #This is the where the Hecke Operators T_m go upto
#m_max= 5

# Since we have
#   m = mu mod q_num
# mu=1, and q_num=2 just means that m runs over all odd numbers

m_list = range(mu, m_max+1, q_num)

# CHARACTER 
# From the table of characters, chi = (-1)^(c/N + b)
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
            lift_exp_m_term = psi_one_half.subs(q_h=q_h^(a/d)*e^(pi*I*b*q_num/d), r_h=r_h^a)
            lift_exp_m += lift_exp_m_term

        if not lift_exp_m==0:
            #print(lift_exp_m)
            lift_exp_m = lift_exp_m.simplify_exp()
        
        lift_exp_m *= d^(-weight)

    lift_exp += m^(weight-1)*lift_exp_m*sh^(2*m*index)

lift_exp = lift_exp.simplify_exp()
lift_exp = lift_exp.expand()

Phi_2p = (lift_exp^2).subs(q_h=q^(1/2), r_h=r^(1/2), s_h=s^(1/2)).simplify_exp().expand()

# vim: ft=python
