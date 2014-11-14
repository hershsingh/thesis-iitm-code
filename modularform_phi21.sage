# CYCLE SHAPE 1^2.2^2.3^2.6^2

# Parameters
num_truncate = 15 #Truncate the expansions num_truncate
n_max = 10 #Truncate the summation for Lift at n_max

q = var('q')
r = var('r')
s = var('s')

# To remove the fractional exponents
qh = var('q_h')
rh = var('r_h')
sh = var('s_h')

#######################
#  Generate the seed  #
#######################

# Generate the expansion for \theta_1
def theta1_qexpansion(num_truncate=10):
    i = var('i') #Dummy index
    theta1_term = (-1)^i*(r^(i+1/2))*q^(((i+1/2)^2)/2)
    theta1_series = sum(theta1_term, i, -num_truncate, num_truncate)
    return theta1_series

# The seed is a modular form of the group \Gamma_0(N)
N=6

# Need to divide the eta functions, PowerSeries works better than
# SymbolicExpression as it is faster and also keeps track of the order
# of the series.
q1 = qexp_eta(QQ[['q']], num_truncate)
q1_prec = q1.prec()
qt = qexp_eta(QQ[['q']], num_truncate).truncate(num_truncate)

# Power Series Ring over the field of Rational Numbers QQ in the
# variable q
S = PowerSeriesRing(QQ, q)

q2 = qt.subs(q=q^2); q2 = S(q2).O(2*q1_prec)
q3 = qt.subs(q=q^3); q3 = S(q3).O(3*q1_prec)
q6 = qt.subs(q=q^6); q6 = S(q6).O(6*q1_prec)

qpp = q2*q3*q6/q1^2

qp = qpp.truncate(num_truncate)
qp = q^(3/8)*qp

eta_product_exp = qp.expand()
eta_product_prec = qpp.prec()
###

# Generate the theta series only uptil the order of the eta product
# expansion.  The higher terms are useless, and it seems to be
# difficult to truncate a symbolic expression later.
theta1_num_truncate = ceil((sqrt(3*eta_product_prec)-1)/2)
theta1_exp = theta1_qexpansion(num_truncate=theta1_num_truncate)
theta1_exp = theta1_exp.simplify_exp().expand()

# This is valid only upto order
# Min(Order(theta1_exp^2),Order(eta_product_exp))^2.  Need to find a
# way to truncate it, to get rid of the higher order junk
phi_two_one = (theta1_exp*eta_product_exp)^2

phi_two_one = phi_two_one.simplify_exp()
phi_two_one = phi_two_one.expand()

# Information about \phi_{2,1}
weight = 2
index = 1

###################
#  Additive Lift  #
###################

# The additive lift in this case is straightforward and the character
# is trivial.  
# Reference: [S. Govindarajan] [2009] [arXiv: 0907.1410] Equation 3.22 
# Phi_2(Z) = \sum_{(n,l,m)>0} \sum_{d|(n,l,m) \\ d=1,5 mod 6} d^(k-1)
# a(mn/d^2,l/d) q^n r^l s^m 
# where (n,l,m)>0 means n,m > 0 and 4nm-l^2>0

def phi_two_one_coeff(m,n):
    if m==0:
        x=flatten(phi_two_one.coeffs(q))
        for i, exp in enumerate(x[1::2]):
            if exp==0:
                return x[2*i].coeff(r^n)
    if n==0:
        x=flatten(phi_two_one.coeffs(r))
        for i, exp in enumerate(x[1::2]):
            if exp==0:
                return x[2*i].coeff(q^m)
        
    return phi_two_one.coeff(q^m).coeff(r^n)

def additive_lift(n_max = 10):

    Phi=0
    for n in range(1,n_max+1,1):
        for m in range(1,n+1,1):

            gcd_mn = gcd(n,m)

            lmax = floor(2*sqrt(n*m))
            l_list = range(-lmax,lmax+1)

            for l in l_list:

                g = gcd(gcd_mn,l)
                d_list = [d for d in divisors(g) \
                        if d in union(range(1,g+1,6),range(5,g+1,6))]

                for d in d_list:

                    fc = phi_two_one_coeff((n*m)/(d^2), l/d)

                    if not fc==0:
                        term_coeff = d^(weight-1)*fc
                        if (n==m):
                            term_vars = q^(n)*r^(l)*s^(m)
                        else:
                            term_vars = \
                                    q^(n)*r^(l)*s^(m) + q^(m)*r^(l)*s^(n)
                        Phi += term_coeff*term_vars
    return Phi

Phi_2 = additive_lift(n_max)

# vim: ft=python
