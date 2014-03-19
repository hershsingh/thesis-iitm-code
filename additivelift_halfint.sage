# Given a Jacobi form of weight k and index 1, the Maass lift is
#   \Phi_k(Z) = \sum_{(n,l,n)>0} \sum_{d|(n,l,m)} d^{k-1} f(nm/d^2, l/d) q^n r^l s^m

q, r, s = var('q r s')

def character(d):
    return 1

def additivelift_halfint():
    Delta=0
    for n in range(1,n_max,2):
        for m in range(1,n+1,2):
            print(n,m) #For debugging

            g1 = gcd(n,m)

            lmax = floor(2*sqrt(n*m))
            l_list = range(-lmax,lmax+1)

            for l in l_list:
                g = gcd(g1,l)
                divisors_g = divisors(g)
                #print "l=", l,divisors_g
                for d in divisors_g:
                    fc = psi_coeff((n*m)/(d^2), l/d)
                    #print d,divisors_g
                    if not fc==0:
                        term_coeff = d^((k-2)/2)*fc
                        if n==m:
                            term_vars = q^(n/2)*r^(l/2)*s^(m/2)
                        else:
                            term_vars = q^(n/2)*r^(l/2)*s^(m/2) + q^(m/2)*r^(l/2)*s^(n/2)
                        Delta += character(d)*term_coeff*term_vars
                        #print character(d)*term_coeff*term_vars
    return Delta

    
def Delta_coeff(Delta_expansion, q_power,s_power):
    # q_power and s_power are (half) integers
    return Delta_expansion.coeff(q^q_power).coeff(s^s_power)

def Delta_latex(Delta_expansion):
    Delta_latex = 0
    for n in range(3,n_max,2):
        Delta_latex += Delta_coeff(Delta_expansion,n/2,n/2)*(q^(n/2)*s^(n/2))
    for n in range(3,n_max,2):
        Delta_latex += Delta_coeff(Delta_expansion,n/2,1/2)*(q^(n/2)*s^(1/2) + q^(1/2)*s^(n/2))
    return Delta_latex

# Global variables (almost)
n_max=7
k=10
N=1

psi_exp = psi_compute(k=10, N=1)
Delta_5 = additivelift_halfint()
Delta_5_tex = Delta_latex(Delta_5)

# vim: ft=python
