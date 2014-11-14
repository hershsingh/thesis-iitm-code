def additivelift_int():
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

    

# vim: ft=python
