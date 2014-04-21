N=24

plists = Partitions(N).list()

num_truncate = 100 # Max order of stuff

# The eta product expansion is valid only till order num_truncate*min(plist)
n_max = num_truncate*min(plist)

# Generate a list of all pairs of coprimes such that their product is less than n_max
coprimes = [[x,y] for x in range(1,n_max) for y in range(1,x) if gcd(x,y)==1 and x*y<n_max]


# Eta expansion
q = var('q')
q1 = qexp_eta(ZZ[['q']], num_truncate).truncate(num_truncate)
eta_exp = q1*q^(1/24) 

# Get the coeff
def etacoeff(i):
    return etaproduct.coeff(q^(i))

# Number of possible multiplicative functions
num_multiplicative = 0 

# Loop over all balanced cycles
for plist in balanced_cycles:

    # Generate the eta product expansion, valid only till order num_truncate*min(plist)
    etaproduct = 1
    for num in plist:
        etaproduct = etaproduct * eta_exp.subs(q=q^num)
    etaproduct = etaproduct.simplify_exp().expand()

    # Check the first few coefficients
    is_multplicative=1
    for [x,y] in coprimes:
        if not etacoeff(x)*etacoeff(y) == etacoeff(x*y):
            is_multplicative=0
            break;

    if is_multplicative==1:
        print "Could be..", plist
        num_multiplicative += 1
    else:
        print "\tNope!", plist

print "No of possible multiplicative functions =", num_multiplicative

# vim: ft=python
