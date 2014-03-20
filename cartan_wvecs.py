# Author: Hersh Singh [hershdeep@gmail.com]
# Date: August 05, 2013
# Description:
# Given the cartan matrix and the dynkin coefficients of the highest weight, return all the weight vectors, their weights
# Todo: dimensionality of each weight space using freudenthal's formula
# Reference: Cahn Chapter 10

from scipy import *

# Cartan Matrix for the given rep
C = array([[2., -1.], [-1., 2.]]) #SU(3)
#C = array([[2., -1., 0.], [-1., 2., -2.], [0., -1., 2.]]) #B3
N = len(C)

# Dynkin Coeffs for the hightest weight
d_highest = array([1, 0])
#d_highest = array([1, 1]) #SU(3) Adjoint rep
#d_highest = array([0, 0, 1]) #B3

#m[j] = 2<v,alpha[j]>/<alpha[j],alpha[j]>
#M[k] = list of roots at level k

M = [[d_highest]]
Mcoord = [[zeros(N)]]


def get_p(Mcoord, k, i):
    #print "\nin get_p"
    p = zeros(N)
    if k==0:
        return p
    Mc = Mcoord[k][i]
    #print Mc

    # for each dynkin coefficient of the current weight vector
    for n in range(N):
        #print "n=",n
        #for each level above the current level
        #print k-1
        for kk in range(k-1, -1, -1):
            element = Mc + (k-kk)*identity(N)[n]
            #print 'looking for', element, 'in',Mcoord[kk]
            #print matchinlist(Mcoord[kk],element)
            if matchinlist(Mcoord[kk], element):
                p[n]=p[n]+1
            else:
                break
    return p

# Returns true if element is found in the list
def matchinlist(list, element):
    return any([array_equal(e,element) for e in list])

# at level k
k = 0
done_flag = 0
while done_flag == 0:
    print ""
    print "Level:", k
    print "Last row of weight vectors:", M[k]
    print "Last row of weight vectors coords:", Mcoord[k]

    M.append([])
    Mcoord.append([])

    for i, v in enumerate(M[k]):
        print "Weight vector: ",i,v

        p = get_p(Mcoord,k,i)
        m = p+v
        print "M,P,V: ", m,p,v
        if (sum(m>0) == 0):
            done_flag = 1
            break
        v_repeat = tile(v, [sum(m > 0), 1])
        Mcoord_repeat = tile(Mcoord[k][i], [sum(m > 0), 1])
        new_wvecs = v_repeat - C[m > 0]
        # using the fact the True,False is typecasted to 1,0 before doing arithmetic with integers
        new_Mcoord = Mcoord_repeat - identity(N)[m > 0]

        # Clean up by removing duplicates
        #print new_wvecs
        for idx,wvec in enumerate(new_wvecs):
            if not matchinlist(M[k+1],wvec):
                M[k+1].append(wvec)
                Mcoord[k+1].append(new_Mcoord[idx])
    k=k+1
