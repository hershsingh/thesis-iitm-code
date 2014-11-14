q1 = qexp_eta(QQ[['q']], num_truncate)
qt = qexp_eta(QQ[['q']], num_truncate).truncate(num_truncate)

S = PowerSeriesRing(QQ, q)
q2 = S(q2)

q2 = qt.subs(q=q^2); q2 = S(q2)
q3 = qt.subs(q=q^3); q3 = S(q3)
q6 = qt.subs(q=q^6); q6 = S(q6)

qp = q2*q3*q6/q1^2
qp = qp.truncate(num_truncate)
qp = q^(1/2)*qp

eta_product_exp = qp
