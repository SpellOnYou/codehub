using Gadfly, StatsBase
#1.	Generate a random data set
#1-(1)
n, m = 8, 3
x = [ones(n) rand(n, m)]
e = rand(n,1)
b = [i for i in 1:2:8]
Y = x*b + e
#1-(2)
n2, m2 = 100, 5
x2 = [ones(n2) rand(n2, m2)]
e2 = rand(n2,1); b2 = sample(1:20,6)
Y2 = x2*b2 + e2

#2. Compute bhat : least squares estimates. Use QR.
#2-(1)
q,r = qr(x)
bhat_qr = inv(r)*q'Y
#2-(2)
q2,r2 = qr(x2)
bhat2_qr = inv(r2)*q2'Y2
