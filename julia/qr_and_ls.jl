using RDatasets;
using GLM;

longley = dataset("datasets", "longley")
lm(Employed ~ GNPDeflator + GNP + Unemployed + ArmedForces + Population + Year, longley)

y = Array(longley[:Employed])
X = Matrix(longley[:, 2:7])
n = length(y)
X = [ones(n) X]

function house(x::Array{Float64})
    n = length(x)
    sigma = dot(x[2:n]', x[2:n])
    nu = [1; x[2:n]]
    beta = 0
    if sigma == 0
        beta = 0
    else
        mu = sqrt(x[1]^2 + sigma)
        if x[1] <= 0
            nu[1] = x[1] - mu
        else
            
            nu[1] = -sigma / (x[1]+ mu)
        end
        beta = 2nu[1]^2 / (sigma + nu[1]^2)
        nu =nu / nu[1]
        nu, beta
    end
end

function qr_householder(A)
    m, n = size(A) # m::low, n::col
    A= copy(A)
    b= zeros(n)
    for j = 1:n
        nu, beta = house(A[j:m,j])
        b[j]=beta
        A[j:m,j:n] = (eye(m-j+1)- beta*nu*nu')*A[j:m,j:n]
        if j < m
            A[j+1:m, j] = nu[2:m-j+1]
        end
    end
    A, b
end

function qr_q(QR,b)
    m, n = size(QR)
    Q = eye(m)
    for j = n:-1:1
        nu =[1.0; QR[j+1:m,j]]
        Q[j:m, j:m] = (eye(m-j+1) - b[j] * nu * nu') * Q[j:m, j:m]
    end
    Q
end
function qr_r(F)
    F = copy(F)
    m,n= size(F)
    R=zeros(m,n)
    for i = 1:n
        R[i:i,i:n] = F[i:i,i:n]
    end
    R
end

function qr_qtc(F, beta, y)
    m, n = size(F)
    Q = copy(y)
    for j = 1:n
        nu = [1.0; F[j+1:m, j]]
        Q[j:m,:] = (eye(m-j+1) - beta[j]*nu*nu')*Q[j:m,:]
    end
    Q[1:n,:]
end

function back_sub(U,b)
    n,m = size(U)
    b = copy(b)
    b[m]=b[m] / U[m,m]
    for i = m-1:-1:1
        b[i] = (b[i] -dot(U[i,i+1:m],b[i+1:m]))/U[i,i]
    end
    b
end

function qr_coef(F,b, y)
    b=copy(b)
    R_1 = qr_r(F)
    qrqtc = qr_qtc(F,b,y)
    back = back_sub(R_1,qrqtc)
    x = back / qrqtc
    x
end

K=rand(10,3)
nu, beta = house(K)
F, beta = qr_householder(K)
Q = qr_q(F, beta)
R = qr_r(F)
qrqtc = qr_qtc(F,beta,y)
backsub = back_sub(R,qrqtc)
qr_coef(F,beta, y)