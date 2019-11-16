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
        vec(nu), beta
    end
end

function house_bidiag(A)
    A = copy(A)
    m,n = size(A)
    for j = 1:n
        nu, beta = house(A[j:m,j])
        A[j:m,j:n] = (eye(m-j+1)- beta*nu*nu')*A[j:m,j:n]
        A[j+1:m,j] = nu[2:m-j+1]
        if j <= n-2
            nu, beta = house(A[j,j+1:n]')
            A[j:m,j+1:n] = A[j:m,j+1:n] * (eye(n-j)-beta*nu*nu')
            A[j,j+2:n] = nu[2:n-j]'
        end
    end
    Bidiagonal(A,true)
end

function givens_local(a,b)
    c= 0; s=1
    if b == 0
        c =1
        s=0
    else
        if abs(b) > abs(a)
            tau = -a/b
            s=1/sqrt(1+tau^2)
            c=s*tau
        else
            tau = -b/a
            c=1/sqrt(1+tau^2)
            s=c*tau
        end
    end
    [c s;-s c]
end

function svd_local(B)
  B=copy(B)
  m,n=size(B)
  T=B' * B
  d = (T[n-1,n-1]-T[n,n])/2
  mu = T[n,n] - T[n,n-1]^2 /(d+sign(d)*sqrt(d^2+T[n,n-1]^2))

  y = T[1,1] -mu
  z = T[1,2]
  for k = 1:n-1
    B[:,k:k+1] = B[:,k:k+1] * givens2(y,z)
    y=B[k,k]
    z=B[k+1,k]
    B[k:k+1,:] = givens2(y,z)' * B[k:k+1,:]
    if k < n-1
      y = B[k,k+1]; z=B[k,k+2]
    end
  end
  B
end
# svd_local(rand(10,7))
