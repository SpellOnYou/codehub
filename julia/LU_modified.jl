# GramISchmidt Orthogonalization
function qrmgs(A)
  V,Q = copy(A),copy(A)
  m, n = size(A)
  for i in 1:n
    Q[1:m,i] = V[1:m,i]/norm(V[1:m,i])
    for j in i+1:n
      rij = dot(Q[1:m,i],V[1:m,j])
      V[1:m,j] -= rij * Q[1:m,i]
    end
  end
  Q
end

result2 = qrmgs(A)

function orthocheck(x)
  A = copy(x)
  for j = 1:size(A)[2]
    for i = 1:length(A[1,:])
      # if i != j
        @printf "%s dot %s is %f \n" i j dot(A[1:4,i],A[1:4,j])
      # end
    end
  end
end
orthocheck(result2)
