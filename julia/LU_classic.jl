# GramISchmidt Orthogonalization
A = rand(4,3)
function qrcgs(A)
  A = copy(A)
  m, n = size(A)
  B = Array{Float64}(m,n)
  for j in 1:n
  vj = A[1:m,j]
    for i in 1:j-1
      rij = dot(B[1:m,i]',vj)
      vj -= rij * B[1:m,i]
    end
  rjj = norm(vj)
  B[1:m,j] = vj/rjj
  end # col end
  B
end #function end

result1 = qrcgs(A)
dot(result1[1:4,1],result1[1:4,2])

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
orthocheck(result1)