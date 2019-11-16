#corresponds to matrix factorizaton , (1), (2)
#(1)simulation
a = rand(4,3); b = rand(3,2)
q1,q2= size(a); p1, p2 = size(b)
#q is item, p is user
r_ui = kron(a,b)
row, col = size(r_ui)


#(2)fit the model
q = rand(q1,q2); p = rand(p1,p2); # random data
# no learning, saved error function for a moment
err_temp = r_ui - kron(q, p)

#sgd
# define object function
# just make square error of each component
lambda = 1e-4
err = Array{Float16}(row, col)
print(err)
for i = 1:row
  for j = 1:col
    se = r_ui[i,j]-q[ceil(Int64,i/p1),ceil(Int64,j/p2)]*p[ceil(Int64,i/q1),ceil(Int64,j/q2)]
    reg = abs2(q[ceil(Int64,i/p1),ceil(Int64,j/p2)])+abs2(p[ceil(Int64,i/q1),ceil(Int64,j/q2)])
    err[i,j] = abs2(se) + lambda*reg
  end
end

#modeling
#https://datajobs.com/data-science-repo/Recommender-Systems-%5BNetflix%5D.pdf
# p.4
# e_ui = r_ui - qi'*pu

epoch = 1000
gamma = 1e-3
q_hat = Array{Float64}(epoch+1, q1, q2)
q_hat[1,:,:] = zeros(Float64,q1,q2)
p_hat = Array{Float64}(epoch+1, p1, p2)
p_hat[1,:,:] = zeros(Float64,p1,p2)
for k = 2:epoch
  for i = 1:row
    for j = 1:col
      #q and p has special inDex because we blocked it
      q_i, q_j = ceil(Int64,i/p1),ceil(Int64,j/p2)
      p_i, p_j = ceil(Int64,i/q1),ceil(Int64,j/q2)
      e_ui = r_ui[i,j] - q_hat[k-1,q_i,q_j]*p_hat[k-1,p_i,p_j]
      q_hat[k,q_i,q_j] = q_hat[k-1,q_i,q_j] + gamma*(e_ui*p_hat[k-1,p_i,p_j] - lambda*q_hat[k-1,q_i,q_j])
      p_hat[k,p_i,p_j] = p_hat[k-1,p_i,p_j] + gamma*(e_ui*q_hat[k-1,q_i,q_j] - lambda*p_hat[k-1,p_i,p_j])
    end
  end
end

#compare the result of the model with simulation data
# print(epoch, lambda, gamma)
# print(r_ui, kron(q_hat[epoch+1,:,:],p_hat[epoch+1,:,:]))
