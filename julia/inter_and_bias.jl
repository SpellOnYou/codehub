#corresponds to matrix factorizaton , (3), (4)
#(1)simulation; b_ui & r_ui

#(1)-2 bias, b_ui
mu=2.0; srand(123);
user = randn(12); item = randn(6)
b_ui = [i+j+mu for i in item for j in user]
b_ui=reshape(b_ui,(12,6))
b1 = length(user); typeof(b1)
b2 = length(item)

#(1)-1 interation, r_ui
q = rand(4,3); p = rand(3,2)
r_ui = kron(q,p)
q1, q2 = size(q); p1,p2 = size(p)

#(2) fit the model

lambda = 1e-4
gamma = 1e-3
epoch = 1000
q_hat=Array{Float64}(epoch+1,q1,q2); q_hat[1,:,:] = zeros(Float64,q1,q2)
p_hat=Array{Float64}(epoch+1,p1,p2); p_hat[1,:,:] = zeros(Float64,p1,p2)
#bi is item and bu is user
bi_hat=Array{Float64}(epoch+1,1,b2)
bi_hat[1,:,:] = zeros(Float64,1,b2)
bu_hat=Array{Float64}(epoch+1,b1,1)
bu_hat[1,:,:] = zeros(Float64,b1,1)
row, col = size(r_ui)
kok = NaN
for k = 2:epoch
  for i = 1:row
    for j = 1:col
        p_i, p_j = ceil(Int64,i/q1),ceil(Int64,j/q2)
        q_i, q_j = ceil(Int64,i/p1),ceil(Int64,j/p2)
        e_ui = r_ui[i,j]-mu-bi_hat[k-1,:,j]-bu_hat[k-1,i,:]-q_hat[k-1,q_i,q_j]*p_hat[k-1,p_i,p_j]
        kok = e_ui
        bi_hat[k,:,j] = bi_hat[k-1,:,j] + gamma*(e_ui- lambda*bi_hat[k-1,:,j])
        bu_hat[k,i,:] = bu_hat[k-1,i,:] + gamma*(e_ui- lambda*bu_hat[k-1,i,:])
        q_hat[k,q_i,q_j] = q_hat[k-1,q_i,q_j] + gamma*(e_ui[1]*p_hat[k-1,p_i,p_j] - lambda*q_hat[k-1,q_i,q_j])
        p_hat[k,p_i,p_j] = p_hat[k-1,p_i,p_j] + gamma*(e_ui[1]*q_hat[k-1,q_i,q_j] - lambda*p_hat[k-1,p_i,p_j])
    end
  end
end

#compare the result of the model with simulation data
# print(epoch, lambda, gamma)
# print(r_ui, kron(q_hat[epoch+1,:,:],p_hat[epoch+1,:,:]))
# print(user, bu_hat[epoch+1,:,:])
# print(item, bi_hat[epoch+1,:,:])
