using Gadfly

#simulation
n, m = 30, 3
x = [ones(n) rand(n, m)]
e = 2*rand(n,1)
b = [i for i in 1:2:8]
Y = x*b +e

#sgd
x_sgd = copy(x)
y_sgd = copy(Y)

bhat_sgd= Array{Float64}(iter+1, m+1)
bhat_sgd[1, :] = [0; 0; 0; 0]
gamma = 1e-5
row , col = size(x_sgd)
iter = 1000
tr = Int64[]
for i = 1:iter
  temp = rand(1:row,1)
  for j =1:col
    bhat_sgd[i+1,j]=bhat_sgd[i,j] - gamma * dot(x_sgd[temp,j],(x_sgd[temp,j]*bhat_sgd[i,j]-y_sgd[temp]))
  end
  x_sgd=x_sgd[shuffle(1:row),:]
  append!(tr,temp)
end
bhat_sgd[iter+1,:]
b
#plot coefs
plot(layer(x=1:iter, y=bhat_sgd[2:end,1], Geom.line, Theme(default_color=color("yellow"))),
  layer(x=1:iter, y=bhat_sgd[2:end,2], Geom.line, Theme(default_color=color("blue"))),
  layer(x=1:iter, y=bhat_sgd[2:end,3], Geom.line, Theme(default_color=color("green"))),
  layer(x=1:iter, y=bhat_sgd[2:end,4], Geom.line, Theme(default_color=color("red")))
)

#plot obj fuc
S = rand(iter,1)
for i = 2:iter
  err = x_sgd[tr[i],:]*bhat_sgd[i]-y_sgd[tr[i]]
  S[i,1] = S[i-1,1]+dot(err,err)
end
plot(x=1:iter, y=S, Geom.line)
