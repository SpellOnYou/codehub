using Gadfly, StatsBase
#1.	Generate a random data set
#1-A
n, m = 30, 3
x = [ones(n) rand(n, m)]
e = 2*rand(n,1) #+[]
b = [i for i in 1:2:8]
Y = x*b + e

#3.	Implement Gradient Descent for the least squares estimation.
#3-(1)
x_gd= copy(x)
y_gd = copy(Y)
r = 10000
bhat_gd=Array{Float64}(r+1,m+1) # +1 is for the intercept
bhat_gd[1,:]=[0;0;0;0]
gamma = 0.01#1e-4#
for i = 1:r
  bhat_gd[i+1,:] =bhat_gd[i,:] - gamma * x_gd' *(x_gd*bhat_gd[i,:]-y_gd)
end
bhat_gd[r+1,:]
b
#3-(2) plot the trace of the coefs
plot(layer(x=1:r, y=bhat_gd[2:end,1], Geom.line, Theme(default_color=color("yellow"))),
  layer(x=1:r, y=bhat_gd[2:end,2], Geom.line, Theme(default_color=color("blue"))),
  layer(x=1:r, y=bhat_gd[2:end,3], Geom.line, Theme(default_color=color("green"))),
  layer(x=1:r, y=bhat_gd[2:end,4], Geom.line, Theme(default_color=color("red")))
)
b
#3-(3) plot the objective function
S = rand(r,1)
for i = 2:r
  e = x*bhat_gd[i,:]-y_gd
  S[i,1] = S[i-1,1]+dot(e',e)
end
plot(x=1:r, y=S, Geom.line)
S[r,1]


#1-B
n2, m2 = 100, 5
x2 = [ones(n2) rand(n2, m2)]
e2 = rand(n2,1); b2 = [11;15;1;6;14;7]#sample(1:20,6)
Y2 = x2*b2 + e2

#3.	Implement Gradient Descent for the least squares estimation.
#3-(1)
x2_gd = copy(x2); y2_gd = copy(Y2); r2 = 10000 #r2=1000
bhat_gd2=Array{Float64}(r2+1,m2+1)
bhat_gd2[1,:]=[0;0;0;0;0;0]; gamma2 = 0.1 #1e-4 #gamma2=0.1
for i = 1:r2
  bhat_gd2[i+1,:] = bhat_gd2[i,:] - gamma2 * x2_gd' *(x2_gd*bhat_gd2[i,:]-y2_gd)
end
#"cyan" "magenta" "yellow"
#3-(2) plot the trace of the coefs
plot(layer(x=1:r2, y=bhat_gd2[2:end,1], Geom.line, Theme(default_color=color("red"))),
  layer(x=1:r2, y=bhat_gd2[2:end,2], Geom.line, Theme(default_color=color("blue"))),
  layer(x=1:r2, y=bhat_gd2[2:end,3], Geom.line, Theme(default_color=color("green"))),
  layer(x=1:r2, y=bhat_gd2[2:end,4], Geom.line, Theme(default_color=color("orange"))),
  layer(x=1:r2, y=bhat_gd2[2:end,4], Geom.line, Theme(default_color=color("cyan"))),
  layer(x=1:r2, y=bhat_gd2[2:end,4], Geom.line, Theme(default_color=color("yellow")))
)
#3-(3) plot the objective function
S2 = rand(r2,1)
for i = 2:r2
  se = x2*bhat_gd2[i,:]-y2_gd
  S2[i,1] = S2[i-1,1]+dot(se',se)
end
plot(x=1:r2, y=S2, Geom.line)
S2[r2,1]
