using Plots
using LinearAlgebra: norm, eigen, dot, cross
using Printf
##########################################################
##########################################################
function plotpaths(n,r,h2)
#generate and plot n gps
#start paths on small circle of radius r
#use step size h2 for euler's method to calculate paths

f(x,y) = 2.0 - x^2-y^2 #the function
grad(x,y) = [-2x, -2y] #the gradient
h=(2*pi)/n #theta spacing between gps

# Create a grid of x and y values for plotting
x = range(-1.5pi, 1.5pi, length=100)
y = range(-1.5pi, 1.5pi, length=100)

# Define the integration limits
x_min, x_max = -.75*pi, .75*pi
y_min, y_max = -.75*pi, .75*pi


# Evaluate the function on the grid
z = [f(i, j) for i in x, j in y]

# Create the contour plot
plt=contour(x, y, z, levels=20, color=:plasma, xlabel="x", ylabel="y")
# ticks as multiples of pi
xticks!(plt, -pi:pi/2:pi, ["-π", "-π/2", "0", "π/2", "π"])
yticks!(plt, -pi:pi/2:pi, ["-π", "-π/2", "0", "π/2", "π"])
# add h/v lines for the integration limits
vline!(plt, [x_min, x_max], color=:black, linestyle=:dash, label="")
hline!(plt, [y_min, y_max], color=:black, linestyle=:dash, label="")

#allocate array of vectors for paths
gps=[Float64[] for a in 1:2*n]

for i=1:n #generate ith gp
gpx=Float64[]
gpy=Float64[]
start=(cos(h*(i))*r,sin(h*(i))*r)
push!(gpx,start[1])
push!(gpy,start[2])

 while abs(gpx[end])<=x_max && abs(gpy[end])<=y_max
  x0=gpx[end]
  y0=gpy[end]
  step=-h2*grad(x0,y0) #euler's method
  x=gpx[end]+step[1]
  y=gpy[end]+step[2]
  push!(gpx,x)
  push!(gpy,y)
 end

#calculate path end points on boundary
 if gpx[end]>x_max
 gpy[end]=(gpy[end]-gpy[end-1])/(gpx[end]-gpx[end-1])*(x_max-gpx[end-1])+gpy[end-1]
 gpx[end]=x_max
 elseif gpx[end]<x_min
 gpy[end]=(gpy[end]-gpy[end-1])/(gpx[end]-gpx[end-1])*(x_min-gpx[end-1])+gpy[end-1]
 gpx[end]=x_min
 elseif gpy[end]>y_max
 gpx[end]=(gpx[end]-gpx[end-1])/(gpy[end]-gpy[end-1])*(y_max-gpy[end-1])+gpx[end-1]
 gpy[end]=y_max
 elseif gpy[end]<y_min
 gpx[end]=(gpx[end]-gpx[end-1])/(gpy[end]-gpy[end-1])*(y_min-gpy[end-1])+gpx[end-1]
 gpy[end]=y_min
 else
 display("error calculating end point")
 end

gps[2*i-1]=gpx
gps[2*i]=gpy
plot!(gpx,gpy, legend = false)
end

return plt #show plot
end
##########################################################
##########################################################
function getpaths(n,r,h2)
#generate and plot n gps
#start paths on small circle of radius r
#use step size h2 for euler's method to calculate paths

f(x,y) = 2.0 - x^2 - y^2 #the function
grad(x,y) = [-2*x, -2*y] #the gradient
h=(2*pi)/n #theta spacing between gps

# Define the integration limits
x_min, x_max = -.75*pi, .75*pi
y_min, y_max = -.75*pi, .75*pi

#allocate array of vectors for paths
gps=[Float64[] for a in 1:2*n]

for i=1:n #generate ith gp
gpx=Float64[]
gpy=Float64[]
start=(cos(h*(i))*r,sin(h*(i))*r)
push!(gpx,start[1])
push!(gpy,start[2])

 while abs(gpx[end])<=x_max && abs(gpy[end])<=y_max
  x0=gpx[end]
  y0=gpy[end]
  step=-h2*grad(x0,y0) #euler's method
  x=gpx[end]+step[1]
  y=gpy[end]+step[2]
  push!(gpx,x)
  push!(gpy,y)
 end

#calculate end points of paths on boundary
 if gpx[end]>x_max
 gpy[end]=(gpy[end]-gpy[end-1])/(gpx[end]-gpx[end-1])*(x_max-gpx[end-1])+gpy[end-1]
 gpx[end]=x_max
 elseif gpx[end]<x_min
 gpy[end]=(gpy[end]-gpy[end-1])/(gpx[end]-gpx[end-1])*(x_min-gpx[end-1])+gpy[end-1]
 gpx[end]=x_min
 elseif gpy[end]>y_max
 gpx[end]=(gpx[end]-gpx[end-1])/(gpy[end]-gpy[end-1])*(y_max-gpy[end-1])+gpx[end-1]
 gpy[end]=y_max
 elseif gpy[end]<y_min
 gpx[end]=(gpx[end]-gpx[end-1])/(gpy[end]-gpy[end-1])*(y_min-gpy[end-1])+gpx[end-1]
 gpy[end]=y_min
 else
 display("error calculating end point")
 end

gps[2*i-1]=gpx
gps[2*i]=gpy
end

return gps
end
##########################################################
##########################################################
function curv(r)
#calculate the curvature of the isosurface at point r
#taken from Tim's code

 f(x,y) = 2.0 - x^2 - y^2 #the function
 grad(x,y) = [-2*x, -2*y] #the gradient
 hess(x,y) = [-2 0; 0 -2] #the hessian
 
   # Compute gradient and hessian at r
    grad_rho = grad(r[1], r[2])
    hessian_rho = hess(r[1], r[2])

    # Compute the shape operator
    S = - hessian_rho / norm(grad_rho)

    # Compute the eigenvalues and eigenvectors of the shape operator
    e = eigen(S)

    # Sort indices based on the absolute inner products of eigenvectors with gradient (in ascending order)
    sorted_indices = sortperm(1:2, by = i -> abs(dot(grad_rho, e.vectors[:, i])))

    # Select the smallest absolute value inner product eigenvalue and its corresponding eigenvector
    principal_curvature = e.values[sorted_indices[1]]

    return principal_curvature
end
#######################################################
#######################################################
function int_area(n,h,intradius)
#calculate the area of the square using gp-only integration


exact=(2*.75*pi)^2 #exact solution area of square
intarea=(pi)*intradius^2 #area of interior circle where gps start

p=getpaths(n,intradius,h) #get array of gradient paths

#allocate array of vector for l1(s)
arcs=[Float64[] for a in 1:n]

#integrate curvature
for i=1:n #loop over gps
 n2=length(p[2*i])
 arc=Vector{Float64}(undef,n2) #array to hold arc lengths l1(s)
 arc[1]=(1/n)*pi*2*intradius #l1(s0)
  for j=2:n2 #loop over points for gp
   int=0.0
    for k=2:j
      x1=p[2i-1][k-1]
      x2=p[2i-1][k]
      y1=p[2i][k-1]
      y2=p[2i][k]
      h=sqrt((x1-x2)^2+(y1-y2)^2)
      int+=curv((x2,y2))*h #integrate curvature rh rule
    end
   arc[j]=arc[1]*exp(int) #take exponential and multiply by l1(s0)
  end
 arcs[i]=arc
end

#integrate l1(s)
area=0.0
for i=1:n
n2=length(arcs[i])
 for k=1:n2-1
  x1=p[2i-1][k]
  x2=p[2i-1][k+1]
  y1=p[2i][k]
  y2=p[2i][k+1]
  h=sqrt((x2-x1)^2+(y2-y1)^2)
  c1=arcs[i][k]
  c2=arcs[i][k+1]
  area+=h*((c1+c2)/2) #trapezoidal rule
 end
end

area+=intarea #add area of interior circle where gps start
err=abs(area-exact)/exact #calculate relative error
return err
end
#######################################################
#######################################################
function errortablen()
@printf "\n Approximate area of square using gp-only integration"
@printf "\n Use radial function 2-x^2-y^2 for method"
@printf "\n Interior circle radius where gps start: .05"
@printf "\n Use step size h=.001 to generate paths and integrate"
@printf "\n n is the number of paths"
@printf "\n Error is relative error"
@printf "\n EOC is estimated order of convergence"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

n=(2,4,8,16,32,64)
er=Array{Float64}(undef,6)
eoc=Array{Float64}(undef,6)
h=.001
r=.05

er[1]=int_area(n[1],h,r)
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:6
er[i]=int_area(n[i],h,r)
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end

end
#######################################################
#######################################################
function errortableh()
@printf "\n Approximate area of square using gp-only integration"
@printf "\n Use radial function 2-x^2-y^2 for method"
@printf "\n Interior circle radius where gps start: .05"
@printf "\n h is step size used to generate paths and integrate"
@printf "\n n=64 is the number of paths"
@printf "\n Error is relative error"
@printf "\n EOC is estimated order of convergence"
@printf "\n ----------------------------------------------------------------"
@printf "\n h           rel. error         EOC"

h=(.064,.032,.016,.008,.004,.002,.001)
er=Array{Float64}(undef,7)
eoc=Array{Float64}(undef,7)
n=64
r=.05

er[1]=int_area(n,h[1],r)
@printf "\n %g           %g          n/a" h[1] er[1]

for i=2:7
er[i]=int_area(n,h[i],r)
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" h[i] er[i] eoc[i]
end

end
