# This method solves an ODE system, the function returns a list of values
# corresponding to the numeric solution for each value of the independent
# variable given.  It also returns a phase plane plot, only if the system
# provided is a second order system.
#
# Parameters:
# d:        Vector of differential equations
#           Example -> ["Eqn1","Eqn2","Eqn3", ...]
#
# vars:     Vector with the names of the variables including the independet
#           variable at the end.
#           Example -> ["x1","x2","x3", ...,"t"]
#
# init:     Vector with the initial values of each variable including the
#           independet variable at the end.
#           Example -> [x1_0,x2_0,x3_0, ...,t_0]

# h:        Step Size.
#
# lims:     Vector with the limits of the solution.
#           Example -> [x_liminf,x_limsup,y_liminf,y_limsup]
#
# method:   Desired Method to find the solution.
#           Chose Between -> ["euler","midpoint",rk4"]
#
# point:    Desired instant to evaluate the solution.
#
# Pablo Millan and Oscar Neira.




library(phaseR)

system_solve <- function(d,vars,init,h,lims,method,point){

  if (method=="euler"){

    values=matrix(init,ncol = 1)

    for (i in 1:length(d)){
      assign(vars[i],init[i])
    }

    npoints=1
    values=cbind(values,numeric(length(d)+1))

    while (values[1,npoints]>lims[1] && values[2,npoints]>lims[3] &&
           values[1,npoints]<lims[2] && values[2,npoints]<lims[4]) {

      for (i in 1:length(d)){
        values[i,npoints+1]=values[i,npoints]+h*eval(parse(text=d[i]))
      }

      for (i in 1:length(d)){
        assign(vars[i],values[i,npoints+1])
      }

      values[length(d)+1,npoints+1]=values[length(d)+1,npoints]+h

      assign(vars[length(d)+1],values[length(d)+1,npoints+1])

      values=cbind(values,numeric(length(d)+1))
      npoints=npoints+1

    }

    values=values[,-npoints-1]
  }


  if (method=="midpoint"){

    values=matrix(init,ncol = 1)

    for (i in 1:length(d)){
      assign(vars[i],init[i]+h/2)
    }

    npoints=1
    values=cbind(values,numeric(length(d)+1))

    while (values[1,npoints]>lims[1] && values[2,npoints]>lims[3] &&
           values[1,npoints]<lims[2] && values[2,npoints]<lims[4]) {

      for (i in 1:length(d)){
        values[i,npoints+1]=values[i,npoints]+h*eval(parse(text=d[i]))
      }

      for (i in 1:length(d)){
        assign(vars[i],values[i,npoints+1]+h/2)
      }

      values[length(d)+1,npoints+1]=values[length(d)+1,npoints]+h

      assign(vars[length(d)+1],values[length(d)+1,npoints+1])

      values=cbind(values,numeric(length(d)+1))
      npoints=npoints+1

    }

    values=values[,-npoints-1]
  }

  if (method=="rk4"){

    values=matrix(init,ncol = 1)

    for (i in 1:length(d)){
      assign(vars[i],init[i]+h/2)
    }

    npoints=1
    values=cbind(values,numeric(length(d)+1))

    while (values[1,npoints]>lims[1] && values[2,npoints]>lims[3] &&
           values[1,npoints]<lims[2] && values[2,npoints]<lims[4]) {

      k1=numeric(length(d))
      k2=numeric(length(d))
      k3=numeric(length(d))
      k4=numeric(length(d))

      for (i in 1:length(d)){
        k1[i]=h*eval(parse(text=d[i]))
        assign(vars[i],values[i,npoints]+k1[i]/2)
      }

      for (i in 1:length(d)){
        k2[i]=h*eval(parse(text=d[i]))
        assign(vars[i],values[i,npoints]+k2[i]/2)
      }

      for (i in 1:length(d)){
        k3[i]=h*eval(parse(text=d[i]))
        assign(vars[i],values[i,npoints]+k3[i])
      }

      for (i in 1:length(d)){
        k4[i]=h*eval(parse(text=d[i]))
      }

      for (i in 1:length(d)){
        values[i,npoints+1]=values[i,npoints]+ k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6
      }

      for (i in 1:length(d)){
        assign(vars[i],values[i,npoints+1])
      }

      values[length(d)+1,npoints+1]=values[length(d)+1,npoints]+h
      assign(vars[length(d)+1],values[length(d)+1,npoints+1])

      values=cbind(values,numeric(length(d)+1))
      npoints=npoints+1

    }

    values=values[,-npoints-1]
  }

  if (length(d)==2){

    Field <- function(t, y_in, parameters){

      dy = numeric(length(y_in))

      for (i in 1:length(dy)){
        for (j in 1:length(dy)){
          assign(vars[j],y_in[j])
        }
        dy[i]=eval(parse(text=parameters[i]))
      }
      return(list(dy))

    }

    Field.flowField <- flowField(Field, xlim=c(lims[1],lims[2]),
                                 ylim = c(lims[3],lims[4]), parameters = d,
                                 points = 20,
                                 add = FALSE, xlab = vars[1], ylab = vars[2] ,
                                 main = "Solution")

    grid()
    lines(values[1,],values[2,], col='red', type = 'l',xlim=c(lims[1],lims[2]),ylim=c(lims[3],lims[4]))
  }

  aux=point-values[3,]
  index=which.min(abs(aux))
  solution_point=values[,index]

  solution <- list("trajectory" = values, "answer" = solution_point)

  return(solution)
}

#sol=system_solve(c("x*y","x*y-1"),c("x","y","t"),c(-2,1,0),0.1,c(-4,4,-3,3),"euler",0.62)
#print(sol)

#sol=f(c("x*y","x*y-1"),c("x","y","t"),c(-2,1,0),0.01,c(-4,4,-3,3),"midpoint",0.5)
#print(sol)

#sol=system_solve(c("x*y","x*y-1"),c("x","y","t"),c(-2,2,0),0.01,c(-4,4,-3,3),"rk4",0.674)
#print(sol)

#sol=f(c("x-y+z","x*y-1","z+3*x"),c("x","y","z","t"),c(0,1,1,0),0.01,c(-6,6,-4,4),"euler",0.5)
#print(sol)
