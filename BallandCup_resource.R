#This clears everything from the working environment.
#I always like to make sure I start off clean, in case
#I opened this .R file in an R session that was already
#doing other things.
rm(list=ls())


#PART ONE

#This section checks if packages are downloaded, downloads if necessary,
#and then loads required packages.
{
  if (!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
  }
  if (!require(gganimate)){
    install.packages("gganimate")
    library(gganimate)
  }
  if (!require(demodelr)){
    install.packages("demodelr")
    library(demodelr)
  }
  if (!require(pracma)){
    install.packages("pracma")
    library(pracma)
  }
  if (!require(polynom)){
    install.packages("polynom")
    library(polynom)
  }
  if (!require(gifski)){
    install.packages("gifski")
    library(gifski)
  }
}

#Keep these bits defining the functions collapsed
#unless you want to see how it works. To collapse
#or uncollapse them, hit the little triangle on line 35
{
#Set a color palette
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  
#Potential function is a fourth order polynomial
#a0+a1*x+a2*x^2+a3*x^3+a4*x^4
#Derivative is a1+2*a2*x+3*a3*x^2+4*a4*x^3

simulator=function(basinleft,threshold,basinright,noise_magnitude,inits,duration){
#Solve the polynomial to get the basins and threshold in
#the right place
a_parms=as.numeric(poly.calc(c(basinleft,threshold,basinright)))

derivative <- c(dx ~ -a1-a2*x-a3*x^2-a4*x^3)
stochastic_func <-  c(dx ~ 1)

# Identify the parameters
parameters <- c(a1=a_parms[1],a2=a_parms[2],a3=a_parms[3],a4=a_parms[4])   # parameters: a named vector

# Identify how long we run the simulation
deltaT_func <- .0001    # timestep length
timesteps_func <- duration/(365*deltaT_func)   # must be a number greater than 1

# Identify the standard deviation of the stochastic noise
D_func <- noise_magnitude/deltaT_func

# Do one simulation of this differential equation
stoch_out_list=list()
for(i in 1:length(inits)){
  stoch_out_list[[i]]<- euler_stochastic(
    deterministic_rate = derivative,
    stochastic_rate = stochastic_func,
    initial_condition = inits[i],
    parameters = parameters,
    deltaT = deltaT_func,
    n_steps = timesteps_func,
    D = D_func
  )
}


potential_function=function(N,a_parms){
  a_parms[1]*N+a_parms[2]/2*N^2+a_parms[3]/3*N^3+a_parms[4]/4*N^4
}

letters=c("A","B","C","D","E","F","G","H","I","J")

sim_df=data.frame(day=rep(NA,length(stoch_out_list[[1]]$t)*length(inits)),
                  N=rep(NA,length(stoch_out_list[[1]]$t)*length(inits)),
                  simID=rep(NA,length(stoch_out_list[[1]]$t)*length(inits)))


for(i in 1:length(inits)){
  sim_df$day[(1+(i-1)*length(stoch_out_list[[1]]$t)):(i*length(stoch_out_list[[1]]$t))]=365*stoch_out_list[[i]]$t
  sim_df$N[(1+(i-1)*length(stoch_out_list[[1]]$t)):(i*length(stoch_out_list[[1]]$t))]=stoch_out_list[[i]]$x
  sim_df$simID[(1+(i-1)*length(stoch_out_list[[1]]$t)):(i*length(stoch_out_list[[1]]$t))]=letters[i]
}

sim_df$Potential=potential_function(sim_df$N,a_parms)

full_range=range(c(sim_df$N,basinleft,basinright))

surface_df=data.frame(N=linspace(full_range[1]-(diff(full_range))/6,full_range[2]+(diff(full_range))/6,1e3))
surface_df$Potential=potential_function(surface_df$N,a_parms)
  output=list()
  output[[1]]=sim_df
  output[[2]]=surface_df
  return(output)
}

landscape_plotter<-function(simulation){
  sim_df<-simulation[[1]]
  surface_df<-simulation[[2]]
  graph0 <- ggplot(surface_df,aes(x=N,y=Potential))+
    geom_line()+
    theme_classic()+
    theme(legend.position="none",axis.title=element_text(size=25),
          axis.text=element_text(size=20),
          plot.subtitle=element_text(size=25),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  graph0
}

landscape_gifer<-function(simulation){
  sim_df<-simulation[[1]]
  surface_df<-simulation[[2]]
graph1 <- ggplot(sim_df,aes(x=N,y=Potential,color=simID))+
  geom_point(alpha=0.7,stroke=0,size=16)+
  scale_colour_manual(values=cbPalette)+
  theme_classic()+
  theme(legend.position="none",axis.title=element_text(size=25),
        axis.text=element_text(size=20),
        plot.subtitle=element_text(size=25),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


graph1.animation = graph1 +
  transition_time(day) +
  labs(subtitle = "Day: {round(frame_time)}") +
  shadow_wake(wake_length = 0.1)+
  geom_line(data=surface_df,aes(x=N,y=Potential),color="black")

#Takes about ten seconds for the gif to come up
animate(graph1.animation, height = 500, width = 600, fps = 10, duration = 20,
        end_pause = 10, res = 50, renderer=gifski_renderer())
}

timeseries_plotter<-function(simulation){
  sim_df<-simulation[[1]]
  graph2 <- ggplot(sim_df,aes(x=day,y=N,color=simID))+
    geom_line(linewidth=2)+
    theme_classic()+
    theme(legend.position="none",axis.title=element_text(size=25),
          axis.text=element_text(size=20),
          plot.subtitle=element_text(size=25))+
    scale_colour_manual(values=cbPalette)
  graph2
}
}

#This code runs some number of stochastic differential equations
#with each simulation able to have its own starting conditions and
#its own instance of the stochasticity. Then it lets you plot things
#in useful ways.

#simulator function arguments
#basinleft sets the N (x-axis) position of the left basin.
#threshold sets the N position of the threshold dividing the two basins
#basinright sets the N position of the right basin. This should be more
#than basinleft.

#For basinleft, threshold, and basinright, you can only move these things
#left by decreasing the value or right by increasing the value. Doing so
#will lead to correlated changes in curvature, potential depth, etc.

#noise_magnitude should be non-negative and determines how
#strong stochasticity is in these simulations. If you make this
#SUPER big, the simulator may crash or do crazy things.
#inits gives a list of initial conditions
#duration sets how long the simulation runs in days (though all of the units
#and numbers here are totally theoretical)

#Takes about one second to simulate, for these parameter values
##they have the same initial conditions and parameter values but have different stochasticity (random noise)
##noise magnitude = the standard deviation of the noise
simulation1 <- simulator(basinleft = 60,
                       threshold = 90,
                       basinright = 130,
                       noise_magnitude = 2,
                       inits = c(x = 60,x = 60,x = 60),
                       duration = 20)

landscape_plotter(simulation1)
timeseries_plotter(simulation1)
#Takes about thirty seconds to make the gif for these parameter values
landscape_gifer(simulation1)

anim_save("simulation1.gif")

#You can run your own example by changing the inputs
#to the simulator() function. The number of simulations
#run at once depends on how many initial conditions you provide (max 10)

#In general, copy line 177 and rename simulation1 to something else so you
#can explore how this works without losing simulation1.

##low resilience system same initial conditions
simulation2 <- simulator(basinleft = 15,
                         threshold = 60,
                         basinright = 75,
                         noise_magnitude = 4,
                         inits = c(x = 60,x = 60,x = 60),
                         duration = 20)

landscape_plotter(simulation2)
timeseries_plotter(simulation2)
#Takes about thirty seconds to make the gif for these parameter values
landscape_gifer(simulation2)

anim_save("simulation2.gif")



##low resilience system with different initial conditions
simulation3 <- simulator(basinleft = 15,
                         threshold = 60,
                         basinright = 75,
                         noise_magnitude = 4,
                         inits = c(x = 60, x = 40,x = 70),
                         duration = 20)

landscape_plotter(simulation3)
timeseries_plotter(simulation3)
#Takes about thirty seconds to make the gif for these parameter values
landscape_gifer(simulation3)

anim_save("simulation3.gif")

##low resilience system with same initial conditions
simulation4 <- simulator(basinleft = 30,
                         threshold = 70,
                         basinright = 130,
                         noise_magnitude = 6,
                         inits = c(x = 60, x = 60,x = 60),
                         duration = 20)

landscape_plotter(simulation4)
timeseries_plotter(simulation4)
#Takes about thirty seconds to make the gif for these parameter values
landscape_gifer(simulation4)

anim_save("simulation4.gif")


##low resilience system with same initial conditions
simulation4 <- simulator(basinleft = 30,
                         threshold = 70,
                         basinright = 130,
                         noise_magnitude = 6,
                         inits = c(x = 20, x = 60,x = 160),
                         duration = 20)

landscape_plotter(simulation4)
timeseries_plotter(simulation4)
#Takes about thirty seconds to make the gif for these parameter values
landscape_gifer(simulation4)

anim_save("simulation4.gif")


#PART TWO
#You can work with the results of a simulation yourself. simulation1[[1]],
#for example, will access a data frame with simulation results in columns
#of $day, $N, $Potential, and $simID. You can use this to extract results
#of a simulation and compare it to other simulations.

#I am going to set up a base, focal simulation, and you can change the threshold positions to create a useful plot
basinleft_focal=30
basinright_focal=130
noise_magnitude_focal = 2
init_focal=c(x=130)
duration_focal=20

#Do this: Populate this vector with 5 threshold values between 30 and 110,
#in increasing order
threshold_positions=c(30, 50, 70, 100, 110)

#Make a vector of NAs to store results
cv_results=rep(NA,length(threshold_positions))

#Loop through the threshold values
for(threshold_index in 1:length(threshold_positions)){
  #For each threshold value, run the simulation
  sim_temp=simulator(basinleft = basinleft_focal,
                     threshold = threshold_positions[threshold_index],
                     basinright = basinright_focal,
                     noise_magnitude = noise_magnitude_focal,
                     inits = init_focal,
                     duration = duration_focal)
  
  #Calculate the variance in N over the simulation
  #and save it in the variance_results vector
  cv_results[threshold_index]=sd(sim_temp[[1]]$N)/mean(sim_temp[[1]]$N)
}

cv_results

#Do this: Take the position of the right basin and subtract
#the vector of threshold positions to get a vector of distances
#to the threshold for the right basin. These should get smaller
#as the threshold moves to the right
rightbasin_distancetothreshold = basinright_focal - threshold_positions

#Do this: Use the "plot(nameofxvector,nameofyvector)" command
#to plot how the distance to threshold of the right basin affects
#the cv, aka, coefficient of variation. This is very closely related to variance
#as it is the standard deviation (square root of variance) divided by the mean
#Before you make the plot, discuss what kind of pattern
#you predict, e.g., increasing, decreasing, etc., and why you predict that.
plot(rightbasin_distancetothreshold, cv_results)

#Do this: Then copy and paste your plot into your google doc, discuss, and write
#a little bit about whether your prediction was right or wrong and why.

