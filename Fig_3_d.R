############### Modeling carbonate sedimentation rates during OAEs #####################

## This 2-layer box model explores the effect of changing microbial metabolisms on
## precipitation/dissolution dynamics in the ocean. The scenario presented here is one
## in which changes in ocean redox shift the net alkalinity:DIC ratio metabolic
## processes away from the O2/H2O system. This results in depressed export of alkalinity
## out of the surface ocean in carbonates, which we relate to the sedimentation rate in
## shallow-water platforms.Compare to Higgins et al., (2009) and Bergmann et al., (2013).

## This script reproduces the plot in Fig. 3d in the main text.

######################### Clear workspace and load packages ############################

## Clears the environment of previous variables

  remove(list = ls())

## Loads the relevant libraries. This code uses the program CO2sys (Lewis et al., 1998)
## as implemented in the R package seacarb (Lavigne and Gattuso, 2010). 

  library(seacarb)
  library(ggplot2)
  library(gridExtra)

############################# Set model parameters #####################################
 tik <-Sys.time()
## Set the water volume in the surface and deep oceans.

  Volume = 1.37e21/1.035                  ## mass in [kg] taken from CERC (1980) and 
                                          ## converted to [L] using a density of 1.035
                                          ## kg/L. 


  f = .029                                ## Sets the partitioning of ocean volume 
                                          ## between the surface ocean versus the deep 
                                          ## ocean. Calculated to match a 1/38 ratio
                                          ## between carbon in the surface and deep 
                                          ## oceans from Zeebe and Ridgwell (2011)
  
  V1 = f*Volume                           ## volume of the surface ocean
  
  V2 = (1-f)*Volume                       ## Volume of the deep ocean


## Set intial DIC and calcium values in the model.  

  DIC_1_ini = 2.1e-3                      ## Sets DIC in the surface ocean. Default 
                                          ## value is 2.1 mmol based on surface values 
                                          ## in Fig. 2 from  Higgins et al., (2009).
 
  DIC_2_ini = 2.35e-3                     ## Sets DIC in the deep ocean and porewaters. 
                                          ## Default value is +250 umol over the 
                                          ## surface ocean value based on Broeker and 
                                          ## Peng, (1982) and Higgins et al., (2009).

  grad_DIC = DIC_2_ini - DIC_1_ini        ## Sets the initial DIC gradient between the 
                                          ## surface ocean and the deep ocean.  
                                       

  Ca = 10.7e-3                            ## Set the calcium concentration of seawater 
                                          ## in molar, assumed to be constant in 
                                          ## the whole ocean. NOTE: not used now, but
                                          ## will be needed for changing ancient ocean
                                          ## states.

  A_i = 30e12                                ## Set the alkalinity flux in terramoles.
                                          ## Taken from Ridgwell (2005) and Munhoven 
                                          ## (2002).
  
  Rx = 4                                  ## This factor sets the "overproduction" of 
                                          ## carbonate in the surface ocean. This is 
                                          ## a multiple of the alkalinity flux A_i to 
                                          ## the surface ocean. The default factor of 
                                          ## 4 is from Broeker and Peng (1982).
  
  inv_r0 = 4                              ## sets the ratio of organic and carbonate
                                          ## fluxes out of the surface oceans, set to
                                          ## be 4:1 after Li et al., (1969).

  F_c = Rx*A_i/2                          ## Sets the initial flux of carbon 
                                          ## precipitated as carbonate out of the 
                                          ## surface layer. The factor of two converts 
                                          ## F_c from moles of alkalinity equivalent 
                                          ## to moles carbon using a 2:1 ratio of 
                                          ## alk:DIC. 
  
  F_d = (Rx-1)*A_i/2                      ## Sets the initial dissolution flux of 
                                          ## carbonate in the deep ocean. 
                                          ## The expression R_x - 1 ensures that the 
                                          ## net flux that escapes dissolution is 
                                          ## A_i at steady state. 
  
  F_o = inv_r0*F_c                        ## Sets the flux of organic carbon from the 
                                          ## surface ocean into the deep ocean.

  J = (F_d + F_o)/grad_DIC                ## Sets the exchange constant between the
                                          ## surface ocean and deep ocean. Following 
                                          ## Bergmann et al., (2013), this is an ad hoc 
                                          ## parameter that includes the net effect of 
                                          ## all advective and diffusive processes.

  C_i = F_c + F_o - J*grad_DIC            ## Sets the influx of carbon into the surface
                                          ## ocean. 

  del_resp_0 = 0                          ## Sets the initial ratio of alkalinity: DIC 
                                          ## for net respiration in the deep ocean. Note 
                                          ## The default is zero, a common approximation
                                          ## for oxygen as the dominant electron 
                                          ## acceptors. See both Bergmann et al., (2013) 
                                          ## and Higgins et al., (2009).

  ALK_1_ini = 2.3e-3                      ## Sets the initial concentration of 
                                          ## alkalinity in the surface ocean in mmolar. 
                                          ## Default value is 2.3 mmol based on surface 
                                          ## values in Fig. 2 from  Higgins et al., 
                                          ## (2009).

  grad_Alk = (2*F_d +F_o*del_resp_0)/J    ## Sets the steady state gradient in 
                                          ## alkalinity between the surface ocean and 
                                          ## the deep ocean.
                


  ALK_2_ini = ALK_1_ini + grad_Alk        ## Sets the starting alkalinity in the deep 
                                          ## ocean. 

  
########################## Check mass balance (debugging) ##############################

## These are implementations of the steady-state Eqns. 1-4. By definition, the values 
## listed as "check" should be zero if the model is working correctly.
  
  check_A_1 = A_i -2*F_c - F_o*del_resp_0 +J*grad_Alk
  check_A_2 = 2*F_d - F_o*del_resp_0 -J*grad_Alk
  check_C_1 = C_i - F_c - F_o  +J*grad_DIC
  check_C_2 = F_d + F_o  -J*grad_DIC

######################## Solve initial carbonate chemistry #############################  

## calculate initial Omega_calcite values using the 'carb' function from CO2sys. The 
## surface and deep oceans are parameterized using different ambient conditions 
## including temperature and pressure. 

  carb_1 <-carb(flag = 15, ALK_1_ini,    ## Solve for starting value in the surface 
                DIC_1_ini, S=35, T=17.5, ## ocean.
                P = 0, k1k2="x")
  
  carb_2 <-carb(flag = 15, ALK_2_ini,    ## Solve for starting value in the deep ocean 
                DIC_2_ini, S=35, T=2.5,
                P = 350, k1k2="x")
  
## Store output values
  
  omega_1_ini = carb_1$OmegaCalcite
  omega_2_ini = carb_2$OmegaCalcite

############## Set simulation time, number of timesteps, and step size #################

  t_max = 3e4                               ## Length of the simulation in years
  nt = 3e3+1                                ## Number of timesteps
  dt = t_max/(nt-1)                         ## Step size in years
  time_vector = seq(0,t_max,dt)             ## Store time values in vector
  
  del_resp_min =    0                       ## Starting del_resp value
  del_resp_max =  1.2                       ## Ending del_resp value
  del_resp_n   =   10                       ## Number of steps
  del_resp_step=  (del_resp_max - 
                     del_resp_min)/(del_resp_n-1)
  
  del_resp_all = seq(from = del_resp_min, 
                     to = del_resp_max, 
                     by = del_resp_step)
  
  del_resp_matrix = matrix(0,nt,1)
   

   
################## Pre-allocate empty vectors for storing outputs ######################

  DIC_1_all = matrix(0,nt,1)
  DIC_1_all[1] = DIC_1_ini

  DIC_2_all = matrix(0,nt,1)
  DIC_2_all[1] = DIC_2_ini

  Alk_1_all = matrix(0,nt,1)
  Alk_1_all[1] = ALK_1_ini

  Alk_2_all = matrix(0,nt,1)
  Alk_2_all[1] = ALK_2_ini

  Omega_1_all = matrix(0,nt,1)
  Omega_1_all[1] = omega_1_ini

  Omega_2_all = matrix(0,nt,1)
  Omega_2_all[1] = omega_2_ini
  
  flux_1       = matrix(Rx, nt,1)
  flux_2       = matrix(-(Rx-1), nt,1)
  
  ss_1_all     = matrix(0, del_resp_n,1)
  ss_2_all     = matrix(0, del_resp_n,1)
  
for (i in 1:del_resp_n) 
{ 
  bin = i
  
  del_resp_max = del_resp_all[bin]
  
 ############# Set onset and recovery intervals for perturbation ########################
   
  if (del_resp_max == 0){ del_resp_vector = matrix(0,nt,1)}
  
  else{
  
   onset = seq(from = del_resp_min, to = del_resp_max,
               by = (del_resp_max - del_resp_min)/(.1*nt-1))
   
   recovery = rev(onset)
   
   n1 = rep(del_resp_min, as.integer(.1*nt)) 
   n2 = onset
   n3 = rep(del_resp_max, as.integer(.8*nt +1))     
   del_resp_vector = c(n1,n2,n3) }
   
####################### Set up and execute the time loop ###############################

  for (time in 2:nt)
  {
      n = time
          
      if (n == 2){start_time <-Sys.time()}
                 
      grad_DIC = DIC_2_all[n - 1] - DIC_1_all[n - 1]
      grad_Alk = Alk_2_all[n - 1] - Alk_1_all[n - 1]
  
      carb_1 <- carb(flag = 15, Alk_1_all[n - 1], DIC_1_all[n - 1], S=35, T=17.5,P = 0,
                     k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                     warn="n", eos="eos80", long=1.e20, lat=1.e20)
  
      carb_2 <- carb(flag = 15, Alk_2_all[n - 1], DIC_2_all[n - 1], S=35, T=2.5,P = 350,
                     k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
                     warn="n", eos="eos80", long=1.e20, lat=1.e20)
  
      Omega_1_all[n] = carb_1$OmegaCalcite
      Omega_2_all[n] = carb_2$OmegaCalcite
  
      k1 = ((Omega_1_all[n]-1)/(omega_1_ini-1))^2
  
      if (carb_2$OmegaCalcite < 1){k2 = ((1- Omega_2_all[n])/(1-omega_2_ini))^2}
      else                        {k2 = -((Omega_2_all[n] -1)/(1-omega_2_ini))^2}
  
      d_DIC_1 = (1/V1)*(C_i + J*grad_DIC - F_o - k1*F_c)*dt
      d_DIC_2 = (1/V2)*(-J*grad_DIC + F_o + k2*F_d)*dt
      d_Alk_1 = (1/V1)*(A_i + J*grad_Alk -F_o*del_resp_vector[n] - 2*k1*F_c)*dt
      d_Alk_2 = (1/V2)*(-J*grad_Alk + F_o*del_resp_vector[n] +2*k2*F_d)*dt
  
      if (time/(nt -1) == .01){print('1% complete')}
      if (time/(nt -1) == .1){print('10% complete')}
      if (time/(nt -1) == .25){print('25% complete')}         
      if (time/(nt -1) == .50){print('50% complete')}          
      if (time/(nt -1) == .75){print('75% complete')}           
      if (time/(nt -1) == .90){print('90% complete')}         
             
      DIC_1_all[n] = DIC_1_all[n - 1] +  d_DIC_1
      DIC_2_all[n] = DIC_2_all[n - 1] +  d_DIC_2
  
      Alk_1_all[n] = Alk_1_all[n - 1] +  d_Alk_1
      Alk_2_all[n] = Alk_2_all[n - 1] +  d_Alk_2
      
      flux_1[n] = Rx*k1
      flux_2[n] = -(Rx-1)*k2
      
      if (n == 2){
      end_time <-Sys.time()
      estimate = round(nt*(end_time-start_time)/60, digits = 1)
      print(paste('Point', bin, 'of', del_resp_n, ': estimated runtime is', 
                  estimate, 'minutes', sep = ' '))
                  }
  }
   
   ss_1_all[bin]     = flux_1[n]
   ss_2_all[bin]     = flux_2[n]
   
}
  
  tok <- Sys.time()
  total_runtime <- round(difftime(tok, tik, units = 'mins'), digits = 1)
################################## Plot results ########################################
  
  ss_check <- abs((ss_1_all+ss_2_all-1))*100
  
  to_plot <- data.frame(del_resp_all,ss_1_all,ss_2_all)
  
  ggplot () + geom_point(data = to_plot, aes(x = del_resp_all, y = ss_1_all), 
                                             color = 'blue') +
              geom_point(data = to_plot, aes(x = del_resp_all, y = ss_2_all),
                                             color = 'red') +
              labs(x = c(expression(Delta)), 
                   y = 'normalized carbonte fluxes',
                   title = 'Steady state solutions') + theme_classic() +
                   theme(plot.title = element_text(hjust = 0.5))
  

  
  print(paste('Maximum error is', round(max(ss_check), digits = 1), 
              'percent', sep = ' '))
  
  print(paste('Total runtime was', total_runtime, 'minute(s)', sep = ' '))