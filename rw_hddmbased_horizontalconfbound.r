RW_hddmbased_horizontalconfbound <- function(params, samples=1000, dt=1e-4, intra_sv=.1, t2time = 1, bound2_height = .5){
  #returns simulated RTs from simulation the whole drift-process
  #required arguments: vector params containing a,v,t
  #
  # In this version, delayed confidence is queried once evidence reached a second horizontal bound (which is put at a+(a*bound2_height) ) 
  # Hopefully this explains the effect of variance on delayed confidence seen in the data.
  # Note that this quantification of confidence is critically different from work like Pleskac & Busemeyer
  #
  # November 18, KD
    
  if(exists('dt') == F)
    dt <- 1e-4
  if(exists('intra_sv') == F)
    intra_sv <- .1
  if(exists('samples') == F) 
    samples = 1
  if(exists('t2time') == F) 
    t2time = 1
  if(exists('bound2_height')== F) 
    bound2_height = .5
   
  nn = 5000
  a = params['a']
  v = params['v']
  t = params['t']
  
  if(exists('params$z')){
    z = params['z']
  } else{
    z = .5 #neutral starting point 
  }
  
  N <- samples
  
  #create delay (ter)
  start_delay = rep(1,N)*t
  
  #create starting_points
  starting_points = rep(1,N)*z*a  
  
  rts = rep(NA,N)
  resps = rep(NA,N)
  evidence1 <- rep(NA,N)
  evidence2 <- rep(NA,N)
  step_size = sqrt(dt)*intra_sv
  rts2 = rep(NA,N)
  
  for(i_sample in 1:N){
    crossed=F
    iter=0
    y_0 <- starting_points[i_sample]
    drift_rate <- v
    
    prob_up = 0.5*(1+sqrt(dt)/intra_sv*drift_rate)
  
    while(crossed == F){
      iter = iter+1
      # Generate nn steps
      position = ((runif(nn) < prob_up)*2 - 1) * step_size
      position[1] = position[1] + y_0
      position = cumsum(position)

      # Find boundary crossings
      cross_idx = which((position <= 0) | (position >= a))
      if(sum(cross_idx)>0){
        crossed <- T
      }else{
        # If not crossed, set last position as starting point for next nn steps to continue drift
        y_0 = position[nn]
      }
    }
  
    #Find the boundary interception
    y2 = position[cross_idx[1]]
    if(cross_idx[1]!=1){
      y1 = position[cross_idx[1]-1]
    }else{
      y1 <- y_0
    }
    m = (y2 - y1)/dt  # slope
    # y = m*x + b
    b = y2 - m*((iter-1)*nn+cross_idx[1])*dt # intercept
    if(y2 < 0){
      rt = ((0 - b) / m)
    }else{
      rt = ((a - b) / m)
    }
    rts[i_sample] = (rt + start_delay[i_sample]) #*sign(y2)
    resps[i_sample] <- sign(y2)
    evidence1[i_sample] <- y2
  }
    
  #Post-decisional evidence accumulation, lasting until evidence reaches a +- a*bound2_height or other bound
  for(i_sample in 1:N){
      crossed=F;iter=0
      y_0 <- evidence1[i_sample]

      while(crossed == F){
        iter = iter+1
        # Generate nn steps
        position = ((runif(nn) < prob_up)*2 - 1) * step_size
        position[1] = position[1] + y_0
        position = cumsum(position)
  
        # Find boundary crossings
        if(resps[i_sample]==1){
          #error awareness = other bound
          cross_idx = which((position >= a+(a*bound2_height)) | (position <= 0))
          #Error awareness is same distance as upper bound
          #cross_idx = which((position >= a+(a*bound2_height)) | (position <= a-(a*bound2_height)))
          #Error awareness is z*a
          #cross_idx = which((position >= a+(a*bound2_height)) | (position <= (z*a)))
        }else{
          #error awareness = other bound
          cross_idx = which((position <= 0-(a*bound2_height)) | (position >= a))
          #Error awareness is same distance as upper bound
          #cross_idx = which((position <= 0-(a*bound2_height)) | (position >= 0+(a*bound2_height)))
          #error awareness is z*a
          #cross_idx = which((position <= (z*a)) | (position >= 0+(a*bound2_height)))
        }
        if(sum(cross_idx)>0){
          crossed <- T
        }else{
          # If not crossed, set last position as starting point for next nn steps to continue drift
          y_0 = position[nn]
        }
      }
    
      #Find the boundary interception
      y2 = position[cross_idx[1]]
      if(cross_idx[1]!=1){
        y1 = position[cross_idx[1]-1]
      }else{
        y1 <- y_0
      }
      m = (y2 - y1)/dt  # slope
      # y = m*x + b
      b = y2 - m*((iter-1)*nn+cross_idx[1])*dt # intercept
      if(y2 < 0){
        rt = ((0 - b) / m)
      }else{
        rt = ((a - b) / m)
      }
      rts2[i_sample] = (rt + rts[i_sample]) #*sign(y2)
      evidence2[i_sample] <- y2  
    }
  
  return(cbind(rts, resps, evidence1, evidence2, rts2))
}
