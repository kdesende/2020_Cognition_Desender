RW_hddmbased <- function(params, samples=1000, dt=1e-4, intra_sv=.1, t2time = 1, evidence_bias,post_drift_modulation){
  #returns simulated RTs from simulation the whole drift-process
  #required arguments: vector params containing a,v,t
  #
  # evidence_bias implements some biased processing of evidence!
  # no=no bias, congruent = response-congruent-only, errors_only = only on error trials, corrects_only
  # March 18, KD
    
  if(exists('dt') == F)
    dt <- 1e-4
  if(exists('intra_sv') == F)
    intra_sv <- .1
  if(exists('samples') == F) 
    samples = 1
  if(exists('t2time') == F) 
    t2time = 1
  if(exists('evidence_bias') == F) 
    evidence_bias = 'no_bias'
  if(exists('post_drift_modulation') == F) 
    post_drift_modulation = 1

  nn = 5000
  a = params['a']
  v = params['v']
  t = params['t']
  
  if(exists('params$z')){
    z = params['z']
  } else{
    z = .5 #relative to .5 
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
    # delay = start_delay[i_sample]/dt
    # drift[i_sample] <- delay*starting_points(i_sample),drift
    # drifts.append(np.concatenate((np.ones(int(delay))*starting_points[i_sample], drift[:int(abs(rt)/dt)])))
    
    
    #Old (wrong?) way
    # y2 = position[cross_idx[1]]
    # rt <- ((iter-1)*nn+cross_idx[1])*dt
    # rts[i_sample] = (rt + start_delay[i_sample]) #*sign(y2)
    resps[i_sample] <- sign(y2)
    evidence1[i_sample] <- y2
    
    #Post-decisional evidence accumulation, lasting n2time/dt
    if(evidence_bias=='post_drift_modulation'){ #Post-decision drift is multiplied by X
      drift_rate <- drift_rate*post_drift_modulation
      prob_up = 0.5*(1+sqrt(dt)/intra_sv*drift_rate)
    }
    
    position = ((runif(t2time/dt) < prob_up)*2 - 1) * step_size
    if(evidence_bias=='congruent'){ #only process resp. congr. information
      position[position<0] <- 0
    }
    if(evidence_bias=='errors_only'){ #only process information on error trials
      position[resps[i_sample]==1] <- 0
    }
    if(evidence_bias=='corrects_only'){ #only process information on correct trials
      position[resps[i_sample]==0] <- 0
    }

    position[1] = position[1] + y2
    position = cumsum(position)
    evidence2[i_sample] <- position[length(position)]
    
    rt2 <- rts[i_sample]+((t2time/dt)*dt)
    rts2[i_sample] <- rt2
  }
  
  return(cbind(rts, resps, evidence1, evidence2, rts2))
}
