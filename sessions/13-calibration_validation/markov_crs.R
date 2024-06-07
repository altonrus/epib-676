#### CRS Markov model in a function ####
markov_crs <- function(v.params) {
  with(as.list(v.params), {
    ## Markov model parameters
    n.t  <- 60                        # time horizon, number of cycles
    v.n  <- c("NED", "Mets", "Death") # the 3 states of the model
    n.s <- length(v.n)                # number of health states 
    
    # Transition probabilities 
    # p.Mets    = 0.10         	# probability to become sicker when sick
    # p.DieMets = 0.05        	  # hazard ratio of death in sick vs healthy
    
    ####### INITIALIZATION ##########################################
    # create the cohort trace
    m.M <- matrix(NA, nrow = n.t + 1 , 
                  ncol = n.s,
                  dimnames = list(0:n.t, v.n)) # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)
    
    m.M[1, ] <- c(1, 0, 0)                     # initialize Markov trace
    
    # create transition probability matrix for NO treatment
    m.P <- matrix(0,
                  nrow = n.s, 
                  ncol = n.s,
                  dimnames = list(v.n, v.n))
    # fill in the transition probability array
    ### From NED
    m.P["NED", "NED"]   <- 1 - (p.Mets)
    m.P["NED", "Mets"]  <- p.Mets
    m.P["NED", "Death"] <- 0            # Not allowed to die from cancer in NED state
    ### From Mets
    m.P["Mets", "NED"]   <- 0
    m.P["Mets", "Mets"]  <- 1 - (p.DieMets)
    m.P["Mets", "Death"] <- p.DieMets
    ### From Death
    m.P["Death", "Death"] <- 1
    
    # check rows add up to 1
    if (!isTRUE(all.equal(as.numeric(rowSums(m.P)), as.numeric(rep(1, n.s))))) {
      stop("This is not a valid transition Matrix")
    }
    
    ############# PROCESS ###########################################
    
    for (t in 1:n.t){                   # throughout the number of cycles
      m.M[t + 1, ] <- m.M[t, ] %*% m.P  # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    
    ####### EPIDEMIOLOGICAL OUTPUT  ###########################################
    #### Overall Survival (OS) ####
    v.os <- 1 - m.M[, "Death"]  # calculate the overall survival (OS) probability
    
    ####### RETURN OUTPUT  ###########################################
    out <- list(Surv = v.os[-c(1:2)])
    
    return(out)
  }
  )
}
