###############################################################################
############ Script to run simulation scenario 1 on the server ################
###############################################################################

########################## Load packages ######################################
library(mvtnorm)
library(rstan)
library(coda)

################### Read command line arguments ###############################
args <- commandArgs(trailingOnly = TRUE)
random_seed <- as.integer(args[1])
n_patients <- as.integer(args[2])
n_epochs <- as.integer(args[3])
n_samples <- as.integer(args[4])

###################### Generate parameter values ##############################

# set random seed for reproducibility
set.seed(random_seed)

# specify relevant parameters
k <- 3
event_duration_rate <- 0.1 # events are iid Exp(0.1) + 10 seconds

# vector to store event rates
event_rate_REM_vec <- rep(NA, n_patients)
event_rate_nonREM_vec <- rep(NA, n_patients)

# fixed effect vector (mu_NR, mu_NA, mu_RN, mu_RA, tau_NR, tau_NA, tau_RN, tau_RA, lambda_R, lambda_N)
fixed_vec <- rep(NA, 10)
fixed_vec[1:4] <- runif(n=4, min=-4, max=-2) # mu
fixed_vec[5:8] <- runif(n=4, min=0, max=1) # tau
fixed_vec[9:10] <- runif(n=2, min=0.01, max=0.03) # lambda

# extra awake fixed effects (mu_AN, mu_AR)
mu_AN <- runif(n=1, min=-3, max=-1)
mu_AR <- runif(n=1, min=-3, max=-1)

# theta_i covariance Sigma (factor model with k=3)
Lambda <- matrix(data=rnorm(n=10*k, mean=0, sd=sqrt(0.2)), nrow=10, ncol=k)
Omega <- diag(rep(0.2, 10))
Sigma <- Lambda %*% t(Lambda) + Omega

# theta_i matrix (gamma_NR, gamma_NA, gamma_RN, gamma_RA, alpha_NR, alpha_NA, alpha_RN, alpha_RA, phi_R, phi_N)
theta_mat <- rmvnorm(n=n_patients, sigma=Sigma)

################ Generate patient data using parameter values #################

### initialize storage ###

# trasition count patient blocks:
# (RR, RA, RN) no event
# (RR, RA, RN) event
# (NN, NA, NR) no event
# (NN, NA, NR) event
transition_counts <- matrix(nrow=n_patients*4, ncol=3)

# other count vectors
epochs_in_REM <- rep(NA, n_patients)
epochs_in_nonREM <- rep(NA, n_patients)
time_in_REM_noevent <- rep(NA, n_patients)
time_in_nonREM_noevent <- rep(NA, n_patients)
events_in_REM <- rep(NA, n_patients)
events_in_nonREM <- rep(NA, n_patients)

# Generate patients in a loop using the parameter values generated above
for (i in 1:n_patients) {
  # define patient's random effect vector
  theta_curr <- theta_mat[i,]
  
  # define vectors to record sleep stages and indicator of events
  sleep_stage_sim <- rep(NA, (n_epochs+1)) # 1 extra index to finish bookkeeping for final epoch
  event_vec <- rep(0, n_epochs*30)
  
  # compute the transition probabilities for rem and non-REM given that there is or is
  # not an event currently happening
  # for all vecs, the probabilities represent (11,12,13)
  REM_transition_prob <- rep(NA,3)
  REM_transition_prob_event <- rep(NA,3)
  nonREM_transition_prob <- rep(NA,3)
  nonREM_transition_prob_event <- rep(NA,3)
  Awake_transition_prob <- rep(NA, 3)
  
  # REM, no event
  REM_transition_prob[1] <- exp(fixed_vec[4] + theta_curr[4]) /
    (1 + exp(fixed_vec[4] + theta_curr[4]) + exp(fixed_vec[3] + theta_curr[3]))
  REM_transition_prob[2] <- 1 /
    (1 + exp(fixed_vec[4] + theta_curr[4]) + exp(fixed_vec[3] + theta_curr[3]))
  REM_transition_prob[3] <- exp(fixed_vec[3] + theta_curr[3]) /
    (1 + exp(fixed_vec[4] + theta_curr[4]) + exp(fixed_vec[3] + theta_curr[3]))
  
  # REM, event
  REM_transition_prob_event[1] <- exp(fixed_vec[4] + theta_curr[4] + fixed_vec[8] + theta_curr[8]) /
    (1 + exp(fixed_vec[4] + theta_curr[4] + fixed_vec[8] + theta_curr[8]) +
       exp(fixed_vec[3] + theta_curr[3] + fixed_vec[7] + theta_curr[7]))
  REM_transition_prob_event[2] <- 1 /
    (1 + exp(fixed_vec[4] + theta_curr[4] + fixed_vec[8] + theta_curr[8]) +
       exp(fixed_vec[3] + theta_curr[3] + fixed_vec[7] + theta_curr[7]))
  REM_transition_prob_event[3] <- exp(fixed_vec[3] + theta_curr[3] + fixed_vec[7] + theta_curr[7]) /
    (1 + exp(fixed_vec[4] + theta_curr[4] + fixed_vec[8] + theta_curr[8]) +
       exp(fixed_vec[3] + theta_curr[3] + fixed_vec[7] + theta_curr[7]))
  
  # nonREM, no event
  nonREM_transition_prob[1] <- exp(fixed_vec[2] + theta_curr[2]) /
    (1 + exp(fixed_vec[2] + theta_curr[2]) + exp(fixed_vec[1] + theta_curr[1]))
  nonREM_transition_prob[2] <- exp(fixed_vec[1] + theta_curr[1]) /
    (1 + exp(fixed_vec[2] + theta_curr[2]) + exp(fixed_vec[1] + theta_curr[1]))
  nonREM_transition_prob[3] <- 1 /
    (1 + exp(fixed_vec[2] + theta_curr[2]) + exp(fixed_vec[1] + theta_curr[1]))
  
  # nonREM, no event
  nonREM_transition_prob_event[1] <- exp(fixed_vec[2] + theta_curr[2] + fixed_vec[6] + theta_curr[6]) /
    (1 + exp(fixed_vec[2] + theta_curr[2] + fixed_vec[6] + theta_curr[6]) +
       exp(fixed_vec[1] + theta_curr[1] + fixed_vec[5] + theta_curr[5]))
  nonREM_transition_prob_event[2] <- exp(fixed_vec[1] + theta_curr[1] + fixed_vec[5] + theta_curr[5]) /
    (1 + exp(fixed_vec[2] + theta_curr[2] + fixed_vec[6] + theta_curr[6]) +
       exp(fixed_vec[1] + theta_curr[1] + fixed_vec[5] + theta_curr[5]))
  nonREM_transition_prob_event[3] <- 1 /
    (1 + exp(fixed_vec[2] + theta_curr[2] + fixed_vec[6] + theta_curr[6]) +
       exp(fixed_vec[1] + theta_curr[1] + fixed_vec[5] + theta_curr[5]))
  
  # Awake
  Awake_transition_prob[1] <- 1 / (1 + exp(mu_AR) + exp(mu_AN))
  Awake_transition_prob[2] <- exp(mu_AR) / (1 + exp(mu_AR) + exp(mu_AN))
  Awake_transition_prob[3] <- exp(mu_AN) / (1 + exp(mu_AR) + exp(mu_AN))
  
  # interevent time rate
  event_rate_REM <- fixed_vec[9] * exp(theta_curr[9])
  event_rate_nonREM <- fixed_vec[10] * exp(theta_curr[10])
  
  # always start the night with Awake
  curr_epoch <- 1
  sleep_stage_sim[1] <- 11
  
  # alternate between simulating (possibly censored) events and sleep stage transitions
  num_events_REM <- 0
  num_events_nonREM <- 0
  next_event_start <- -1
  next_event_end <- -1
  
  while (curr_epoch <= n_epochs) {
    # if current stage is Awake, do any accounting from the previous stage if needed,
    # then draw the next stage, and if relevant, the next event
    if (sleep_stage_sim[curr_epoch] == 11) {
      # If the event ended in the previous epoch:
      # First, simulate an event from the end of the previous one based on the previous sleep stage.
      # Then, evaluate whether anything further needs to be done. (If the event starts after the start
      # of the next epoch and the sleep stage changed.)
      if (curr_epoch > 1 & ceiling(next_event_end / 30) < curr_epoch) {
        if (sleep_stage_sim[curr_epoch-1] == 12) {
          next_event_start <- next_event_end + ceiling(rexp(n=1, rate=event_rate_REM))
          next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
          # check if the next event starts in the previous epoch, and increment the counter accordingly
          if (ceiling(next_event_start / 30) < curr_epoch) {
            num_events_REM <- num_events_REM + 1
            event_vec[next_event_start:(min(next_event_end, (curr_epoch-1)*30))] <- 1
          }
          # check that the new event actually crossed into the new epoch or beyond
          while (ceiling(next_event_end / 30) < curr_epoch) {
            next_event_start <- next_event_end + ceiling(rexp(n=1, rate=event_rate_REM))
            next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
            if (ceiling(next_event_start / 30) < curr_epoch) {
              num_events_REM <- num_events_REM + 1
              event_vec[next_event_start:(min(next_event_end, (curr_epoch-1)*30))] <- 1
            }
          }
        }
        else if (sleep_stage_sim[curr_epoch-1] == 13) {
          next_event_start <- next_event_end + ceiling(rexp(n=1, rate=event_rate_nonREM))
          next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
          # check if the next event starts in the previous epoch, and increment the counter accordingly
          if (ceiling(next_event_start / 30) < curr_epoch) {
            num_events_nonREM <- num_events_nonREM + 1
            event_vec[next_event_start:(min(next_event_end, (curr_epoch-1)*30))] <- 1
          }
          # check that the new event actually crossed into the new epoch or beyond
          while (ceiling(next_event_end / 30) < curr_epoch) {
            next_event_start <- next_event_end + ceiling(rexp(n=1, rate=event_rate_nonREM))
            next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
            if (ceiling(next_event_start / 30) < curr_epoch) {
              num_events_nonREM <- num_events_nonREM + 1
              event_vec[next_event_start:(min(next_event_end, (curr_epoch-1)*30))] <- 1
            }
          }
        }
      }
      
      # draw next stage
      sleep_stage_sim[curr_epoch + 1] <-sample(x=c(11,12,13), size=1, prob=Awake_transition_prob)
      curr_epoch <- curr_epoch + 1
      if (sleep_stage_sim[curr_epoch] == 12) {
        next_event_start <- (curr_epoch-1)*30 + ceiling(rexp(n=1, rate=event_rate_REM))
        next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
      }
      else if (sleep_stage_sim[curr_epoch] == 13) {
        next_event_start <- (curr_epoch-1)*30 + ceiling(rexp(n=1, rate=event_rate_nonREM))
        next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
      }
    }
    
    # Otherwise, things are more complicated
    else {
      # if event hasn't happened yet, draw the next stage and increment the epoch
      if (ceiling(next_event_start / 30) > curr_epoch) {
        if (sleep_stage_sim[curr_epoch] == 12) {
          sleep_stage_sim[curr_epoch + 1] <- sample(x=c(11,12,13), size=1, prob=REM_transition_prob)
        }
        else if (sleep_stage_sim[curr_epoch] == 13) {
          sleep_stage_sim[curr_epoch + 1] <- sample(x=c(11,12,13), size=1, prob=nonREM_transition_prob)
        }
        curr_epoch <- curr_epoch + 1
        if (sleep_stage_sim[curr_epoch] != sleep_stage_sim[curr_epoch-1]) {
          # if they switch stages, only draw a new event if there was not currently an event happening
          if (sleep_stage_sim[curr_epoch] == 12) {
            next_event_start <- (curr_epoch-1)*30 + ceiling(rexp(n=1, rate=event_rate_REM))
            next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
          }
          else if (sleep_stage_sim[curr_epoch] == 13) {
            next_event_start <- (curr_epoch-1)*30 + ceiling(rexp(n=1, rate=event_rate_nonREM))
            next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
          }
        }
      }
      
      # if event is occurring in the current epoch, account for this in the transition
      else if (ceiling(next_event_start / 30) <= curr_epoch) { 
        # if the event starts here, increment the appropriate counter
        if (ceiling(next_event_start / 30) == curr_epoch) {
          event_vec[next_event_start:(min(next_event_end, (curr_epoch)*30))] <- 1
          if (sleep_stage_sim[curr_epoch] == 12) {
            num_events_REM <- num_events_REM + 1
          }
          else if (sleep_stage_sim[curr_epoch] == 13) {
            num_events_nonREM <- num_events_nonREM + 1
          }
        }
        # otherwise, still update event_vec
        else {
          event_vec[((curr_epoch-1) * 30 + 1):(min(next_event_end, (curr_epoch)*30))] <- 1
        }
        # regardless, choose the next stage based on an event currently happening
        if (sleep_stage_sim[curr_epoch] == 12) {
          sleep_stage_sim[curr_epoch + 1] <-
            sample(x=c(11,12,13), size=1, prob=REM_transition_prob_event)
        }
        else if (sleep_stage_sim[curr_epoch] == 13) {
          sleep_stage_sim[curr_epoch + 1] <-
            sample(x=c(11,12,13), size=1, prob=nonREM_transition_prob_event)
        }
        curr_epoch <- curr_epoch + 1
      }
      
      # If the event ended in the previous epoch:
      # First, simulate an event from the end of the previous one based on the previous sleep stage.
      # Then, evaluate whether anything further needs to be done. (If the event starts after the start
      # of the next epoch and the sleep stage changed.)
      if (ceiling(next_event_end / 30) < curr_epoch) {
        if (sleep_stage_sim[curr_epoch-1] == 12) {
          next_event_start <- next_event_end + ceiling(rexp(n=1, rate=event_rate_REM))
          next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
          # check if the next event starts in the previous epoch, and increment the counter accordingly
          if (ceiling(next_event_start / 30) < curr_epoch) {
            num_events_REM <- num_events_REM + 1
            event_vec[next_event_start:(min(next_event_end, (curr_epoch-1)*30))] <- 1
          }
          # check that the new event actually crossed into the new epoch or beyond
          while (ceiling(next_event_end / 30) < curr_epoch) {
            next_event_start <- next_event_end + ceiling(rexp(n=1, rate=event_rate_REM))
            next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
            if (ceiling(next_event_start / 30) < curr_epoch) {
              num_events_REM <- num_events_REM + 1
              event_vec[next_event_start:(min(next_event_end, (curr_epoch-1)*30))] <- 1
            }
          }
        }
        else if (sleep_stage_sim[curr_epoch-1] == 13) {
          next_event_start <- next_event_end + ceiling(rexp(n=1, rate=event_rate_nonREM))
          next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
          # check if the next event starts in the previous epoch, and increment the counter accordingly
          if (ceiling(next_event_start / 30) < curr_epoch) {
            num_events_nonREM <- num_events_nonREM + 1
            event_vec[next_event_start:(min(next_event_end, (curr_epoch-1)*30))] <- 1
          }
          # check that the new event actually crossed into the new epoch or beyond
          while (ceiling(next_event_end / 30) < curr_epoch) {
            next_event_start <- next_event_end + ceiling(rexp(n=1, rate=event_rate_nonREM))
            next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
            if (ceiling(next_event_start / 30) < curr_epoch) {
              num_events_nonREM <- num_events_nonREM + 1
              event_vec[next_event_start:(min(next_event_end, (curr_epoch-1)*30))] <- 1
            }
          }
        }
        if (ceiling(next_event_start / 30) >= curr_epoch &
            sleep_stage_sim[curr_epoch] != sleep_stage_sim[curr_epoch-1]) {
          if (sleep_stage_sim[curr_epoch] == 12) {
            next_event_start <- (curr_epoch-1)*30 + ceiling(rexp(n=1, rate=event_rate_REM))
            next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
          }
          else if (sleep_stage_sim[curr_epoch] == 13) {
            next_event_start <- (curr_epoch-1)*30 + ceiling(rexp(n=1, rate=event_rate_nonREM))
            next_event_end <- ceiling(next_event_start + rexp(n=1, rate=event_duration_rate)) + 10
          }
        }
      }
    }
  }
  
  # process simulated data to get needed counts
  curr_stage_vec <- sleep_stage_sim[1:(n_epochs-1)]
  next_stage_vec <- sleep_stage_sim[2:n_epochs]
  stage_seconds_vec <- rep(sleep_stage_sim[1:n_epochs], each=30)
  event_epoch_vec <- rep(0, (n_epochs-1))
  for (j in 1:(n_epochs-1)) {
    if (sum(event_vec[((j-1)*30 + 1):(j*30)]) > 0) {
      event_epoch_vec[j] <- 1
    }
  }
  
  
  # transition counts:
  # RR no event
  transition_counts[(i-1)*4 + 1, 1] <- sum(curr_stage_vec == 12 &
                                             next_stage_vec == 12 &
                                             event_epoch_vec == 0)
  # RA no event
  transition_counts[(i-1)*4 + 1, 2] <- sum(curr_stage_vec == 12 &
                                             next_stage_vec == 11 &
                                             event_epoch_vec == 0)
  # RN no event
  transition_counts[(i-1)*4 + 1, 3] <- sum(curr_stage_vec == 12 &
                                             next_stage_vec == 13 &
                                             event_epoch_vec == 0)
  # RR event
  transition_counts[(i-1)*4 + 2, 1] <- sum(curr_stage_vec == 12 &
                                             next_stage_vec == 12 &
                                             event_epoch_vec == 1)
  # RA event
  transition_counts[(i-1)*4 + 2, 2] <- sum(curr_stage_vec == 12 &
                                             next_stage_vec == 11 &
                                             event_epoch_vec == 1)
  # RN event
  transition_counts[(i-1)*4 + 2, 3] <- sum(curr_stage_vec == 12 &
                                             next_stage_vec == 13 &
                                             event_epoch_vec == 1)
  # NN no event
  transition_counts[(i-1)*4 + 3, 1] <- sum(curr_stage_vec == 13 &
                                             next_stage_vec == 13 &
                                             event_epoch_vec == 0)
  # NA no event
  transition_counts[(i-1)*4 + 3, 2] <- sum(curr_stage_vec == 13 &
                                             next_stage_vec == 11 &
                                             event_epoch_vec == 0)
  # NR no event
  transition_counts[(i-1)*4 + 3, 3] <- sum(curr_stage_vec == 13 &
                                             next_stage_vec == 12 &
                                             event_epoch_vec == 0)
  # NN event
  transition_counts[(i-1)*4 + 4, 1] <- sum(curr_stage_vec == 13 &
                                             next_stage_vec == 13 &
                                             event_epoch_vec == 1)
  # NA event
  transition_counts[(i-1)*4 + 4, 2] <- sum(curr_stage_vec == 13 &
                                             next_stage_vec == 11 &
                                             event_epoch_vec == 1)
  # NR event
  transition_counts[(i-1)*4 + 4, 3] <- sum(curr_stage_vec == 13 &
                                             next_stage_vec == 12 &
                                             event_epoch_vec == 1)
  
  
  # other count vectors
  epochs_in_REM[i] <- sum(curr_stage_vec == 12)
  epochs_in_nonREM[i] <- sum(curr_stage_vec == 13)
  events_in_REM[i] <- num_events_REM
  events_in_nonREM[i] <- num_events_nonREM
  time_in_REM_noevent[i] <- sum(stage_seconds_vec == 12 &
                                  event_vec == 0) + num_events_REM
  time_in_nonREM_noevent[i] <- sum(stage_seconds_vec == 13 &
                                  event_vec == 0) + num_events_nonREM
  event_rate_REM_vec[i] <- event_rate_REM
  event_rate_nonREM_vec[i] <- event_rate_nonREM
}

# save count vectors
#scen1_num_events_REM <- readRDS("scen1_num_events_REM.rds")
#scen1_num_events_nonREM <- readRDS("scen1_num_events_nonREM.rds")
#scen1_num_epochs_REM <- readRDS("scen1_num_epochs_REM.rds")
#scen1_num_epochs_nonREM <- readRDS("scen1_num_epochs_nonREM.rds")
#scen1_num_events_REM[,random_seed] <- events_in_REM
#scen1_num_events_nonREM[,random_seed] <- events_in_nonREM
#scen1_num_epochs_REM[,random_seed] <- epochs_in_REM
#scen1_num_epochs_nonREM[,random_seed] <- epochs_in_nonREM
#saveRDS(scen1_num_events_REM, "scen1_num_events_REM.rds")
#saveRDS(scen1_num_events_nonREM, "scen1_num_events_nonREM.rds")
#saveRDS(scen1_num_epochs_REM, "scen1_num_epochs_REM.rds")
#saveRDS(scen1_num_epochs_nonREM, "scen1_num_epochs_nonREM.rds")


# make X and Z matrix to use in stan model
XZ_single <- t(matrix(c(1,0,0,0,
                        1,0,1,0,
                        0,1,0,0,
                        0,1,0,1), ncol=4))

X <- kronecker(matrix(rep(1,n_patients), nrow=n_patients), XZ_single)
Z <- kronecker(diag(n_patients), XZ_single)


############################## Fit stan model #################################

# load model
joint_model_factor <-
  stan_model("joint_fit_factormodel.stan")

# fit model
options(mc.cores = parallel::detectCores())
joint_fit_factor <- sampling(joint_model_factor,
                             list(N=n_patients,
                                  k=k, # k=3 factors as in data generation
                                  Y=transition_counts,
                                  X=X,
                                  Z=Z,
                                  Y_seconds=matrix(data=c(time_in_REM_noevent,
                                                          time_in_nonREM_noevent),
                                                   ncol=2),
                                  Y_eventcounts=matrix(data=c(events_in_REM,
                                                               events_in_nonREM),
                                                        ncol=2)),
                             iter=n_samples)


########################## Summarize accuracy and output ######################


# extract posterior samples
post_samples_event <- extract(joint_fit_factor)

# ESS and Rhat -- fixed effects

min_ess_fixed <- 10000
max_rhat_fixed <- 0

# lambda_N
print("lambda_N ESS and rhat:")
ess_curr <- effectiveSize(mcmc(data=matrix(post_samples_event$lambda_N, ncol=1)))
rhat_curr <- Rhat(matrix(post_samples_event$lambda_N, ncol=4))
print(ess_curr)
print(rhat_curr)
if (ess_curr < min_ess_fixed) {
  min_ess_fixed <- ess_curr
}
if (rhat_curr > max_rhat_fixed) {
  max_rhat_fixed <- rhat_curr
}


# lambda_R
print("lambda_R ESS and rhat:")
ess_curr <- effectiveSize(mcmc(data=matrix(post_samples_event$lambda_R, ncol=1)))
rhat_curr <- Rhat(matrix(post_samples_event$lambda_R, ncol=4))
print(ess_curr)
print(rhat_curr)
if (ess_curr < min_ess_fixed) {
  min_ess_fixed <- ess_curr
}
if (rhat_curr > max_rhat_fixed) {
  max_rhat_fixed <- rhat_curr
}

# mus and taus
print("gamma and alpha ESS and rhat values:")
for (i in 1:4) {
  for (j in 1:2) {
    ess_curr <- effectiveSize(mcmc(data=matrix(post_samples_event$mu_tau[,i,j], ncol=1)))
    rhat_curr <- Rhat(matrix(post_samples_event$mu_tau[,i,j], ncol=4))
    print(ess_curr)
    print(rhat_curr)
    if (ess_curr < min_ess_fixed) {
      min_ess_fixed <- ess_curr
    }
    if (rhat_curr > max_rhat_fixed) {
      max_rhat_fixed <- rhat_curr
    }
  }
}

# ESS and Rhat -- random effects

# minimum ess and maximum rhat for random effects
min_ess_random <- 10000
max_rhat_random <- 0

# phi_N and phi_R
for (i in 1:n_patients) {
  ess_phi_N <- effectiveSize(mcmc(matrix(data=post_samples_event$phi_N[,i], ncol=1)))
  ess_phi_R <- effectiveSize(mcmc(matrix(data=post_samples_event$phi_R[,i], ncol=1)))
  rhat_phi_N <- Rhat(matrix(data=post_samples_event$phi_N[,i], ncol=4))
  rhat_phi_R <- Rhat(matrix(data=post_samples_event$phi_R[,i], ncol=4))
  
  if (ess_phi_N < min_ess_random) {
    min_ess_random <- ess_phi_N
  }
  if (ess_phi_R < min_ess_random) {
    min_ess_random <- ess_phi_R
  }
  if (rhat_phi_N > max_rhat_random) {
    max_rhat_random <- rhat_phi_N
  }
  if (rhat_phi_R > max_rhat_random) {
    max_rhat_random <- rhat_phi_R
  }
}

# alphas and gammas
for (i in 1:(n_patients*4)) {
  for (j in 1:2) {
    curr_ess <- effectiveSize(mcmc(matrix(data=post_samples_event$gamma_alpha[,i,j], ncol=1)))
    curr_rhat <- Rhat(matrix(data=post_samples_event$gamma_alpha[,i,j], ncol=4))
    if (curr_ess < min_ess_random) {
      min_ess_random <- curr_ess
    }
    if (curr_rhat > max_rhat_random) {
      max_rhat_random <- curr_rhat
    }
  }
}

# print out the min ess and max rhat for the random effects
print("Min ESS for random effects:")
print(min_ess_random)
print("Max Rhat for random effects:")
print(max_rhat_random)


#### Cov(theta_i) MSE and coverage

# theta_i matrix (gamma_NR, gamma_NA, gamma_RN, gamma_RA, alpha_NR, alpha_NA, alpha_RN, alpha_RA, phi_R, phi_N)
# stan version (gamma_RA, gamma_NA, alpha_RA, alpha_NA, gamma_RN, gamma_NR, alpha_RN, alpha_NR, phi_R, phi_N)
# matching up indices for theta in stan model vs. here
stan_indices <- c(4,2,8,6,3,1,7,5,9,10) # use when calling the Sigma generated here (not the stan version)

# form array of covariance samples
# n_samples_total is based on four chains
n_samples_total <- (n_samples / 2) * 4
cov_mat_samples <- array(dim=c(10,10,n_samples_total))
for (i in 1:n_samples_total) {
  cov_mat_samples[,,i] <- post_samples_event$Lambda[i,,] %*% t(post_samples_event$Lambda[i,,]) + diag(post_samples_event$Sigma[i,])
}
# compute MSE
cov_mse <- 0
for (i in 1:10) {
  for (j in 1:i) {
    cov_mse <- cov_mse + (mean(cov_mat_samples[i,j,]) - Sigma[stan_indices[i],stan_indices[j]])^2
  }
}
cov_mse <- cov_mse / 55 # (dividing by the number of lower-triangular elements)
# save
#scen1_cov_mat_mse <- readRDS("scen1_cov_mat_mse.rds")
#scen1_cov_mat_mse[random_seed] <- cov_mse
#saveRDS(scen1_cov_mat_mse, "scen1_cov_mat_mse.rds")


### Cov(theta_i) coverage
cov_coverage <- 0
for (i in 1:10) {
  for (j in 1:i) {
    if (quantile(cov_mat_samples[i,j,], probs=c(0.025)) < Sigma[stan_indices[i], stan_indices[j]] &
        quantile(cov_mat_samples[i,j,], probs=c(0.975)) > Sigma[stan_indices[i], stan_indices[j]]) {
      cov_coverage <- cov_coverage + 1
    }
  }
}
cov_coverage <- cov_coverage / 55 # (dividing by the number of lower-triangular elements)


# create values to save with output
cov_mean <- apply(cov_mat_samples, c(1,2), mean)
cov_lower <- apply(cov_mat_samples, c(1,2), quantile, probs=c(0.025))
cov_upper <- apply(cov_mat_samples, c(1,2), quantile, probs=c(0.975))
cov_true <- Sigma[stan_indices,stan_indices]

### estimated theta_i (posterior means)
theta_samples <- array(dim=c(10,n_patients,n_samples_total))
for (l in 1:n_samples_total) {
  for (i in 1:n_patients) {
    theta_samples[1:8,i,l] <- as.vector(post_samples_event$gamma_alpha[l,((i-1)*4 + 1):(i*4),1:2])
    theta_samples[9,i,l] <- post_samples_event$phi_R[l,i]
    theta_samples[10,i,l] <- post_samples_event$phi_N[l,i]
  }
}
theta_postmean <- apply(theta_samples, c(2,1), mean)
theta_lower <- apply(theta_samples, c(2,1), quantile, probs=0.025)
theta_upper <- apply(theta_samples, c(2,1), quantile, probs=0.975)
# save
#scen1_theta_postmean <- readRDS("scen1_theta_postmean.rds")
#scen1_theta_lower <- readRDS("scen1_theta_lower.rds")
#scen1_theta_upper <- readRDS("scen1_theta_upper.rds")
#scen1_theta_postmean[,,random_seed] <- theta_postmean
#scen1_theta_lower[,,random_seed] <- theta_lower
#scen1_theta_upper[,,random_seed] <- theta_upper
#saveRDS(scen1_theta_postmean, "scen1_theta_postmean.rds")
#saveRDS(scen1_theta_lower, "scen1_theta_lower.rds")
#saveRDS(scen1_theta_upper, "scen1_theta_upper.rds")

### true theta_i
theta_true <- theta_mat[,stan_indices]
# save
#scen1_theta_true <- readRDS("scen1_theta_true.rds")
#scen1_theta_true[,,random_seed] <- theta_true
#saveRDS(scen1_theta_true, "scen1_theta_true.rds")

### fixed effect MSE
# fixed effect vector (mu_NR, mu_NA, mu_RN, mu_RA, tau_NR, tau_NA, tau_RN, tau_RA, lambda_R, lambda_N)
# stan version -- (mu_RA, mu_NA, tau_RA, tau_NA, mu_RN, mu_NR, tau_RN, tau_NR, lambda_R, lambda_N)
# ---> same indices as above
fixed_effect_means <- c(as.vector(apply(post_samples_event$mu_tau, c(2,3), mean)),
                        mean(post_samples_event$lambda_R),
                        mean(post_samples_event$lambda_N))
true_fixed_effects <- fixed_vec[stan_indices]
fixed_mse <- mean((fixed_effect_means - true_fixed_effects)^2)
# save
#scen1_fixed_effect_mse <- readRDS("scen1_fixed_effect_mse.rds")
#scen1_fixed_effect_mse[random_seed] <- fixed_mse
#saveRDS(scen1_fixed_effect_mse, "scen1_fixed_effect_mse.rds")

### fixed effect coverage
fixed_effect_lower <- c(as.vector(apply(post_samples_event$mu_tau, c(2,3), quantile, probs=0.025)),
                        quantile(post_samples_event$lambda_R, probs=0.025),
                        quantile(post_samples_event$lambda_N, probs=0.025))
fixed_effect_upper <- c(as.vector(apply(post_samples_event$mu_tau, c(2,3), quantile, probs=0.975)),
                        quantile(post_samples_event$lambda_R, probs=0.975),
                        quantile(post_samples_event$lambda_N, probs=0.975))
fixed_coverage <- sum(fixed_effect_upper > true_fixed_effects & fixed_effect_lower < true_fixed_effects) / 10
# save
#scen1_fixed_effect_coverage <- readRDS("scen1_fixed_effect_coverage.rds")
#scen1_fixed_effect_coverage[random_seed] <- fixed_coverage
#saveRDS(scen1_fixed_effect_coverage, "scen1_fixed_effect_coverage.rds")



#################### organize quantities to save ########################
sim_results <- list(events_in_REM,
                    events_in_nonREM,
                    epochs_in_REM,
                    epochs_in_nonREM,
                    cov_mse,
                    cov_coverage,
                    cov_mean,
                    cov_lower,
                    cov_upper,
                    cov_true,
                    theta_postmean,
                    theta_lower,
                    theta_upper,
                    theta_true,
                    fixed_mse,
                    fixed_coverage,
                    fixed_effect_lower,
                    fixed_effect_upper,
                    true_fixed_effects,
                    fixed_effect_means,
                    min_ess_fixed,
                    max_rhat_fixed,
                    min_ess_random,
                    max_rhat_random)

names(sim_results) <- c("events_in_REM",
                        "events_in_nonREM",
                        "epochs_in_REM",
                        "epochs_in_nonREM",
                        "cov_mse",
                        "cov_coverage",
                        "cov_mean",
                        "cov_lower",
                        "cov_upper",
                        "cov_true",
                        "theta_postmean",
                        "theta_lower",
                        "theta_upper",
                        "theta_true",
                        "fixed_mse",
                        "fixed_coverage",
                        "fixed_effect_lower",
                        "fixed_effect_upper",
                        "true_fixed_effects",
                        "fixed_effect_means",
                        "min_ess_fixed",
                        "max_rhat_fixed",
                        "min_ess_random",
                        "max_rhat_random")


############ save results ##################
results_filename <- paste0("scen1_results_seed",
                           random_seed,
                           ".rds")

saveRDS(sim_results, results_filename)





