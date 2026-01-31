## Simulation Study
## Replication Materials


# This scripts generate simulated samples for different cases
# One case is defined by a combination of N (number of units) and T0 (number of pre-treatment periods)

# %% set environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls(all=TRUE))

set.seed(123)
setwd("/Users/wuhang/Desktop/metric/sc/replication_files/simulation")
library(R.matlab)
library(doParallel)
library(foreach)
library(MASS)
source("functions/sampling.R")

# %% pre-set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sims<-100 # number of simulations per case

# define cases
TT0 <- c(50) # vector of target pre-treatment periods
NN  <- c(20) # vector of target units
S <- 10 # number of post-treatment periods
cases <- expand.grid(TT0 = TT0, NN = NN) # Create all combinations of TT0 and NN
ncases <- nrow(cases)

# parameters for DGP
r <- 2 # number of latent factors
FE <- TRUE # include fixed effects
p <- 0 # number of observed covariates
AR1 <- 0.5 # AR(1) coefficient for latent factors and time fixed effects



# %% Generate simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# register multiple cores
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)
cat("Cores: ", cores, "; cases: ", ncases, "\n", sep = "")
begin.time <- Sys.time()
# start simulation
sim_data <- list()
sim_data_mat <- list()  # for MATLAB-compatible structure
for (case in 1:ncases) {
    N <- cases$NN[case]
    T0 <- cases$TT0[case]
    case_name <- paste0("N", N, "_T0", T0)

    cat("Simulate: sims = ", sims,
        ", N = ", N, ", T0 = ", T0, ", S = ", S, "\n", sep = "")
    
    # generate treatment matrix
    repeat { # randomly generate A matrix until at least r+1 units are treated, for feasibility of gyc
        A = get_treatment_matrix(N=N, S=S, r=2)
        N_treated = sum(colSums(A) > 0)
        if (N-N_treated >= r+1) break
    }
    cat("Treated units proportion: ", N_treated / N, "\n", sep = "")
    
    # simulate data in parallel
    onecase <- foreach(i = 1:sims,
                       .combine = 'c',
                       .multicombine = TRUE,
                       .inorder = FALSE,
                       .packages = c("MASS")
    ) %dopar% {
        # generate a random sample: 
        panel <- simulate_data_fm(N = N, T0 = T0, S = S, p = p, r = r, AR1 = AR1,
                          D.sd = 0, te = NULL,
                          A = A,
                          mu = 5,
                          FE = FE)
        list(panel)
    }
    sim_data[[case_name]] <- onecase

    # Convert simulation result to MATLAB-style arrays
    onecase_list <- lapply(onecase, function(df) {
    required_cols <- c('id', 'time', 'Y', 'D','att') # remain only id, time, D, Y for matlab
    df_sel <- df[, required_cols]
    df_sel[] <- lapply(df_sel, as.numeric)
    as.matrix(df_sel)
    })
    rows_NT <- nrow(onecase_list[[1]])
    num_col <- ncol(onecase_list[[1]])
    num_rep <- length(onecase_list)
    onecase_array <- array(NA, dim = c(rows_NT, num_col, num_rep))
    for (i in seq_along(onecase_list)) {
      onecase_array[,,i] <- onecase_list[[i]]
    }
    sim_data_mat[[case_name]] <- onecase_array

    # Save RData intermediate
    save(sim_data, file = "sim_data.RData") # store intermediate results to avoid loss by crashes 
}
stopCluster(cl) # stop parallel computing

time<-Sys.time()-begin.time
print(time)

# %%Save simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# check structure of simulated data
if (interactive()) {
  cat("Number of cases:", length(sim_data), "\n")
  cat("Case names:", paste(names(sim_data), collapse = ", "), "\n")
  cat("Number of simulations in first case:", length(sim_data[[1]]), "\n")
  cat("First 5 rows of first simulated panel (first case, first sim):\n")
  print(print(head(sim_data[[1]][[1]], 2)))
}

# save simulated data
save(sim_data,file="data/sim_data.RData")
do.call(writeMat, c(list("data/sim_data.mat"), sim_data_mat))

print("Simulation data saved to 'sim_data.RData' and 'sim_data.mat'.")