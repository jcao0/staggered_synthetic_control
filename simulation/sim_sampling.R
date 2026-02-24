## Simulation Study
## Replication Materials

# This scripts generate simulated samples for different cases
# One case is defined by a combination of: 
#   - N (number of units) 
#   - T (number of pre-treatment periods)
#   - r (number of latent factors)

# %% set environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls(all=TRUE))
set.seed(123)

# install.packages(c("R.matlab", "doParallel", "foreach", "doRNG","MASS"))
library(R.matlab)
library(doParallel)
library(foreach)
library(doRNG)
library(MASS)
source("functions/sampling.R")

# %% pre-set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Define cases  
NN  <- c(33) # vector of target unit numbers 
TT <- c(15,42,157) # vector of target pre-treatment periods numbers
RR <- c(3,6) # vector of target latent factors numbers

# number of simulations per case
sims<-1000

## other parameters for DGP
S <- 7 # number of post-treatment periods
FE <- TRUE # include fixed effects
p <- 0 # number of observed covariates
AR1 <- 0.5 # AR(1) coefficient for latent factors and time fixed effects


# %% Generate simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create all combinations of cases
cases <- expand.grid(TT = TT, NN = NN, RR = RR) # Create all combinations of TT, N, and r
ncases <- nrow(cases)

# Generate treatment matrix for each N-S combination
list_A = list()
for(N in NN) {
    repeat { # randomly generate A matrix until at least r units are treated for feasibility of gyc, and longest treatment time is S and shortest treatment time is at least 2 for results comparability
        A = generate_treatment_matrix(N=N, S=S)
        S_treat = colSums(A) # treatment times for each unit
        S_first = max(S_treat[which(colSums(A) > 0)]) # treatment periods of first treated unit
        S_last = min(S_treat[which(colSums(A) > 0)]) # treatment periods of last treated unit
        N_treated = sum(colSums(A) > 0)
        if (N-N_treated >= min(RR) & S_first == S & S_last >= 2) break
    }
    list_A[[as.character(N)]] = A
}


# Register multiple cores
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)
cat("Cores: ", cores, "; cases: ", ncases, "\n", sep = "")
begin.time <- Sys.time()

# Start simulation
sim_data <- list()
sim_data_mat <- list()  # for MATLAB-compatible structure
for (case in 1:ncases) {
    N <- cases$NN[case]
    T <- cases$TT[case]
    r <- cases$RR[case]
    
    A = list_A[[as.character(N)]]


    te <- c(1:S)
    if( S >=7 ) te[7:length(te)] = 7 # cap long-time treatment effect time at 7


    N_treated = sum(colSums(A) > 0)
    case_name <- paste0("N", N, "_T", T, "_r", r)

    cat("Case", case, ": sims = ", sims,
        ", N = ", N, ", T = ", T, ", r = ", r, "\n", sep = "")
    cat("# of treated units: ", N_treated, "\n", sep = "")

    # simulate data in parallel
    onecase <- foreach(i = 1:sims,
                       .combine = 'c',
                       .multicombine = TRUE,
                       .inorder = FALSE,
                       .packages = c("MASS"),
                       .options.RNG = 123 # set seed
    ) %dorng% {
        # generate a random sample: 
        panel <- simulate_data_fm(N = N, T = T, S = S, p = p, r = r, AR1 = AR1,
                          D.sd = 0, te = te,
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
    # save(sim_data, file = "sim_data.RData") # store intermediate results to avoid loss by crashes 
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
