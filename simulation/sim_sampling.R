## Simulation Study
## Replication Materials


# This scripts generate simulated samples


setwd("/Users/wuhang/Desktop/metric/sc/script/simulation")
library(R.matlab)
library(doParallel)
library(foreach)
library(MASS)
source("functions/sampling.R")

# %% pre-set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sims<-100

TT0 <- c(50)
NN  <- c(20)
S <- 10
# Create all combinations of TT0 and NN
cases <- expand.grid(TT0 = TT0, NN = NN)
ncases <- nrow(cases)

r <- 2
p <- 0



# %% Generate simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(123)
# register multiple cores
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)
cat("Cores: ", cores, "; cases: ", ncases, "\n", sep = "")
begin.time <- Sys.time()
# start simulation
results_list <- list()
results_list_mat <- list()  # for MATLAB-compatible structure
for (case in 1:ncases) {
    N <- cases$NN[case]
    T0 <- cases$TT0[case]
    case_name <- paste0("N", N, "_T0", T0)

    cat("Simulate: sims = ", sims,
        ", N = ", N, ", T0 = ", T0, ", S = ", S, "\n", sep = "")
    
    # generate treatment matrix
    repeat { # randomly generate A matrix until at least some units are treated
        A = get_treatment_matrix(N=N, S=S, r=2)
        r_treated = mean(colSums(A) > 0)
        if (r_treated > 0) break
    }
    cat("Treated units proportion: ", r_treated, "\n", sep = "")
    
    # simulate data in parallel
    onecase <- foreach(i = 1:sims,
                       .combine = 'c',
                       .multicombine = TRUE,
                       .inorder = FALSE,
                       .packages = c("MASS")
    ) %dopar% {
        # generate a random sample: 
        panel <- simulate_data_fm(N = N, T0 = T0, S = S, p = p, r = r,
                          D.sd = 0, te = NULL,
                          A = A,
                          mu = 5,
                          FE = TRUE)
        list(panel)
    }
    results_list[[case_name]] <- onecase

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
    results_list_mat[[case_name]] <- onecase_array

    # Save RData intermediate
    # save(results_list, file = "sim_data.RData") # store intermediate results to avoid loss by crashes 
}
stopCluster(cl) # stop parallel computing

time<-Sys.time()-begin.time
print(time)

# %%Save simulated data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# check structure of saved data
if (interactive()) {
  cat("Number of cases:", length(results_list), "\n")
  cat("Case names:", paste(names(results_list), collapse = ", "), "\n")
  cat("Number of simulations in first case:", length(results_list[[1]]), "\n")
  cat("First 5 rows of first simulated panel (first case, first sim):\n")
  print(print(head(results_list[[1]][[1]], 2)))
}
a = results_list[[1]] # list of sims
a = results_list[[1]][[1]] # dataframe of one sim

# save simulated data
save(results_list,file="data/sim_data.RData")
do.call(writeMat, c(list("data/sim_data.mat"), results_list_mat))

print("Simulation data saved to 'sim_data.RData' and 'sim_data.mat'.")

#  %%%%%%% archive%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Convert results_list to MATLAB-style arrays
load("sim_data.RData")
results_list_mat <- list()  # for MATLAB-compatible structure
# convert
results_list_mat <- foreach(case_name = names(results_list), 
                            .combine = 'list',
                            .multicombine = TRUE,
                            .packages = c(),
                            .export = c()) %do% {  # use %dopar% if you want parallel
  cat("Processing case: ", case_name, "\n", sep = "")
  onecase <- results_list[[case_name]]
  # Convert each simulation to a numeric matrix
    onecase_list <- lapply(onecase, function(df) {
    required_cols <- c('id', 'time', 'Y', 'D','att') # remain only id, time, D, Y for matlab
    df_sel <- df[, required_cols]
    df_sel[] <- lapply(df_sel, as.numeric)
    as.matrix(df_sel)
  })

  # Build 3D array
  rows_NT <- nrow(onecase_list[[1]])
  num_col <- ncol(onecase_list[[1]])
  num_rep <- length(onecase_list)
  onecase_array <- array(NA, dim = c(rows_NT, num_col, num_rep))
  for (i in seq_along(onecase_list)) {
    onecase_array[,,i] <- onecase_list[[i]]
  }

  onecase_array
}
names(results_list_mat) <- names(results_list) # Set the names of the list to match your cases
a = results_list_mat[[1]] # array of one case: rows x cols x sims

# A test for parallel simulation
sims <- 500
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)
onecase <- foreach(i = 1:sims,
                   .combine = 'c',          # combine lists
                   .multicombine = TRUE,    # allow multicombine
                   .packages = c("MASS")) %dopar% {
  # simulate a dataframe
  df <- data.frame(x = rnorm(5), y = runif(5))
  
  list(df)  # wrap in list so c() combines correctly
}
a = onecase[[1]] # first dataframe
