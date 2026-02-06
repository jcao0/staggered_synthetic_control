## Simulation Study
## Replication Materials


## This scripts produce simulation results of alternative methods for event-time ATT
## Alternative methods include:
## 1. Generalized Synthetic Control (gsc) by Xu (2017)
## 2. Augmented Synthetic Control (asy) by Ben-Michael, Feller, and Rothstein (2021)
## 3. Difference-in-Differences (did) method by Callaway and Sant'Anna (2021)


# %% set environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls(all=TRUE))
set.seed(123)
setwd("/Users/wuhang/Desktop/metric/sc/submission/replication_files/simulation")


library(dplyr)
library(doParallel)
library(parallel)
library(foreach)
library(gsynth)
library(augsynth)
library(did)
source("functions/padding_matrix.R")

# load simulated data
load("data/sim_data.RData")

# %% pre-set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# methods to compare:
gsc = TRUE  # Generalized Synthetic Control
asy = TRUE  # Augmented Synthetic Control
did = FALSE  # Callaway and Sant'Anna (2021) DID method

# set small number for testing cases if needed
n_test_case = NULL # set NULL to run all cases

# set smallnumber of testing simulations if needed
n_test_sim <- NULL # set NULL to run all simulations per case





# %%Start estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# results would be stored in a nested list: method -> case -> att_e array + time array
# dim of att_e array =  event-time; stats; sims, stats = (att_hat, bias, se, CI_low, CI_up, CI_coverage)
# dim of time array = nsims; 1

# create a list for results storage
results <- list()
if(gsc==TRUE){
  results[['gsc']] <- list()
}
if(asy==TRUE){
  results[['asy']] <- list()
}
if(did==TRUE){
  results[['did']] <- list()
}

# extract simulated case names from data
WW  = names(sim_data) # names of cases
if (!is.null(n_test_case)) {
  WW <- WW[n_test_case]
}
ncases<-length(WW)
nsims <- if(!is.null(n_test_sim)) n_test_sim else length(sim_data[[1]]) # number of simulations per case


# start estimation
# register multiple cores
#cores<-detectCores()
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)
begin.time<-Sys.time()
cat("\nCores: ",cores,"; # of cases: ",ncases,"; # of sims: ",nsims,sep="")

# loop
for (case in 1:ncases) {
    case_name <- WW[case]

    # extract case data
    data_case <- sim_data[[case_name]]

    # extract case parameters
    T <- as.numeric(sub(".*T([0-9]+).*", "\\1", case_name))
    N <- as.numeric(sub("N([0-9]+)_.*", "\\1", case_name))
    r <- as.numeric(sub(".*_r([0-9]+)", "\\1", case_name))
    S <- data_case[[1]]$time %>% unique() %>% length() - T 
    
    ## annouce case
    cat("\nCase: ", case_name, "\n",sep="")
    

    # start simulations for one case
    onecase<-foreach (i=1:nsims,.combine='list', .multicombine=TRUE,.inorder=FALSE,
                      .packages=c("gsynth", "augsynth", "did","dplyr"), .verbose = FALSE) %dopar% { # parallel computing
        
        ## annouce simulations
        if (i %% 200 == 0) {
          cat(sprintf("  [Case %d/%d] Simulation %d/%d completed...\n", case, ncases, i, nsims))
          flush.console()
        }

        ## Prepare data:
        panel <- data_case[[i]]
        required_cols <- c('id', 'time', 'Y', 'D','att')
        panel <- panel[, required_cols]


        ## GSC
        if(gsc==TRUE){
        t0_gsc <- Sys.time()
        if (exists("gsc_SE")) {
          out_gsc <- gsynth(Y ~ D, data = panel, EM = TRUE, 
                        index=c("id","time"), inference="parametric",
                        se = TRUE, nboots = 1000, r = 2, CV = FALSE, # the number of factors are set to be 2
                        force = "two-way", parallel = TRUE)
          t1_gsc <- Sys.time()
          time_gsc <- as.numeric(difftime(t1_gsc, t0_gsc, units = "secs")) # used time for 1 simulation
          
          att_e.gsc = out_gsc$est.att # event time att
          att_e.gsc = att_e.gsc[as.numeric(rownames(att_e.gsc))>0,c('ATT') ]
        } else {
            out_gsc <- gsynth(Y ~ D, data = panel, 
                          index=c("id","time"),
                          se = FALSE, r = 2, CV = FALSE, # the number of factors are set to be 2
                          force = "two-way", parallel = TRUE)
                          att_e.gsc = out_gsc$eff # event time att
            t1_gsc <- Sys.time()
            time_gsc <- as.numeric(difftime(t1_gsc, t0_gsc, units = "secs")) # used time for 1 simulation

            att_e.gsc = out_gsc$att # event time att
            att_e.gsc = tail(att_e.gsc, S) # keep last S row
          }
        } else {
          att_e.gsc <- NULL
          time_gsc <- NULL
        }

        ## ASY: do not calculate ATT, as they only calulate att_e for the shortest treated periods across all treated units.
        if(asy==TRUE){
        t0_asy <- Sys.time()   
        
        ppool_syn <- multisynth(Y ~ D, unit = id, time = time,
                        panel)            
        t1_asy <- Sys.time()
        time_asy <- as.numeric(difftime(t1_asy, t0_asy, units = "secs"))
        
        att_e.asy = summary(ppool_syn)
        att_e.asy <- as.data.frame(att_e.asy$att)
        att_e.asy = att_e.asy[att_e.asy$Level =='Average' & is.na(att_e.asy$Time) ==FALSE & att_e.asy$Time>=0, ] # only post-treatment periods for shortest treated periods
        cols <- c('Estimate')
        att_e.asy <- att_e.asy[, cols]
        att_e.asy[] <- lapply(att_e.asy, as.numeric)
        att_e.asy <- unlist(att_e.asy, use.names = FALSE)

        } else {
          att_e.asy <- NULL
          time_asy <- NULL
        }
        ## DID
        if(did==TRUE){
        panel <- panel  |>  group_by(id)  |>  mutate(first.treat = ifelse(any(D == 1), min(time[D == 1]), 0))  |>  ungroup() # prepare data for DID
        
        t0_did <- Sys.time()
        out_did <- att_gt(
          yname = "Y",
          gname = "first.treat",  # <-- This should be the group/treatment timing variable
          idname = "id",
          tname = "time",
          xformla = ~1,
          data = panel,
          est_method = "data",
          control_group = "notyettreated"
        )
        did_summary_att_e <- aggte(out_did, type = "dynamic")
        t1_did <- Sys.time()
        time_did <- as.numeric(difftime(t1_did, t0_did, units = "secs"))

        att_e.did <- data.frame(
          att = did_summary_att_e$att.egt,
          time = did_summary_att_e$egt
        )
        att_e.did <- att_e.did[att_e.did$time >= 0, ]  # extract only post-treatment periods ATT
        att_e.did <- att_e.did[, c("att")]
        
        } else {
          att_e.did <- NULL
          time_did <- NULL
        }

        out = list(
          gsc_att_e = att_e.gsc,     # matrix
          gsc_time  = time_gsc,

          asy_att_e = att_e.asy,     # matrix
          asy_time  = time_asy,

          did_att_e = att_e.did,     # matrix
          did_time  = time_did
        )
        return(out)      
  }
  # collect results and store in arrays : dim =  event-time; sims
  event_names  <- c(1:S)
  sim_names   <- c(1:nsims)

  if (gsc==TRUE){
    # extract results and converts to arrays
    gsc_att_e_list <- lapply(onecase, `[[`, "gsc_att_e")
    gsc_att_e_pad <- lapply(gsc_att_e_list, pad_vector, nrow = S)
    gsc_att_e_array <- t(simplify2array(gsc_att_e_pad))
    dimnames(gsc_att_e_array) <- list(sim = sim_names,
                                      event = event_names
                                      )

    gsc_time_list  <- lapply(onecase, `[[`, "gsc_time")
    gsc_time_array <- do.call(rbind, gsc_time_list)
    rownames(gsc_time_array) <- sim_names

    # Store results
    results[['gsc']][[case_name]] <- list(
    att_e = gsc_att_e_array,
    time  = gsc_time_array
    )
  } 

  if (asy==TRUE){
    # extract results and converts to arrays
    asy_att_e_list <- lapply(onecase, `[[`, "asy_att_e")
    asy_att_e_pad <- lapply(asy_att_e_list, pad_vector, nrow = S)
    asy_att_e_array <- t(simplify2array(asy_att_e_pad))
    dimnames(asy_att_e_array) <- list(sim = sim_names,
                                      event = event_names)

    asy_time_list  <- lapply(onecase, `[[`, "asy_time")
    asy_time_array <- do.call(rbind, asy_time_list)
    rownames(asy_time_array) <- sim_names

    # Store results
    results[['asy']][[case_name]] <- list(
    att_e = asy_att_e_array,
    time  = asy_time_array
    )
  }

  if (did==TRUE){
    # extract results and converts to arrays
    did_att_e_list <- lapply(onecase, `[[`, "did_att_e")
    did_att_e_pad <- lapply(did_att_e_list, pad_vector, nrow = S)
    did_att_e_array <- t(simplify2array(did_att_e_pad))
    dimnames(did_att_e_array) <- list(sim = sim_names,
                                      event = event_names)

    did_time_list  <- lapply(onecase, `[[`, "did_time")
    did_time_array <- do.call(rbind, did_time_list)
    rownames(did_time_array) <- sim_names

    # Store results
    results[['did']][[case_name]] <- list(
    att_e = did_att_e_array,
    time  = did_time_array
    )
  }
  

  # Intermediate save after each case for safety
  # if (!dir.exists("output")) {
  #   dir.create("output") # create output folder if it does not exist
  # }
  # save(results,file="output/sim_result_alternative_methods.RData")
  
} # end of all cases
stopCluster(cl) # stop parallel computing
cat("\nEstimation completed.")
cat("\nRun time: ");print(Sys.time()-begin.time)
cat("\nRun the code block following to save results\n")


# %% Save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# results would be saved in "output" folder: 
# one csv file per case-method, one row per simulation

if (!dir.exists("output")) {
  dir.create("output") # create output folder if it does not exist
}

for (case in 1:ncases) {
  case_name <- WW[case]
  # if (!dir.exists(paste0("output/", case_name))) {
  #   dir.create(paste0("output/", case_name)) # create output folder if it does not exist
  # }

  for(method in names(results)){
    method_results <- results[[method]]
    att_e_array <- method_results[[case_name]]$att_e
    time_array  <- method_results[[case_name]]$time
    result_array = cbind(att_e_array, time_array)
    result_array <- as.data.frame(result_array)
    colnames(result_array)[1:ncol(att_e_array)] <- paste('event_time', colnames(att_e_array), sep = "_")
    colnames(result_array)[(ncol(att_e_array)+1):ncol(result_array)] <- c("time_sec")

    # save results: one folder per case, one file per method, one row per simulation
    # output_filename <- paste0("output/", case_name, "/", method , "_results_" , case_name, ".csv")
    files_previous <- list.files("output/", pattern = method) # delete previous results containing this method for any case
    if (length(files_previous) > 0 & case == 1) file.remove(paste0("output/", files_previous)) 

    output_filename <- paste0("output/sim_results_", case_name, "_", method , ".csv")
    write.csv(result_array, file = output_filename, row.names = FALSE)
  }
}

cat("\nSimulation results of alternative methods saved in 'output' folder.\n")




