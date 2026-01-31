## Simulation Study
## Replication Materials


## This scripts produce simulation results of alternative methods for event-time ATT
## Alternative methods include:
## 1. Generalized Synthetic Control (gsc) by Xu (2017)
## 2. Augmented Synthetic Control (asy) by Ben-Michael, Feller, and Rothstein (2021)
## 3. Difference-in-Differences (did) method by Callaway and Sant'Anna (2021)


# %% set environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list=ls(all=TRUE))
set.seed(123)
setwd("/Users/wuhang/Desktop/metric/sc/script/simulation")


library(dplyr)
library(doParallel)
library(parallel)
library(foreach)
library(gsynth)
library(augsynth)
library(did)



# %% pre-set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# load simulated data
load("data/sim_data.RData")
# select test cases
# results_list <- results_list[c("N15_T15")]

# for each model, consider 9 cases: varying N and T, fix S=10
WW  = names(sim_data) # names of cases
ncases<-length(WW)

S <- 10 # post-treatment periods
nsims <- length(sim_data[[1]]) # number of simulations per case
# nsims <-2 # set small number for testing

# methods to compare:
gsc = TRUE  # Generalized Synthetic Control
asy = TRUE  # Augmented Synthetic Control
did = FALSE  # Callaway and Sant'Anna (2021) DID method





# %% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# padding function for to unify matrix sizes: due to inconsistent event-time lengths
pad_matrix <- function(mat, nrow, ncol) {
  padded_mat <- matrix(NA, nrow = nrow, ncol = ncol)
  padded_mat[1:nrow(mat), 1:ncol(mat)] <- mat
  return(padded_mat)
}




# %%Start estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

# register multiple cores
#cores<-detectCores()
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)
begin.time<-Sys.time()
cat("Cores: ",cores,"; # of cases: ",ncases,"; # of sims: ",nsims,sep="")

# loop
for (case in 1:ncases) {
    case_name <- WW[case]
    T <- as.numeric(sub(".*T([0-9]+).*", "\\1", case_name))
    N <- as.numeric(sub("N([0-9]+)_.*", "\\1", case_name))
    ## annouce case
    cat("\nCase: ", case_name, "\n",sep="")
    
    # extract case data
    data_case <- sim_data[[case_name]]

    # start simulations for one case
    onecase<-foreach (i=1:nsims,.combine='list', .multicombine=TRUE,.inorder=FALSE,
                      .packages=c("gsynth", "augsynth", "did","dplyr"), .verbose = FALSE) %dopar% { # parallel computing
        
        ## annouce simulations
        if (i %% 200 == 0) {
          cat(sprintf("  [Case %d/%d] Simulation %d/%d completed...\n", case, ncases, i, length(data_case)))
          flush.console()
        }

        ## Prepare data:
        panel <- data_case[[i]]
        required_cols <- c('id', 'time', 'Y', 'D','att')
        panel <- panel[, required_cols]
        att_true <- unique(panel$att)


        ## GSC
        if(gsc==TRUE){
        t0_gsc <- Sys.time()
        out_gsyn <- gsynth(Y ~ D, data = panel, EM = TRUE, 
                        index=c("id","time"), inference="parametric",
                        se = TRUE, nboots = 1000, r = 2, CV = FALSE, # the number of factors are set to be 2
                        force = "two-way", parallel = TRUE)
                
        t1_gsc <- Sys.time()
        time_gsc <- as.numeric(difftime(t1_gsc, t0_gsc)) # used time for 1 simulation 

        att_e.gsy = out_gsyn$est.att # event time
        att_e.gsy = att_e.gsy[as.numeric(rownames(att_e.gsy))>0,c('ATT', 'CI.lower','CI.upper') ]
        att_e.gsy <- cbind(eff = 1:nrow(att_e.gsy), att_e.gsy)
        att_e.gsy <- cbind(att_e.gsy, bias = att_e.gsy[, "ATT"] - att_e.gsy[, "eff"])
        att_e.gsy <- cbind(att_e.gsy, sq_error = (att_e.gsy[, "ATT"] - att_e.gsy[, "eff"])^2)
        att_e.gsy <- cbind(att_e.gsy, in_CI = (att_e.gsy[, "CI.lower"] <= att_e.gsy[, "eff"]) & (att_e.gsy[, "CI.upper"] >= att_e.gsy[, "eff"]))
        att_e.gsy <- att_e.gsy[, c("eff", "ATT", "bias", "sq_error", "CI.lower", "CI.upper", "in_CI")]

        } else {
          att_e.gsy <- NULL
          time_gsc <- NULL
        }

        ## ASY: do not calculate ATT, as they only calulate att_e for the shortest treated periods across all treated units.
        if(asy==TRUE){
        t0_asy <- Sys.time()   
        
        ppool_syn <- multisynth(Y ~ D, unit = id, time = time,
                        panel)            
        t1_asy <- Sys.time()
        time_asy <- as.numeric(difftime(t1_asy, t0_asy))
        
        att_e.asy = summary(ppool_syn)
        att_e.asy <- as.data.frame(att_e.asy$att)
        att_e.asy = att_e.asy[att_e.asy$Level =='Average' & is.na(att_e.asy$Time) ==FALSE & att_e.asy$Time>=0, ] # only post-treatment periods for shortest treated periods
        cols <- c('Time', 'Estimate', 'lower_bound', 'upper_bound')
        att_e.asy <- att_e.asy[, cols]
        att_e.asy[] <- lapply(att_e.asy, as.numeric)
        att_e.asy <- within(att_e.asy, {
          eff <- Time + 1
          sq_error <- (Estimate - eff)^2
          in_CI <- (lower_bound <= eff) & (upper_bound >= eff)
        })
        att_e.asy$bias <- att_e.asy$Estimate - att_e.asy$eff
        att_e.asy <- att_e.asy[, c("eff", "Estimate", "bias", "sq_error", "lower_bound", "upper_bound", "in_CI")]
        att_e.asy <- as.matrix(att_e.asy)

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
        time_did <- as.numeric(difftime(t1_did, t0_did))

        att_e.did <- data.frame(
          att = did_summary_att_e$att.egt,
          time = did_summary_att_e$egt,
          CI.lower = did_summary_att_e$att.egt - did_summary_att_e$crit.val.egt * did_summary_att_e$se.egt,
          CI.upper = did_summary_att_e$att.egt + did_summary_att_e$crit.val.egt * did_summary_att_e$se.egt
        )
        att_e.did <- att_e.did[att_e.did$time >= 0, ]  # extract only post-treatment periods ATT
        att_e.did$eff <- att_e.did$time + 1
        att_e.did <- within(att_e.did, {
          bias <- att - eff
          sq_error <- (att - eff)^2
          in_CI <- (CI.lower <= eff) & (CI.upper >= eff)
        })
        att_e.did <- att_e.did[, c("eff", "att", "bias", "sq_error", "CI.lower", "CI.upper", "in_CI")]
        att_e.did <- as.matrix(att_e.did)
        
        } else {
          att_e.did <- NULL
          time_did <- NULL
        }

        out = list(
          gsc_att_e = att_e.gsy,     # matrix
          gsc_time  = time_gsc,

          asy_att_e = att_e.asy,     # matrix
          asy_time  = time_asy,

          did_att_e = att_e.did,     # matrix
          did_time  = time_did
        )
        return(out)      
  }
  # collect results and store in arrays : dim =  event-time; stats; sims
  stat_names  <- c("eff", "ATT", "bias", "sq_error", "CI.lower", "CI.upper", "in_CI")
  event_names <- c(1:S)
  sim_names   <- c(1:nsims)

  if (gsc==TRUE){
    # extract results and converts to arrays
    gsc_att_e_list <- lapply(onecase, `[[`, "gsc_att_e")
    gsc_att_e_pad <- lapply(gsc_att_e_list, pad_matrix, nrow = S, ncol = 7)
    gsc_att_e_array <- simplify2array(gsc_att_e_pad)
    dimnames(gsc_att_e_array) <- list(event = event_names,
                                      stat = stat_names,
                                      sim = sim_names)

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
    asy_att_e_pad <- lapply(asy_att_e_list, pad_matrix, nrow = S, ncol = 7)
    asy_att_e_array <- simplify2array(asy_att_e_pad)
    dimnames(asy_att_e_array) <- list(event = event_names,
                                      stat = stat_names,
                                      sim = sim_names)

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
    did_att_e_pad <- lapply(did_att_e_list, pad_matrix, nrow = S, ncol = 7)
    did_att_e_array <- simplify2array(did_att_e_pad)
    dimnames(did_att_e_array) <- list(event = event_names,
                                      stat = stat_names,
                                      sim = sim_names)

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
cat("\n");print(Sys.time()-begin.time)



# %% Results summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# summary would be stored in a nested list: method -> case -> data.frame of summary statistics


# create a list for summary storage
summary = list()
for (method in names(results)) {
  summary[[method]] <- list()
  method_results <- results[[method]]
  for (case_name in names(method_results)) {
    att_e_array <- method_results[[case_name]]$att_e
    # Calculate summary statistics across simulations
    bias_att_e <- apply(att_e_array[,"bias", ], 1, mean, na.rm = TRUE)
    var_att_e   <- apply(att_e_array[, "ATT", ], 1, var, na.rm = TRUE)
    mse_att_e <- apply(att_e_array[, "sq_error", ], 1, mean, na.rm = TRUE)
    coverage_att_e <- apply(att_e_array[, "in_CI", ], 1, mean, na.rm = TRUE)
    elapsed_time <- sum(method_results[[case_name]]$time, na.rm = TRUE); elapsed_time <- rep(elapsed_time, times = nrow(att_e_array)); elapsed_time[is.na(coverage_att_e)] <- NA
    

    # Print or save the summary statistics as needed
    summary_df <- data.frame(
      event_time = as.numeric(rownames(att_e_array)),
      bias = bias_att_e,
      var = var_att_e,
      mse = mse_att_e,
      coverage_rate = coverage_att_e,
      `time(s)`  = elapsed_time,
      check.names = FALSE
    )
    summary[[method]][[case_name]] <- summary_df
    print(paste("Method:", method, "Case:", case_name))
    print(head(summary_df,2))
  }
}


# %% Save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (!dir.exists("output")) {
  dir.create("output") # create output folder if it does not exist
}

for (case in 1:ncases) {
  case_name <- WW[case]
  summary_all_methods <- NULL

  # load SSC results for combination
  filename_ssc <- paste0("output/summary_ssc_" , case_name, ".csv")
  if (file.exists(filename_ssc)) {
    summary_all_methods <- read.csv(filename_ssc, check.names = FALSE)
    colnames(summary_all_methods)[2:ncol(summary_all_methods)] <- paste('ssc',colnames(summary_all_methods)[2:ncol(summary_all_methods)],sep = "_")
  }

   # extract results of alternative methods and combine 
  for(method in names(summary)){
    summary_df <- summary[[method]][[case_name]][, -1] # remove event_time column to avoid duplication
    colnames(summary_df) <- paste(method,colnames(summary_df),sep = "_")
    if (is.null(summary_all_methods)) {
      summary_all_methods <- summary_df
    }else{ summary_all_methods <- cbind(summary_all_methods, summary_df) }
  }

  # Save results across all methods for one case
  output_filename <- paste0("output/summary_all_methods_" , case_name, ".csv")
  write.csv(summary_all_methods, file = output_filename, row.names = FALSE)
}




