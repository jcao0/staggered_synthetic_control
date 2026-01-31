## Simulation Study
## Replication Materials

# This scripts define functions to generate simulated samples

## Function to determin hyperparameters for treatment assignment
## Function to create treatment matrix
## Function to create simulated samples (Revised based on  Yiqing Xu (2017)


# function to get logistic parameters for treatment assignment--------
get_logistic_theta <- function(r_treat, r_treat_period, seed=NULL) {
  #' Estimate logistic model parameters for treatment assignment
  #'
  #' This function estimates the parameters of a logistic model for treatment assignment,
  #' given the desired proportion of treated units and the proportion treated per period.
  #'
  #' @param r_treat Numeric. Desired proportion of units eventually treated.
  #' @param r_treat_period Numeric. Desired proportion of treated units per period.
  #' @return Numeric vector of estimated logistic parameters (theta).
  #' @examples
  #' theta <- get_logistic_theta(0.5, 0.1)
    if (!is.null(seed)) {
        set.seed(seed)
    }
  mu_grand <- r_treat       # chosen mean of negative binomial with p(Z), percent of units eventually treated
  mu <- r_treat_period      # chosen mean of p(Z), percent of treated units per period
  M <- 5e5
  lambda <- array(runif(3 * M, min = -sqrt(3), max = sqrt(3)), dim = c(3, M))
  lambda <- apply(lambda, c(2), sum) # returns an N x 1000 matrix, sum over columns for each simulation
  Z <- lambda
  F <- function(theta) {
    a <- theta[1]
    b <- theta[2]
    pZ <- plogis(a + b * Z)
    c(
      mean((1 - pZ)^10) - (1 - mu_grand),
      mean(pZ) - mu
    )
  }
  # Initial guess
  start <- c(0, 0)
  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    install.packages("nleqslv")
  }
  library(nleqslv)
  out <- nleqslv::nleqslv(start, F)
  return(out$x)
}

# function to generate treatment matrix A --------
get_treatment_matrix <- function(N, S, theta=c(-0.0026, 1.037), r=0, seed=NULL) {
    #' Generate treatment assignment matrix based on logistic model parameters
    #'
    #' @param N Numeric. Number of units.
    #' @param S Numeric. Number of post-treatment periods.
    #' @param theta Numeric vector. Logistic model parameters for treatment assignment. Default corresponds to approx. 90% units eventually treated and 50% treated per period.
    #' @param r Numeric. Number of latent factors influencing treatment assignment.
    #' @param seed Numeric. Random seed for reproducibility.
    #' @return Treatment assignment matrix A of dimension (S x N).
    #' @examples
    #' A <- get_treatment_matrix(N=20, S=10, theta=c(-1.7,-0.6), r=2)
    
    if (!is.null(seed)) {
        set.seed(seed)
    }
    # generate latent factors for treatment assignment
    ss <- sqrt(3) # to ensure variance =1
    lambda<-matrix(runif(N*r,min=-ss,max=ss),N,r)  # loadings

    # generate treatment matrix N#S
    A <-array(0,dim=c(S,N))    # regressor matrix, must be T_total by N  by  p 
    for (s in 1:S) {
        pi = plogis(theta[1] + theta[2] * rowSums(lambda))
        if(s==1){
            A[s,] <- rbinom(N,1,pi)
        } else {
            A[s,] <- ifelse(A[s-1,]==1, 1, rbinom(N,1,pi))} # once treated, always treated
    }
    return(A)
}

# function to generate simulated data)--------
simulate_data_fm <- function(N, T, S, te=NULL, D.sd=0, A=NULL, theta = c(-1.7,-0.6),
                   p=0, beta=NULL, fsize=1, r=0, FE=FALSE, AR1=0, mu=0, fixF=FALSE, fixL=FALSE, F=NULL, L=FULL,
                    seed=NULL){
    ## N：number of units
    ## T, S: pre-treatment and post-treatment periods
    ## te: event-time treatment effect, NULL means effect increases by 1 each period
    ## D.sd: treatment effect heterogeneity
    ## A: treatment design matrix (S by N), if NULL, generate according to a logit model
    ## theta: parameters of logit model for the treatment assignment model. Defalt c(-1.7,-0.6) corresponds to about 75% units treated eventually and 20% treated each period
    ## p: number of covariates (p=0, no covariates)
    ## beta: true coefficients for covariates X
    ## fsize: relative importance of factors vs. covariates
    ## r: number of unobserved factors
    ## FE: to include unit and time fixed effects
    ## AR1: Autoregressive(1) coefficient, apply for factors and time effects
    ## mu: grand mean
    ## fixF: factors as given
    ## fixL: loadings as given
    ## F: given factors if fixF=TRUE
    ## L: given loadings if fixL=TRUE
    ## seed: random seed


    if (is.null(seed)==FALSE) {
        set.seed(seed)
    }
    
    N<-N
    T_total<-T+S
    ## ###########################
    ## Data generating process
    ## ###########################
    
    rr<-r #  number of unobserved factors

    ## loadings
    ss<-sqrt(3) # to ensure variance =1
    if (rr > 0) {
        ## loadings
        if (fixL==FALSE) { # L not given
            lambda<-matrix(runif(N*rr,min=-ss,max=ss),N,rr)  # loadings     
        } else {lambda<-L}
        
        ## factors
        if (fixF==FALSE) { # F not given
            factor<-matrix(rnorm(T_total*rr),T_total,rr)   # factors
            if(AR1>0){
                factor <- apply(factor, 2, stats::filter, filter = AR1, method = "recursive")
            }
        } else {factor<-F}
    }

    ## fixed effects
    if (FE==TRUE) {
        ## unit fixed effects
        if (fixL==FALSE) {
            alpha<-runif(N,min=-ss,max=ss)
        } 

        ## time fixed effects
        if (fixF==FALSE) {
            xi<-rnorm(T_total,0,1) 
            if(AR1>0){
                xi <- stats::filter(xi, filter = AR1, method = "recursive")
            }
        }   
    }

    ## error
    e <- matrix(rnorm(T_total*N),T_total,N) # disturbances

    ## time varying covariates: always generated by the same first factors
    truebeta<-beta
    if (p!=0) {
        X<-array(0,dim=c(T_total,N,p))    # regressor matrix, must be T_total by N  by  p 
        for (j in 1:p) {
            X[,,j] <- matrix(rnorm(T_total*N),T_total,N)  +
                0.5 * factor[,1] %*% T_total(lambda[,1]) +
                0.25 * matrix(1,T_total,1) %*% T_total(lambda[,1]) +
                0.25 * factor[,1] %*% matrix(1,1,N)+1
        }
    }  

    ## treatment assignment
    D_pre <- matrix(0, nrow = T, ncol = N)  # pre-treatment matrix 
    if (is.null( A )==TRUE) {
        if(theta==NULL){
            theta = c(-1.7,-0.6) # default parameters, corresponds to about 75% units treated eventually and 20% treated each period
        }
        l <- r + as.numeric(FE)  # number of latent factors for treatment assignment
        A = get_treatment_matrix(N=N, S=S, theta = theta, r=l)    # treatment matrix, S by N
        }
    D = rbind(D_pre, A )  # full treatment matrix
    N_tr<- sum(colSums(D)>0) # number of treated units

    ## event-time treatment indicator
    D.event<-matrix(0,T_total,N)
    for (s in 1:S) {
        D.event[T+s,]<-D.event[(T+s-1),]+D[T+s,]
    }

    ##  event-time treatment effect
    if(is.null(te)==TRUE) {
        te<-c(1:S) # effect increases by 1 each period
    }
    eff <- matrix(0, T_total, N)
    for (s in 1:S) {
    eff[T + s, D.event[T + s, ] > 0] <- te[D.event[T + s, D.event[T + s, ] > 0]] + rnorm(sum(D.event[T + s, ] > 0),0,D.sd)
    }
    att <- mean(eff[eff != 0]) # sampling ATT 


    ## outcome variable

    Y0<- e  +  matrix(mu,T_total,N) # error + grand mean  
    if (r>0) {
        Y0 <- Y0 + fsize*factor%*%t(lambda)
    } 
    if (FE==TRUE) {
        Y0 <- Y0 + matrix(alpha,T_total,N,byrow=TRUE) + matrix(xi,T_total,N,byrow=FALSE)
    }
    if (p!=0) { # covariates
        for (k in 1:p) {
            Y0<-Y0+X[,,k]*truebeta[k]
        }
    }
    Y1 <- Y0 + eff
 
    Y<-(matrix(1,T_total,N)-D)*Y0+D*Y1 
    # if (AR1>0) {
    #     Y.lag<-matrix(NA,T_total,N); Y.lag[2:T_total,]<-Y[1:(T_total-1),]
    # }

    ## substract error
    Y.bar <- Y - e

    ## panel structure
    panel<-as.data.frame(cbind(rep(101:(100+N),each=T_total),rep(1:T_total,N),
                               c(Y),c(Y0),c(Y1),c(D),c(eff),c(e),c(Y.bar),
                               rep(mu,T_total*N),rep(1,T_total*N),rep(N_tr,T_total*N),rep(att,T_total*N)))
    cname<-c("id","time","Y","Y0","Y1","D","eff","error","Ybar","mu","X0","N_tr","att")
    if (p!=0) {
        for (i in 1:p) {
            panel<-cbind(panel,c(X[,,i]))
            cname<-c(cname,paste("X",i,sep=""))
        }
    }
    if (rr > 0) {
        for (i in 1:rr) {
            panel<-cbind(panel,rep(factor[,i],N))
            cname<-c(cname,paste("F",i,sep=""))
        }
    }
    if (rr > 0) {
        for (i in 1:rr) {
            panel<-cbind(panel,rep(lambda[,i],each=T_total))
            cname<-c(cname,paste("L",i,sep=""))
        }
    } 
    if (FE==TRUE) {
        panel<-cbind(panel,rep(alpha,each=T_total),rep(xi,N))
        cname<-c(cname,"alpha","xi")
    }  
    # if (AR1>0) {
    #     panel<-cbind(panel,c(Y.lag),c(eff.acc))
    #     cname<-c(cname,"Y_lag","eff_acc")
    # } 
    colnames(panel)<-cname
    return(panel)
}

# example usage:
# theta <- get_logistic_theta(r_treat = 0.8, r_treat_period = 0.5)
# A <- get_treatment_matrix(N=20, S=10, theta=c(0.0368, 1.725), r=2)
# sim_data <- simulate_data_fm(N=20, T=5, S=10, p=0, r=2, AR1=0.5, D.sd=0, A=A)