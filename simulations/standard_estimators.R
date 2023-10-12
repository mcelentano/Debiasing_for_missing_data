# run from the simulations directory

library(glmnet)
library(boot)
library(parallel)

source("../utils/standard_estimators_utils.R")

run.replicate <- function(n,p,theta.d.0,theta.y.0,theta.d,theta.y,epsilon,seed,filename.tracking) {
  
  # set seed
  set.seed(seed)
  cat(n,seed,"\n",file=filename.tracking,append=TRUE)
  
  # Generate data
  X = matrix(rnorm(n*p),nrow=n,ncol=p)
  z.y = rnorm(n)
  y.lin = theta.y.0 + X %*% theta.y + z.y  
  y.nlin = theta.y.0 + X %*% theta.y + (X[,1]^2 - 1) + z.y  # subtract 1 so mean of nonlinear term is 0
  d = rbinom( n , 1 , epsilon + (1-epsilon) * inv.logit(theta.d.0 + X %*% theta.d) )
  
  # Data splits
  folds.cf1 <- 
    list(
      1:n
    )
  folds.cf2 <- 
    list(
      1:round(n/2),
      (round(n/2)+1):n
    )
  folds.cf3 <- 
    list(
      1:round(n/3),
      (round(n/3) + 1):round(2*n/3),
      (round(2*n/3)+1):n
    )
  
  # Models
  lin.out.models <- 
    lapply(
      list(folds.cf1,folds.cf2,folds.cf3),
      function(folds) {
        lapply(
          folds,
          function(fold) {
            fold.obs <- fold[as.logical(d[fold])]
            model = lm(y.lin[fold.obs] ~ 1 + X[fold.obs,])
            return( coef(model) )  
          }
        )
      }
      )
  lin.out.reg.models <- 
    lapply(
      list(folds.cf1,folds.cf2,folds.cf3),
      function(folds) {
        lapply(
          folds,
          function(fold) {
            fold.obs <- fold[as.logical(d[fold])]
            model = cv.glmnet(
              X[fold.obs,],y.lin[fold.obs],
              family = "gaussian",
              alpha = 0, # ridge
              standardize = FALSE
            )
            return( coef(model) )  
          }
        )
      }
    )
  nlin.out.models <- 
    lapply(
      list(folds.cf1,folds.cf2,folds.cf3),
      function(folds) {
        lapply(
          folds,
          function(fold) {
            fold.obs <- fold[as.logical(d[fold])]
            model = lm(y.nlin[fold.obs] ~ 1 + X[fold.obs,])
            return( coef(model) )
          }
        )
      }
    )
  prop.models <- 
    lapply(
      list(folds.cf1,folds.cf2,folds.cf3),
      function(folds) {
        lapply(
          folds,
          function(fold) {
            model = glm(d[fold] ~ 1 + X[fold,],family = "binomial")
            return( coef(model) )
          }
        )
      }
    )
  
  # Some variable definitions to help with readability
  cf1.idx <- 1
  cf2.idx <- 2
  cf3.idx <- 3
  
  ##### Various estimates #####
  
  # Intercept of well-specified model
  lin.out.intercept <- 
    lin.out.models[[cf1.idx]][[1]][1]
  
  # G-computation estimates without regularization, well-specified
  lin.out.cf1.g.est <- # G computation of well-specified model without sampling splitting
    g.computation(X,
                  out.models=lin.out.models[[cf1.idx]],g.folds=folds.cf1,
                  out.idxs=c(1),g.idxs=c(1))
  lin.out.cf2.g.est <- # G computation of well-specified model with sampling splitting
    g.computation(X,
                  out.models=lin.out.models[[cf2.idx]],g.folds=folds.cf2,
                  out.idxs=c(1,2),g.idxs=c(2,1))
  
  # G-computation estimates without regularization, misspecified
  nlin.out.cf1.g.est <- # G computation of well-specified model without sampling splitting
    g.computation(X,
                  out.models=nlin.out.models[[cf1.idx]],g.folds=folds.cf1,
                  out.idxs=c(1),g.idxs=c(1))
  nlin.out.cf2.g.est <- # G computation of well-specified model with sampling splitting
    g.computation(X,
                  out.models=nlin.out.models[[cf2.idx]],g.folds=folds.cf2,
                  out.idxs=c(1,2),g.idxs=c(2,1))
  
  # G-computation estimates with regularization
  lin.out.cf1.g.reg.est <- # G computation of well-specified model without sampling splitting
    g.computation(X,
                  out.models=lin.out.reg.models[[cf1.idx]],g.folds=folds.cf1,
                  out.idxs=c(1),g.idxs=c(1))
  lin.out.cf2.g.reg.est <- # G computation of well-specified model with sampling splitting
    g.computation(X,
                  out.models=lin.out.reg.models[[cf2.idx]],g.folds=folds.cf2,
                  out.idxs=c(1,2),g.idxs=c(2,1))
  
  # AIPW estimates with estimated propensity score and well specified outcome
  aipw.est_prop.lin_out.cf1.est <- # AIPW for well-specified model without sampling splitting
    aipw(y.lin,d,X,epsilon,
         out.models=lin.out.models[[cf1.idx]],
         prop.models=prop.models[[cf1.idx]],
         aipw.fold=folds.cf1,
         out.idxs=c(1),prop.idxs=c(1),aipw.idxs=c(1))
  aipw.est_prop.lin_out.cf2.est <- # AIPW for well-specified model, outcome and propensity fit on same fold, aipw terms on second fold
    aipw(y.lin,d,X,epsilon,
         out.models=lin.out.models[[cf2.idx]],
         prop.models=prop.models[[cf2.idx]],
         aipw.fold=folds.cf2,
         out.idxs=c(1,2),prop.idxs=c(1,2),aipw.idxs=c(2,1))
  aipw.est_prop.lin_out.cf3.est <- # AIPW for well-specified model, outcome fit, propensity fit, nd aipw terms on all different folds
    aipw(y.lin,d,X,epsilon,
         out.models=lin.out.models[[cf3.idx]],
         prop.models=prop.models[[cf3.idx]],
         aipw.fold=folds.cf3,
         out.idxs=c(1,1,2,2,3,3),prop.idxs=c(2,3,1,3,1,2),aipw.idxs=c(3,2,3,1,2,1))
  
  # AIPW estimates with estimated propensity score and misspecified outcome
  aipw.est_prop.nlin_out.cf1.est <- # AIPW for well-specified model without sampling splitting
    aipw(y.nlin,d,X,epsilon,
         out.models=nlin.out.models[[cf1.idx]],
         prop.models=prop.models[[cf1.idx]],
         aipw.fold=folds.cf1,
         out.idxs=c(1),prop.idxs=c(1),aipw.idxs=c(1))
  aipw.est_prop.nlin_out.cf2.est <- # AIPW for well-specified model, outcome and propensity fit on same fold, aipw terms on second fold
    aipw(y.nlin,d,X,epsilon,
         out.models=nlin.out.models[[cf2.idx]],
         prop.models=prop.models[[cf2.idx]],
         aipw.fold=folds.cf2,
         out.idxs=c(1,2),prop.idxs=c(1,2),aipw.idxs=c(2,1))
  aipw.est_prop.nlin_out.cf3.est <- # AIPW for well-specified model, outcome fit, propensity fit, nd aipw terms on all different folds
    aipw(y.nlin,d,X,epsilon,
         out.models=nlin.out.models[[cf3.idx]],
         prop.models=prop.models[[cf3.idx]],
         aipw.fold=folds.cf3,
         out.idxs=c(1,1,2,2,3,3),prop.idxs=c(2,3,1,3,1,2),aipw.idxs=c(3,2,3,1,2,1))
  
  # AIPW estimates with estimated propensity score, well-specified outcome, and outcome regularization
  aipw.est_prop.reg_lin_out.cf1.est <- # AIPW for well-specified model without sampling splitting
    aipw(y.lin,d,X,epsilon,
         out.models=lin.out.reg.models[[cf1.idx]],
         prop.models=prop.models[[cf1.idx]],
         aipw.fold=folds.cf1,
         out.idxs=c(1),prop.idxs=c(1),aipw.idxs=c(1))
  aipw.est_prop.reg_lin_out.cf2.est <- # AIPW for well-specified model, outcome and propensity fit on same fold, aipw terms on second fold
    aipw(y.lin,d,X,epsilon,
         out.models=lin.out.reg.models[[cf2.idx]],
         prop.models=prop.models[[cf2.idx]],
         aipw.fold=folds.cf2,
         out.idxs=c(1,2),prop.idxs=c(1,2),aipw.idxs=c(2,1))
  aipw.est_prop.reg_lin_out.cf3.est <- # AIPW for well-specified model, outcome fit, propensity fit, nd aipw terms on all different folds
    aipw(y.lin,d,X,epsilon,
         out.models=lin.out.reg.models[[cf3.idx]],
         prop.models=prop.models[[cf3.idx]],
         aipw.fold=folds.cf3,
         out.idxs=c(1,1,2,2,3,3),prop.idxs=c(2,3,1,3,1,2),aipw.idxs=c(3,2,3,1,2,1))
  
  # AIPW estimates with oracle propensity score
  aipw.oracle_prop.lin_out.cf1.est <- # AIPW for well-specified model, oracle propensity fit, outcome fit and aipw terms same fold
    aipw(y.lin,d,X,epsilon,
         out.models=lin.out.models[[cf1.idx]],
         prop.models=list(c(theta.d.0,theta.d)),
         aipw.folds=folds.cf1,
         out.idxs=c(1),prop.idxs=c(1),aipw.idxs=c(1))
  aipw.oracle_prop.nlin_out.cf1.est <- # AIPW for misspecified model, oracle propensity fit, outcome fit and aipw terms on same fold
    aipw(y.nlin,d,X,epsilon,
         out.models=nlin.out.models[[cf1.idx]],
         prop.models=list(c(theta.d.0,theta.d)),
         aipw.folds=folds.cf1,
         out.idxs=c(1),prop.idxs=c(1),aipw.idxs=c(1))
  aipw.oracle_prop.lin_out.cf2.est <- # AIPW for well-specified model, oracle propensity fit, outcome fit and aipw terms on different folds
    aipw(y.lin,d,X,epsilon,
         out.models=lin.out.models[[cf2.idx]],
         prop.models=list(c(theta.d.0,theta.d)),
         aipw.folds=folds.cf2,
         out.idxs=c(1,2),prop.idxs=c(1,1),aipw.idxs=c(2,1))
  aipw.oracle_prop.nlin_out.cf2.est <- # AIPW for misspecified model, oracle propensity fit, outcome fit and aipw terms on different folds
    aipw(y.nlin,d,X,epsilon,
         out.models=nlin.out.models[[cf2.idx]],
         prop.models=list(c(theta.d.0,theta.d)),
         aipw.folds=folds.cf2,
         out.idxs=c(1,2),prop.idxs=c(1,1),aipw.idxs=c(2,1))
  aipw.oracle_prop.oracle_out.nlin_out.cf1.est <- # AIPW for misspecified model, oracle propensity fit, oracled outcome fit, single fold
    aipw(y.nlin,d,X,epsilon,
         out.models=list(c(theta.y.0,theta.y)),
         prop.models=list(c(theta.d.0,theta.d)),
         aipw.folds=folds.cf1,
         out.idxs=c(1),prop.idxs=c(1),aipw.idxs=c(1))
  
  # IPW estimates with estimated propensity scores
  ipw.est_prop.lin_out.cf1.est <- # IPW for well-specified model, estimated propensity fit, propensity fit and ipw terms on same fold
    ipw(y.lin,d,X,epsilon,
        prop.models=prop.models[[cf1.idx]],
        ipw.folds=folds.cf1,
        prop.idxs=c(1),ipw.idxs=c(1))
  ipw.est_prop.lin_out.cf2.est <- # IPW for well-specified model, estimated propensity fit, propensity fit and ipw terms on different folds
    ipw(y.lin,d,X,epsilon,
        prop.models=prop.models[[cf2.idx]],
        ipw.folds=folds.cf2,
        prop.idxs=c(1,2),ipw.idxs=c(2,1))
  ipw.oracle_prop.lin_out.est <- # IPW for well-specified model, oracle propensity
    ipw(y.lin,d,X,epsilon,
        prop.models=list(c(theta.d.0,theta.d)),
        ipw.folds=folds.cf1,
        prop.idxs=c(1),ipw.idxs=c(1))
  ipw.oracle_prop.nlin_out.est <- # IPW for misspecified model, oracle propensity
    ipw(y.nlin,d,X,epsilon,
        prop.models=list(c(theta.d.0,theta.d)),
        ipw.folds=folds.cf1,
        prop.idxs=c(1),ipw.idxs=c(1))
  
  return(
    data.frame(
      
      seed = seed,
      n = n,
      p = p,
      
      lin.out.intercept = lin.out.intercept,
      
      lin.out.cf1.g.est = lin.out.cf1.g.est,
      lin.out.cf2.g.est = lin.out.cf2.g.est,
      
      nlin.out.cf1.g.est = nlin.out.cf1.g.est,
      nlin.out.cf2.g.est = nlin.out.cf2.g.est,
      
      lin.out.cf1.g.reg.est = lin.out.cf1.g.reg.est,
      lin.out.cf2.g.reg.est = lin.out.cf2.g.reg.est,
      
      aipw.est_prop.lin_out.cf1.est = aipw.est_prop.lin_out.cf1.est,
      aipw.est_prop.lin_out.cf2.est = aipw.est_prop.lin_out.cf2.est,
      aipw.est_prop.lin_out.cf3.est = aipw.est_prop.lin_out.cf3.est,
      
      aipw.est_prop.nlin_out.cf1.est = aipw.est_prop.nlin_out.cf1.est,
      aipw.est_prop.nlin_out.cf2.est = aipw.est_prop.nlin_out.cf2.est,
      aipw.est_prop.nlin_out.cf3.est = aipw.est_prop.nlin_out.cf3.est,
      
      aipw.est_prop.reg_lin_out.cf1.est = aipw.est_prop.reg_lin_out.cf1.est,
      aipw.est_prop.reg_lin_out.cf2.est = aipw.est_prop.reg_lin_out.cf2.est,
      aipw.est_prop.reg_lin_out.cf3.est = aipw.est_prop.reg_lin_out.cf3.est,
      
      aipw.oracle_prop.lin_out.cf1.est = aipw.oracle_prop.lin_out.cf1.est,
      aipw.oracle_prop.nlin_out.cf1.est = aipw.oracle_prop.nlin_out.cf1.est,
      aipw.oracle_prop.lin_out.cf2.est = aipw.oracle_prop.lin_out.cf2.est,
      aipw.oracle_prop.nlin_out.cf2.est = aipw.oracle_prop.nlin_out.cf2.est,
      aipw.oracle_prop.oracle_out.nlin_out.cf1.est = aipw.oracle_prop.oracle_out.nlin_out.cf1.est,
      
      ipw.est_prop.lin_out.cf1.est = ipw.est_prop.lin_out.cf1.est,
      ipw.est_prop.lin_out.cf2.est = ipw.est_prop.lin_out.cf2.est,
      ipw.oracle_prop.lin_out.est = ipw.oracle_prop.lin_out.est,
      ipw.oracle_prop.nlin_out.est = ipw.oracle_prop.nlin_out.est
      
    )
  )
  
}

# Experiment meta data
R.version <- version
glmnet.version <- packageVersion("glmnet")
system.version <- Sys.info()
sim.time <- Sys.time()
sim.date <- Sys.Date()

# Output file
filename.tracking = paste("../sim_tracking/standard_estimators",sim.time,sep="")

# Set parameters
N.replicates=1000
delta = 100/7 # value used in Figure 1 of Jiang, Mukherjee, Sen, Sur
theta.d.0 = 0
theta.y.0 = 0
epsilon = .1 # overlap parameter

# experiment specs
nsmall = 100
nbig = 10000
ns.length = 10
ns <- round( exp( seq( log(nsmall) , log(nbig) , length.out = ns.length) ) )
specs <- mapply(c,
                rep(1:N.replicates,times=ns.length),
                rep(ns,each=N.replicates),
                SIMPLIFY=FALSE)

# Run experiment
experiment.data <-
  mclapply(
    specs,
    function(spec){
      run.replicate(
        spec[2], # n
        round(spec[2]/delta), # p
        theta.d.0,
        theta.y.0,
        theta.d = c(1,rep(0,round(spec[2]/delta)-1)),
        theta.y = c(1,rep(0,round(spec[2]/delta)-1)),
        epsilon,
        spec[1],
        filename.tracking
      )
    },
    mc.cores=1 # Change to number of cores
  )
experiment.df <- do.call(rbind,experiment.data)
save(experiment.df,R.version,glmnet.version,system.version,sim.time,
     file=paste("../data/standard_estimators",sim.date,".Rda",sep=""))
write.csv(
  experiment.df,
  paste("../data/standard_estimators",sim.date,".csv",sep="")
)

