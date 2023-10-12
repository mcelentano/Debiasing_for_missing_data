g.computation <-
  function(X,
           out.models,g.folds,
           out.idxs,g.idxs){
    
    g.terms <- mapply(
      function(out.idx,g.idx) {
        out.model <- out.models[[out.idx]]
        g.fold <- g.folds[[g.idx]]
        return( out.model[1] + X[g.fold,] %*% out.model[-1] )
      },
      out.idxs,g.idxs,
      SIMPLIFY=FALSE
    ) # returns list of arrays
    g.terms <- do.call(c,g.terms) # combines arrays
    
    return(mean(g.terms))
  }

aipw <-
  function(y,d,X,epsilon,
           out.models,prop.models,aipw.folds,
           out.idxs,prop.idxs,aipw.idxs){
    
    aipw.terms <- mapply(
      function(out.idx,prop.idx,aipw.idx){
        out.model <- out.models[[out.idx]]
        prop.model <- prop.models[[prop.idx]]
        aipw.fold <- aipw.folds[[aipw.idx]]
        
        pi.hat =  epsilon + (1-epsilon) * inv.logit( prop.model[1] + X[aipw.fold,] %*% prop.model[-1] )
        aipw.terms = (d[aipw.fold] / pi.hat)*y[aipw.fold] + (1 - d[aipw.fold] / pi.hat)*( out.model[1] + X[aipw.fold,] %*% out.model[-1] )
      },
      out.idxs,prop.idxs,aipw.idxs,
      SIMPLIFY=FALSE
    ) # returns list of arrays
    aipw.terms <- do.call(c,aipw.terms) # combines arrays
    
    return(mean(aipw.terms))
    
  }

ipw <- function(y,d,X,epsilon,
                prop.models,ipw.folds,
                prop.idxs,ipw.idxs){
  
  ipw.terms <- mapply(
    function(prop.idx,ipw.idx){
      
      prop.model <- prop.models[[prop.idx]]
      ipw.fold <- ipw.folds[[ipw.idx]]
      
      pi.hat =  epsilon + (1-epsilon) * inv.logit( prop.model[1] + X[ipw.fold,] %*% prop.model[-1] )
      ipw.terms = (d[ipw.fold] / pi.hat)*y[ipw.fold]
    },
    prop.idxs,ipw.idxs,
    SIMPLIFY=FALSE
  ) # returns list of arrays
  ipw.terms <- do.call(c,ipw.terms) # combines arrays
  
  return(mean(ipw.terms))
  
}