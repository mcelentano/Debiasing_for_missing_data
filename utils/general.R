# Define functions we need
d.inv.logit <- function(eta) {
  return(
    exp(eta)/(1+exp(eta))^2
  )
}
dd.inv.logit <- function(eta) {
  return(
    exp(eta)*(exp(eta)-1)/(exp(eta)+1)^3
  )
}
prop.score <- function(eta,epsilon) {
  return(
    epsilon + (1-epsilon) * inv.logit(eta)
  )
}
d.prop.score <- function(eta,epsilon) {
  return(
    (1-epsilon) * d.inv.logit(eta)
  )
}
dd.prop.score <- function(eta,epsilon) {
  return(
    (1-epsilon) * dd.inv.logit(eta)
  )
}
properr <- function(eta.d,eta.d.hat,epsilon) {
  return(
    1 - prop.score(eta.d,epsilon) / prop.score(eta.d.hat,epsilon)
  )
}
d.properr.d.etad <- function(eta.d,eta.d.hat,epsilon) {
  return(
    -d.prop.score(eta.d,epsilon) / prop.score(eta.d.hat,epsilon)
  )
}
d.properr.d.etadhat <- function(eta.d,eta.d.hat,epsilon) {
  return(
    prop.score(eta.d,epsilon) * d.prop.score(eta.d.hat,epsilon) / prop.score(eta.d.hat,epsilon)^2
  )
}
dd.properr.dd.etad <- function(eta.d,eta.d.hat,epsilon) {
  return(
    -dd.prop.score(eta.d,epsilon) / prop.score(eta.d.hat,epsilon)
  )
}
dd.properr.d.etad.d.etadhat <- function(eta.d,eta.d.hat,epsilon) {
  return(
    d.prop.score(eta.d,epsilon)*d.prop.score(eta.d.hat,epsilon) / prop.score(eta.d.hat,epsilon)^2
  )
}
dd.properr.dd.etadhat <- function(eta.d,eta.d.hat,epsilon) {
  return(
    prop.score(eta.d,epsilon) * dd.prop.score(eta.d.hat,epsilon) / prop.score(eta.d.hat,epsilon)^2 -
      2 * prop.score(eta.d,epsilon) * d.prop.score(eta.d.hat,epsilon)^2 / prop.score(eta.d.hat,epsilon)^3
  )
}
sqrtm <- function(K) {
  SVD <- svd(K)
  return(
    SVD$u %*% diag(sqrt(SVD$d)) %*% t(SVD$v)
  )
}