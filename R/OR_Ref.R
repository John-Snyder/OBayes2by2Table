#' Reference Analysis of the 2x2 table
#'
#' Perform Bayesian analysis on the 2x2 table using a reference prior
#'
#' @param n1 Upper left element of the table
#' @param n2 Upper right element of the table
#' @param n3 Lower left element of the table
#' @param n4 Lower right element of the table
#' @param conf.int Should a credible interval be computed using numerical integration?
#' @param conf.level confidence level for the returned credible interval. NOTE: Equal tailed
#' @param post.sample Should Posterior sampling be conducted?
#' @param sampling.depth Depth of sampling.  1 for theta only, 2 to add eta1, 3 to add eta3, 4 to add eta2
#' @param num.samples Number of Posterior Samples
#' @export
#' @examples
#' n1=5
#' n2=0
#' n3=8
#' n4=15
#'
#' res <- OR_Ref(n1,n2,n3,n4,conf.int = FALSE,num.samples = 500)
OR_Ref=function(n1,n2,n3,n4,conf.int = TRUE,conf.level = 0.95,
                post.sample = TRUE,sampling.depth=4,num.samples=1000){

  Final.Results <- NULL
  Lam1.post <- rgamma(n=5000,shape = n1+1/2, scale = 1)
  Lam2.post <- rgamma(n=5000,shape = n2+1/2, scale = 1)
  Lam3.post <- rgamma(n=5000,shape = n3+1/2, scale = 1)
  Lam4.post <- rgamma(n=5000,shape = n4+1/2, scale = 1)

  JeffThetaSamp <- Lam1.post*Lam4.post/(Lam2.post*Lam3.post)
  Quantiles <- quantile(JeffThetaSamp,c(.00001,(1-conf.level)/2,.5,1-(1-conf.level)/2,1-.00001))

  ConnInfo <- GetConnectionRatios(approx.point.low=Quantiles[1]/10,
                                  approx.point.high=min(Quantiles[5]*10,25000),
                                  n1=n1,n2=n2,n3=n3,n4=n4)
  NormConst <- integrate(Post.Theta.Vec,lower = 0,upper = Inf,n1=n1,n2=n2,n3=n3,n4=n4,
                         approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                         approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4],abs.tol = 0L)$value

  if(conf.int)
  {
    cat("Computing CI via numerical integration","\n")

   LB <- uniroot(QuantileDiff,
                  lower=Quantiles[1],
                  upper=Quantiles[3],
                  targetP=(1-conf.level)/2,n1=n1,n2=n2,n3=n3,n4=n4,
                  approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                  approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4],
                  NormalizingConstant=NormConst)$root

    MED <- uniroot(QuantileDiff,
                   lower=LB,
                   upper=Quantiles[4],
                   targetP=.5,n1=n1,n2=n2,n3=n3,n4=n4,
                   approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                   approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4],
                   NormalizingConstant=NormConst)$root

    UB <- uniroot(QuantileDiff,
                  lower=MED,
                  upper=ifelse(sum(c(n1,n2,n3,n4)==0)>1,150000,Quantiles[4]),
                  targetP=1-(1-conf.level)/2,n1=n1,n2=n2,n3=n3,n4=n4,
                  approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                  approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4],
                  NormalizingConstant=NormConst)$root

    Final.Results$CI <- c(LB,MED,UB)
  }

  if(post.sample)
  {
    Final.Results$Theta.Samples <-Theta_ROU(n1,n2,n3,n4,nsim=num.samples,ConnInfo=ConnInfo)

    if(sampling.depth>=2) Final.Results$Eta1.Samples <- Eta1_ROU(Final.Results$Theta.Samples,n1,n2,n3,n4)

    if(sampling.depth>=3) Final.Results$Eta3.Samples <- Eta3_sim(theta.sim = Final.Results$Theta.Samples,
                                                                 eta1.sim = Final.Results$Eta1.Samples,n1,n2,n3,n4)

    if(sampling.depth>=4) Final.Results$Eta2.Samples <- Eta2_sim(theta.sim = Final.Results$Theta.Samples,
                                                                 eta1.sim = Final.Results$Eta1.Samples,
                                                                 eta3.sim = Final.Results$Eta3.Samples,n1,n2,n3,n4)
  }
  return(Final.Results)
}





