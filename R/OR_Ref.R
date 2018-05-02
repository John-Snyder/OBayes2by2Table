#' @title 
#'  Reference Analysis of the 2x2 table
#'
#'@description 
#' Performs Objective Bayesian analysis on the odds ratio \eqn{\theta} of a 2x2 contingency 
#' table using the reference prior of Snyder and Sun 2018.  
#'
#'
#' @param x a 2x2 contingency table in matrix form
#' @param conf.int Should a credible interval be computed using numerical integration?
#' @param conf.level confidence level for the returned credible interval. NOTE: Equal tailed
#' @param post.sample Should Posterior sampling be conducted?
#' @param sampling.depth Depth of sampling.  1 for theta only, 2 to add eta1, 3 to add eta3, 4 to add eta2
#' @param num.samples Number of Posterior Samples
#' 
#' @details This methdology takes an unstructured  2x2 table of the form
#'  \tabular{rr}{\eqn{n1} \tab \eqn{n2} \cr \eqn{n3} \tab \eqn{n4}} where each 
#' \eqn{ni~Poisson(\lambdai)}.  From this, we primarilly seek to perform inference on the odds ratio 
#' \eqn{\theta = \lambda1\lambda4/\lambda2\lambda3}.
#' 
#' The Bayesian reference structure also contains several other nuisance parameters that may
#' be of interest to the researcher. We have
#' 
#' \deqn{\theta = \frac{\lambda1\lambda4}{\lambda2\lambda3},
#'       \eta1 = \frac{\lambda1}{\lambda1+\lambda3},
#'       \eta2 = \lambda1+\lambda3,
#'       \eta3 = \frac{\lambda2}{\lambda1}}{%
#'       \theta = \lambda1\lambda4/\lambda2\lambda3,
#'       \eta1 = \lambda1/(\lambda1+\lambda3),
#'       \eta2 = \lambda1+\lambda3,
#'       \eta3 = \lambda2/\lambda1}
#' 
#' We first note that despite being ``nuisance" parameters, the three \eqn{\etai} have
#'  highly informative interpretations.  \eqn{\eta_1} represents the proportion of 
#'  negative/positive outcomes in the control/treatment group.  \eqn{\eta_2} represents
#'  the total average number of negative/positive outcomes, and is a parameter that is
#'  related to the "scale" of the table.  Finally, \eqn{\eta_3} represents the relative
#'  risk of a positive outcome in the control group. 
#' @return 
#' \item{CI}{The Credible interval and median}
#' \item{Theta.Samples}{If \code{post.sample = TRUE}, the posterior Samples of the odds ratio theta}
#' \item{Eta1.Samples}{If \code{post.sample = TRUE} and \code{sampling.depth} is at least 2, 
#'                      the posterior Samples of \eqn{\eta1}}
#' \item{Eta2.Samples}{If \code{post.sample = TRUE} and \code{sampling.depth} is 4, 
#'                      the posterior Samples of \eqn{\eta2}}
#' \item{Eta3.Samples}{If \code{post.sample = TRUE} and \code{sampling.depth} is at least 3, 
#'                      the posterior Samples of \eqn{\eta3}}
#' @export
#' @examples
#' 
#' # Snyder and Sun (2018), Data were collected to establish a relationship between 
#' #  mobile phone and laptop operating system types for an undergraduate statistics class.
#' MyTable <- rbind(c(5,0),
#'                  c(8,15))
#' dimnames(MyTable) <-  list(Phone = c("Android", "Iphone"),Computer = c("Windows", "Mac"))
#' MyTable
#' 
#' # Fisher's exact test yields an infinitely wide interval and no sample estimate                 
#' fisher.test(MyTable)
#' 
#' # Results from a reference Bayesian Analysis are sensible.
#' res <- OR_Ref(MyTable,conf.int = TRUE,num.samples = 1000)
#' res$CI
#' quantile(res$Theta.Samples,c(.025,.5,.975))
#' 
#' # Agresti (1990, p. 61f; 2002, p. 91) Fisher's Tea Tasting Experiment
#' # A colleague of Fisher claimed to be able to distinguish if milk or
#' #  tea was added to a cup first.  She was given 8 cups of
#' #  tea, in which she was told that four of which had milk added first.  
#' #  The null hypothesis is that there is no association between the true 
#' #  order of pouring and the woman's guess, the alternative that there 
#' #  is a positive association (that the odds ratio is greater than 1).
#' 
#' TeaTasting <-
#' matrix(c(3, 1, 1, 3),
#'      nrow = 2,
#'      dimnames = list(Guess = c("Milk", "Tea"),
#'      Truth = c("Milk", "Tea")))
#' 
#' fisher.test(TeaTasting, alternative = "greater")
#' # => p = 0.2429, association could not be established
#' 
#' res <- OR_Ref(TeaTasting,conf.int = FALSE,num.samples = 2500)
#' sum(res$Theta.Samples>1)/sum(res$Theta.Samples<1)
#' # Bayes Factor indicates substantial evidence of her ability
#' 
#' # Fisher (1962, 1970), Criminal convictions of like-sex twins
#' Convictions <-
#'   matrix(c(2, 10, 15, 3),
#'          nrow = 2,
#'          dimnames =
#'            list(c("Dizygotic", "Monozygotic"),
#'                 c("Convicted", "Not convicted")))
#' Convictions
#' fisher.test(Convictions, alternative = "less")
#' #Fisher's exact test is infinitely wide on log scale.
#' 
#' res <- OR_Ref(Convictions,conf.int = FALSE,num.samples = 5000)
#' quantile(res$Theta.Samples,c(.025,.5,.975))
#' 
OR_Ref=function(x,conf.int = TRUE,
                conf.level = 0.95,
                post.sample = TRUE,
                sampling.depth=1,
                num.samples=1000){

  n1 <- x[1,1]
  n2 <- x[1,2]
  n3 <- x[2,1]
  n4 <- x[2,2]
  
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
    cat("Computing CI via numerical integration...","\n")

   LB <- uniroot(QuantileDiff,
                  lower=Quantiles[1],
                  upper=Quantiles[3],
                  targetP=(1-conf.level)/2,n1=n1,n2=n2,n3=n3,n4=n4,
                  approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                  approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4],
                  NormalizingConstant=NormConst)$root
   cat("2.5% quantile finished","\n")
   if(1<0){
     MED <- uniroot(QuantileDiff,
                   lower=Quantiles[2],
                   upper=Quantiles[4],
                   targetP=.5,n1=n1,n2=n2,n3=n3,n4=n4,
                   approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                   approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4],
                   NormalizingConstant=NormConst)$root
   }

    UB <- uniroot(QuantileDiff,
                  lower=Quantiles[2],
                  upper=ifelse(sum(c(n1,n2,n3,n4)==0)>1,150000,Quantiles[4]),
                  targetP=1-(1-conf.level)/2,n1=n1,n2=n2,n3=n3,n4=n4,
                  approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                  approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4],
                  NormalizingConstant=NormConst)$root
    cat("97.5% quantile finished","\n")
    Final.Results$CI <- c(LB,UB)
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





