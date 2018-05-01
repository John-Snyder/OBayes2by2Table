library(hypergeo)
library(fmsb)
library(Rmpfr);PrecBits<-237
library(HI)
library(utils)

H.theta <- function(theta)
{
  #if(theta>1e15){return(log(theta)/sqrt(theta))}
  ifelse(theta<1,
         (pi*Re(hypergeo(1/2,1/2,1,1-theta,maxiter=1e6))),
         (pi/sqrt(theta)*Re(hypergeo(1/2,1/2,1,(theta-1)/theta,
                                     tol=1e-10 , maxiter=1e6))))
}

H.theta.int <- function(theta)
{
  integrate(f = function(t){
    1/(sqrt(t*(1-t))*sqrt(theta*(1-t)+t))
  },lower = 0,upper = 1,subdivisions = 1e8,stop.on.error=FALSE)$value
}

G.theta <- function(theta)
{
  integrate(f = function(t){
    log(t*(1-t))/(sqrt(t*(1-t))*sqrt(theta*(1-t)+t))
  },lower = 0,upper = 1,subdivisions = 1e6,stop.on.error=FALSE)$value
}

pi.R.theta=function(theta)
{
  (1/theta)*exp((1/2)*G.theta(theta)/H.theta(theta))
}

GetConnectionRatios <- function(approx.point.low=.01 - .00165165165,
                                approx.point.high=200 + .16546516546,
                                n1=n1,n2=n2,n3=n3,n4=n4)
{
  if(approx.point.low>0){
    ConnFact1 <- integrate(function(x,n1,n2,n3,n4){
      x^(n1+n2-1/2) * (1-x)^(n3+n4-1/2) * (1-((approx.point.low-1)/approx.point.low)*x)^(-n2-n4-1)},
      lower = 0,upper = 1,subdivisions = 1e8,n1=n1,n2=n2,n3=n3,n4=n4,abs.tol = 0L)$value/
      (beta(n1+n2+.5,n3+n4+.5)*approx.point.low^(n2+2)*H.theta(approx.point.low)) *
      exp((1/2)*G.theta(approx.point.low)/H.theta(approx.point.low))

    ConnFact2 <- ifelse(n1>n4,
                        -approx.point.low^(n4-3/4)/log(approx.point.low),
                        -approx.point.low^(n1-5/4)/log(approx.point.low))

    ConnectRatio.Low <- ConnFact1/ConnFact2
  }else{ConnectRatio.Low <- 1}


  ConnFact3 <- integrate(function(x,n1,n2,n3,n4){
    x^(n1+n2-1/2) * (1-x)^(n3+n4-1/2) * (1-((approx.point.high-1)/approx.point.high)*x)^(-n2-n4-1)},
    lower = 0,upper = 1,subdivisions = 1e8,n1=n1,n2=n2,n3=n3,n4=n4,abs.tol = 0L)$value/
    (beta(n1+n2+.5,n3+n4+.5)*approx.point.high^(n2+2)*H.theta(approx.point.high)) *
    exp((1/2)*G.theta(approx.point.high)/H.theta(approx.point.high))

  ConnFact4 <- ifelse(n3>n2,
                      1/(approx.point.high^(n2+7/4)*log(approx.point.high)),
                      1/(approx.point.high^(n3+5/4)*log(approx.point.high)))


  ConnectRatio.High <- ConnFact3/ConnFact4
  return(c(approx.point.low,ConnectRatio.Low,approx.point.high,ConnectRatio.High))
}

Posterior_Theta.n_o2 <- function(theta,n1,n2,n3,n4,
                                 approx.point.low,ConnectRatio.Low,
                                 approx.point.high,ConnectRatio.High)
{
  if(theta<approx.point.low)
  {
    ifelse(n1>n4,
           -ConnectRatio.Low * theta^(n4-3/4)/log(theta),
           -ConnectRatio.Low * theta^(n1-5/4)/log(theta))
  }else if(theta>approx.point.low & theta < approx.point.high){
    (integrate(function(x,n1,n2,n3,n4){
      x^(n1+n2-1/2) * (1-x)^(n3+n4-1/2) * (1-((theta-1)/theta)*x)^(-n2-n4-1)},
      lower = 0,upper = 1,subdivisions = 1e8,n1=n1,n2=n2,n3=n3,n4=n4,abs.tol = 0L)$value/
       (beta(n1+n2+.5,n3+n4+.5)*theta^(n2+2)*H.theta(theta))) *
      exp((1/2)*G.theta(theta)/H.theta(theta))
  }else if(theta>approx.point.high){
    ifelse(n3>n2,
           ConnectRatio.High/(theta^(n2+7/4)*log(theta)),
           ConnectRatio.High/(theta^(n3+5/4)*log(theta)))
  }
}

Post.Theta.Vec <- function(thetaS,n1,n2,n3,n4,
                           approx.point.low,ConnectRatio.Low,
                           approx.point.high,ConnectRatio.High)
{
  thetaS.Matrix <- matrix(thetaS,ncol=1)
  apply(thetaS.Matrix,1,Posterior_Theta.n_o2,n1=n1,n2=n2,n3=n3,n4=n4,
        approx.point.low=approx.point.low,ConnectRatio.Low=ConnectRatio.Low,
        approx.point.high=approx.point.high,ConnectRatio.High=ConnectRatio.High)
}

Percentile.Theta <- function(q,n1,n2,n3,n4,
                             approx.point.low,ConnectRatio.Low,
                             approx.point.high,ConnectRatio.High,
                             NormalizingConstant)
{
  integrate(Post.Theta.Vec,lower = 0,upper = q,n1=n1,n2=n2,n3=n3,n4=n4,
            approx.point.low=approx.point.low,ConnectRatio.Low=ConnectRatio.Low,
            approx.point.high=approx.point.high,ConnectRatio.High=ConnectRatio.High,abs.tol = 0L)$value/NormalizingConstant
}

QuantileDiff <- function(q,targetP,n1,n2,n3,n4,
                         approx.point.low,ConnectRatio.Low,
                         approx.point.high,ConnectRatio.High,
                         NormalizingConstant)
{
  Percentile.Theta(q,n1=n1,n2=n2,n3=n3,n4=n4,
                   approx.point.low=approx.point.low,ConnectRatio.Low=ConnectRatio.Low,
                   approx.point.high=approx.point.high,ConnectRatio.High=ConnectRatio.High,
                   NormalizingConstant) - targetP
}

eta1.given.theta<-function(eta1,theta,n1,n2,n3,n4)
{
  eta1^(n1+n2-1/2)*(1-eta1)^(n2+n4-1/2)*(eta1 + theta*(1-eta1))^(-(n2+n4+1))
}

POST.eta1theta=function(theta,eta1,n1,n2,n3,n4)
{
  ((theta^(n4-1))/H.theta(theta))*exp((1/2)*G.theta(theta)/H.theta(theta))*
    ((eta1^(n1+n2-1/2))*((1-eta1)^(n3+n4-1/2)))/((theta*(1-eta1) + eta1)^(n2+n4+1))
}

POST.eta1_given_theta=function(eta1,theta,n1,n2,n3,n4)
{
  ((eta1^(n1+n2-1/2))*((1-eta1)^(n3+n4-1/2)))/((theta*(1-eta1) + eta1)^(n2+n4+1))
}

POST.theta_given_eta1=function(theta,eta1,n1,n2,n3,n4)
{
  ((theta^(n4-1))/H.theta(theta))*exp((1/2)*G.theta(theta)/H.theta(theta))/
    ((theta*(1-eta1) + eta1)^(n2+n4+1))
}



Theta_ROU <-function(n1,n2,n3,n4,nsim=100,plots=0,ConnInfo)
{
  print(sprintf("Generating %s posterior samples from theta",nsim))
  ThetaSamples<-rep(NA,nsim)

    a.inf.theta<-0
  a.sup.theta<-optimize(function(th){
    sqrt(Posterior_Theta.n_o2(th,n1=n1,n2=n2,n3=n3,n4=n4,
                         approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                         approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4]))
  },
  interval = c(0,1000),
  maximum = TRUE)$objective

  b.inf.theta<-0
  b.sup.theta<-optimize(function(th){
    th*sqrt(Posterior_Theta.n_o2(th,n1=n1,n2=n2,n3=n3,n4=n4,
                              approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                              approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4]))
  },
  interval = c(0,1000),
  maximum = TRUE)$objective

  #sample theta
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  c<-1
  while(is.na(ThetaSamples[nsim]))
  {
    u<-runif(n = 1,min = a.inf.theta,max = a.sup.theta)
    v<-runif(n = 1,min = b.inf.theta,max = b.sup.theta)

    if(v/u > 1e-15 & v/u < 1e8)
    {
      Value<-NA
      try(Value<-sqrt(Posterior_Theta.n_o2(v/u,n1=n1,n2=n2,n3=n3,n4=n4,
                                           approx.point.low=ConnInfo[1],ConnectRatio.Low=ConnInfo[2],
                                           approx.point.high=ConnInfo[3],ConnectRatio.High=ConnInfo[4])),silent = TRUE)
      if(is.na(Value)){
        Value<-sqrt(Posterior_Theta.n(v/u,n1=n1,n2=n2,n3=n3,n4=n4))}

      if(u <= Value){
        ThetaSamples[c]<-v/u
        c<-c+1
        setTxtProgressBar(pb, c)
        #print(c)
      }
    }
  }
  return(ThetaSamples)
}

Eta1_ROU <- function(theta.sim,n1,n2,n3,n4)
{
  nsim <- length(theta.sim)
  eta1.S <- rep(NA,nsim)
  print(sprintf("Generating %s posterior samples from eta1",nsim))
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for(i in 1:nsim)
  {
    #print(i)
    #====== Sim eta1 | theta ======#
    theta.cur <- theta.sim[i]

    start <- Sys.time()
    a.inf.theta<-0
    a.sup.theta<-optimize(function(x) sqrt(POST.eta1_given_theta(x,theta=theta.cur,n1=n1,n2=n2,n3=n3,n4=n4)),
                          interval = c(0,1),
                          maximum = TRUE)$objective
    b.inf.theta<-0
    b.sup.theta<-optimize(function(x) x*sqrt(POST.eta1_given_theta(x,theta=theta.cur,n1=n1,n2=n2,n3=n3,n4=n4)),
                          interval = c(0,1),
                          maximum = TRUE)$objective

    repeat{
      u<-runif(n = 1,min = a.inf.theta,max = a.sup.theta)
      v<-runif(n = 1,min = b.inf.theta,max = b.sup.theta)

      if(v/u<1){
        if(u <= sqrt(POST.eta1_given_theta(eta1=v/u,theta=theta.cur,n1=n1,n2=n2,n3=n3,n4=n4))){
          eta1.S[i]<-v/u
          #c<-c+1
          #setTxtProgressBar(pb, c)
          #print(c)
          break
        }
      }
    }
    #print(sprintf("iteration %s: eta1 opt time=%s, eta1 ROU time=%s",
    #              i,eta1.opt.time,eta1.ROU.time))
    setTxtProgressBar(pb, i)
  }

  return(eta1.S)
  close(pb)
}

Eta3_sim <- function(theta.sim,eta1.sim,n1,n2,n3,n4)
{

  nsim <- length(theta.sim)
  eta3.S <- rep(NA,nsim)
  print(sprintf("Generating %s posterior samples from eta3",nsim))
  for(i in 1:nsim)
  {
    gamma.star <- rbeta(n = 1,shape1 = n2+n4+.5,shape2 = n1+n3)
    eta3.S[i] <- (1/(eta1.sim[i]+theta.sim[i]*(1-eta1.sim[i])))*
      ((gamma.star)/(1-gamma.star))
  }
  return(eta3.S)
}


Eta2_sim <- function(theta.sim,eta1.sim,eta3.sim,n1,n2,n3,n4)
{

  nsim <- length(theta.sim)
  eta2.S <- rep(NA,nsim)
  print(sprintf("Generating %s posterior samples from eta2",nsim))
  for(i in 1:nsim)
  {
    eta2.S[i] <- rgamma(n = 1,
                           shape = n1+n2+n3+n4+.5,
                           rate = eta1.sim[i]*(1+eta3.sim[i]) +
                             (1-eta1.sim[i])*(1+theta.sim[i]*eta3.sim[i])
    )
  }
  return(eta2.S)
}







