###############################################################################
## Statistical Inference for Hicks–Moorsteen Productivity Indices
## Author: Shirong Zhao
## The programming codes used in this paper involve 
## some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
###############################################################################
coverage.simple.mean <- function(x1,y1,x2,y2) {
  L=100
  np=ncol(x1)
  nq=ncol(y1)
  n=nrow(x1)
  na=floor(n/2)
  nb=n-na
  kappa=2/(np+nq+1)
  bc.fac=1/(2**kappa - 1)
  nk=floor(n**(2*kappa))
  # evaluate efficiency using VRS-DEA
  lam11= FEAR::dea(XOBS=t(x1),YOBS=t(y1),XREF=t(x1),YREF=t(y1), 
                   METRIC=2,ORIENTATION=2, RTS=1)
  lam12= FEAR::dea(XOBS=t(x2),YOBS=t(y1),XREF=t(x2),YREF=t(y2), 
                   METRIC=2,ORIENTATION=2, RTS=1)
  lam21= FEAR::dea(XOBS=t(x1),YOBS=t(y2),XREF=t(x1),YREF=t(y1), 
                   METRIC=2,ORIENTATION=2, RTS=1)
  lam22= FEAR::dea(XOBS=t(x2),YOBS=t(y2),XREF=t(x2),YREF=t(y2), 
                   METRIC=2,ORIENTATION=2, RTS=1)
  the11= FEAR::dea(XOBS=t(x1),YOBS=t(y1),XREF=t(x1),YREF=t(y1), 
                   METRIC=2,ORIENTATION=1, RTS=1)
  the12= FEAR::dea(XOBS=t(x1),YOBS=t(y2),XREF=t(x2),YREF=t(y2), 
                   METRIC=2,ORIENTATION=1, RTS=1)
  the21= FEAR::dea(XOBS=t(x2),YOBS=t(y1),XREF=t(x1),YREF=t(y1), 
                   METRIC=2,ORIENTATION=1, RTS=1)
  the22= FEAR::dea(XOBS=t(x2),YOBS=t(y2),XREF=t(x2),YREF=t(y2), 
                   METRIC=2,ORIENTATION=1, RTS=1)
  #
  HM1i=-0.5*(log(lam21)-log(lam11)-log(the21)+log(the11))
  HM2i=-0.5*(log(lam22)-log(lam12)-log(the22)+log(the12))
  HMi=HM1i+HM2i
  HM=exp(mean(HMi))
  sig=HM*sd(HMi)
  # compute bias corrections via generalized jackknife:
  # involving some earlier codes from Paul Wilson
  tbar=0
  ind=c(1:n)
  for (j in 1:L) {
    if (j==1) {
      ind1=c(1:n)
      x1.b=x1
      y1.b=y1
      x2.b=x2
      y2.b=y2
    } else {
      ind1=sample(ind,size=n)
      x1.b[1:n,]=x1[ind1,]
      y1.b[1:n,]=y1[ind1,]
      x2.b[1:n,]=x2[ind1,]
      y2.b[1:n,]=y2[ind1,]
    }
    #
    x1a=matrix(x1.b[1:na,],ncol=np)
    y1a=matrix(y1.b[1:na,],ncol=nq)
    x1b=matrix(x1.b[(na+1):n,],ncol=np)
    y1b=matrix(y1.b[(na+1):n,],ncol=nq)
    #
    x2a=matrix(x2.b[1:na,],ncol=np)
    y2a=matrix(y2.b[1:na,],ncol=nq)
    x2b=matrix(x2.b[(na+1):n,],ncol=np)
    y2b=matrix(y2.b[(na+1):n,],ncol=nq)
    #
    lam11a= FEAR::dea(XOBS=t(x1a),YOBS=t(y1a),XREF=t(x1a),YREF=t(y1a), 
                     METRIC=2,ORIENTATION=2, RTS=1)
    lam12a= FEAR::dea(XOBS=t(x2a),YOBS=t(y1a),XREF=t(x2a),YREF=t(y2a), 
                       METRIC=2,ORIENTATION=2, RTS=1)
    lam21a= FEAR::dea(XOBS=t(x1a),YOBS=t(y2a),XREF=t(x1a),YREF=t(y1a), 
                       METRIC=2,ORIENTATION=2, RTS=1)
    lam22a= FEAR::dea(XOBS=t(x2a),YOBS=t(y2a),XREF=t(x2a),YREF=t(y2a), 
                       METRIC=2,ORIENTATION=2, RTS=1)
    the11a= FEAR::dea(XOBS=t(x1a),YOBS=t(y1a),XREF=t(x1a),YREF=t(y1a), 
                       METRIC=2,ORIENTATION=1, RTS=1)
    the12a= FEAR::dea(XOBS=t(x1a),YOBS=t(y2a),XREF=t(x2a),YREF=t(y2a), 
                       METRIC=2,ORIENTATION=1, RTS=1)
    the21a= FEAR::dea(XOBS=t(x2a),YOBS=t(y1a),XREF=t(x1a),YREF=t(y1a), 
                       METRIC=2,ORIENTATION=1, RTS=1)
    the22a= FEAR::dea(XOBS=t(x2a),YOBS=t(y2a),XREF=t(x2a),YREF=t(y2a), 
                       METRIC=2,ORIENTATION=1, RTS=1)
    #
    HM1ia=-0.5*(log(lam21a)-log(lam11a)-log(the21a)+log(the11a))
    HM2ia=-0.5*(log(lam22a)-log(lam12a)-log(the22a)+log(the12a))
    HMia=HM1ia+HM2ia
    HMa=exp(mean(HMia))
    #
    lam11b= FEAR::dea(XOBS=t(x1b),YOBS=t(y1b),XREF=t(x1b),YREF=t(y1b), 
                     METRIC=2,ORIENTATION=2, RTS=1)
    lam12b= FEAR::dea(XOBS=t(x2b),YOBS=t(y1b),XREF=t(x2b),YREF=t(y2b), 
                     METRIC=2,ORIENTATION=2, RTS=1)
    lam21b= FEAR::dea(XOBS=t(x1b),YOBS=t(y2b),XREF=t(x1b),YREF=t(y1b), 
                     METRIC=2,ORIENTATION=2, RTS=1)
    lam22b= FEAR::dea(XOBS=t(x2b),YOBS=t(y2b),XREF=t(x2b),YREF=t(y2b), 
                     METRIC=2,ORIENTATION=2, RTS=1)
    the11b= FEAR::dea(XOBS=t(x1b),YOBS=t(y1b),XREF=t(x1b),YREF=t(y1b), 
                     METRIC=2,ORIENTATION=1, RTS=1)
    the12b= FEAR::dea(XOBS=t(x1b),YOBS=t(y2b),XREF=t(x2b),YREF=t(y2b), 
                     METRIC=2,ORIENTATION=1, RTS=1)
    the21b= FEAR::dea(XOBS=t(x2b),YOBS=t(y1b),XREF=t(x1b),YREF=t(y1b), 
                     METRIC=2,ORIENTATION=1, RTS=1)
    the22b= FEAR::dea(XOBS=t(x2b),YOBS=t(y2b),XREF=t(x2b),YREF=t(y2b), 
                     METRIC=2,ORIENTATION=1, RTS=1)
    #
    HM1ib=-0.5*(log(lam21b)-log(lam11b)-log(the21b)+log(the11b))
    HM2ib=-0.5*(log(lam22b)-log(lam12b)-log(the22b)+log(the12b))
    HMib=HM1ib+HM2ib
    HMb=exp(mean(HMib))
    #
    tbar=tbar+(HMa+HMb)/2-HM
    #
  }
  # tbar contains the bias
  HM.bias=(1/L)*bc.fac*tbar
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  if (np+nq<4) {
    ts=sig/sqrt(n)
    bounds1=matrix((HM-HM.bias-ts*crit),nrow=3,ncol=2)
    # make a list of results to return to calling routine and then quit:
    res=list(bounds1=bounds1,
             sig=c(sig),
             estimate=c(HM,HM.bias,HM-HM.bias)
             )
  } else {
    HM.nk=exp(mean(HMi[1:nk]))
    ts.nk=sig/sqrt(nk)
    bounds2=matrix((HM.nk-HM.bias-ts.nk*crit),nrow=3,ncol=2)
    res=list(bounds2=bounds2,
             sig=c(sig),
             estimate=c(HM.nk,HM.bias,HM.nk-HM.bias)
    )
  }
  return(res)
}