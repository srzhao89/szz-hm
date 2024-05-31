###############################################################################
## Statistical Inference for Hicks–Moorsteen Productivity Indices
## Author: Shirong Zhao
## The programming codes used in this paper involve 
## some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
###############################################################################
coverage.agg <- function(x1,y1,x2,y2) {
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
  rev1=as.vector(y1)
  rev2=as.vector(y2)
  cos1=as.vector(y1)
  cos2=as.vector(y2)
  #
  U1=lam21*rev2
  U2=lam11*rev1
  U3=the21*cos2
  U4=the11*cos1
  #
  U5=lam22*rev2
  U6=lam12*rev1
  U7=the22*cos2
  U8=the12*cos1
  #
  U9=rev2
  U10=rev1
  U11=cos2
  U12=cos1
  #
  mu1=mean(U1)
  mu2=mean(U2)
  mu3=mean(U3)
  mu4=mean(U4)
  mu5=mean(U5)
  mu6=mean(U6)
  mu7=mean(U7)
  mu8=mean(U8)
  mu9=mean(U9)
  mu10=mean(U10)
  mu11=mean(U11)
  mu12=mean(U12)
  #
  xi=exp(-0.5*(log(mu1)-log(mu2)-log(mu3)+log(mu4)+
           log(mu5)-log(mu6)-log(mu7)+log(mu8))+
           log(mu9)-log(mu10)-log(mu11)+log(mu12)
  )
  #
  Ti=cbind(U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,U12)
  Sigma=cov(Ti)
  g1 =-1/(2*mu1)
  g2 = 1/(2*mu2)
  g3 = 1/(2*mu3)
  g4 =-1/(2*mu4)
  g5 =-1/(2*mu5)
  g6 = 1/(2*mu6)
  g7 = 1/(2*mu7)
  g8 =-1/(2*mu8)
  g9 = 1/mu9
  g10=-1/mu10
  g11=-1/mu11
  g12= 1/mu12
  g=c(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12)
  V.Sigma=t(g)%*%Sigma%*%g
  sig=xi*as.vector(sqrt(V.Sigma))
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
    ###################
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
    rev1a=as.vector(y1a)
    rev2a=as.vector(y2a)
    cos1a=as.vector(y1a)
    cos2a=as.vector(y2a)
    #
    U1a=lam21a*rev2a
    U2a=lam11a*rev1a
    U3a=the21a*cos2a
    U4a=the11a*cos1a
    #
    U5a=lam22a*rev2a
    U6a=lam12a*rev1a
    U7a=the22a*cos2a
    U8a=the12a*cos1a
    #
    mu1a=mean(U1a)
    mu2a=mean(U2a)
    mu3a=mean(U3a)
    mu4a=mean(U4a)
    mu5a=mean(U5a)
    mu6a=mean(U6a)
    mu7a=mean(U7a)
    mu8a=mean(U8a)
    #
    xia=exp(-0.5*(log(mu1a)-log(mu2a)-log(mu3a)+log(mu4a)+
              log(mu5a)-log(mu6a)-log(mu7a)+log(mu8a))+
              log(mu9)-log(mu10)-log(mu11)+log(mu12)
    )
    ################
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
    rev1b=as.vector(y1b)
    rev2b=as.vector(y2b)
    cos1b=as.vector(y1b)
    cos2b=as.vector(y2b)
    #
    U1b=lam21b*rev2b
    U2b=lam11b*rev1b
    U3b=the21b*cos2b
    U4b=the11b*cos1b
    #
    U5b=lam22b*rev2b
    U6b=lam12b*rev1b
    U7b=the22b*cos2b
    U8b=the12b*cos1b
    #
    mu1b=mean(U1b)
    mu2b=mean(U2b)
    mu3b=mean(U3b)
    mu4b=mean(U4b)
    mu5b=mean(U5b)
    mu6b=mean(U6b)
    mu7b=mean(U7b)
    mu8b=mean(U8b)
    #
    xib=exp(-0.5*(log(mu1b)-log(mu2b)-log(mu3b)+log(mu4b)+
              log(mu5b)-log(mu6b)-log(mu7b)+log(mu8b))+
              log(mu9)-log(mu10)-log(mu11)+log(mu12)
    )
    #
    tbar=tbar+(xia+xib)/2-xi
    #
  }
  # tbar contains the bias
  xi.bias=(1/L)*bc.fac*tbar
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  if (np+nq<4) {
    ts=sig/sqrt(n)
    bounds1=matrix((xi-xi.bias-ts*crit),nrow=3,ncol=2)
    # make a list of results to return to calling routine and then quit:
    res=list(bounds1=bounds1,
             sig=c(sig),
             estimate=c(xi,xi.bias,xi-xi.bias)
             )
  } else {
    #
    U1.nk=U1[1:nk]
    U2.nk=U2[1:nk]
    U3.nk=U3[1:nk]
    U4.nk=U4[1:nk]
    #
    U5.nk=U5[1:nk]
    U6.nk=U6[1:nk]
    U7.nk=U7[1:nk]
    U8.nk=U8[1:nk]
    #
    U9.nk=U9[1:nk]
    U10.nk=U10[1:nk]
    U11.nk=U11[1:nk]
    U12.nk=U12[1:nk]
    #
    mu1.nk=mean(U1.nk)
    mu2.nk=mean(U2.nk)
    mu3.nk=mean(U3.nk)
    mu4.nk=mean(U4.nk)
    mu5.nk=mean(U5.nk)
    mu6.nk=mean(U6.nk)
    mu7.nk=mean(U7.nk)
    mu8.nk=mean(U8.nk)
    mu9.nk=mean(U9.nk)
    mu10.nk=mean(U10.nk)
    mu11.nk=mean(U11.nk)
    mu12.nk=mean(U12.nk)
    #
    xi.nk=-0.5*(log(mu1.nk)-log(mu2.nk)-log(mu3.nk)+log(mu4.nk)+
                log(mu5.nk)-log(mu6.nk)-log(mu7.nk)+log(mu8.nk))+
                log(mu9.nk)-log(mu10.nk)-log(mu11.nk)+log(mu12.nk)
    ts.nk=sig/sqrt(nk)
    bounds2=matrix((xi.nk-xi.bias-ts.nk*crit),nrow=3,ncol=2)
    res=list(bounds2=bounds2,
             sig=c(sig),
             estimate=c(xi.nk,xi.bias,xi.nk-xi.bias)
    )
  }
  return(res)
}