###############################################################################
## Statistical Inference for Hicks–Moorsteen Productivity Indices
## Author: Shirong Zhao
## The programming codes used in this paper involve 
## some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
###############################################################################
rm(list=ls())
require(readxl)
require(FEAR)
source("./Functions/coverage.simple.mean.R")
source("./Functions/coverage.agg.R")
#### Set Seed ####
if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(90001)
################################################
np=2
nq=1
######################################################
pwt100 <- read_excel("Data/pwt100.xlsx", sheet = "Data")
######################################################
bhz<-c("Albania","Argentina","Armenia","Australia", "Austria", "Azerbaijan","Belarus",
       "Belgium", "Bolivia (Plurinational State of)","Brazil","Bulgaria","Canada","Chile",
       "China","Colombia","Costa Rica","Croatia","Czech Republic","Denmark","Dominican Republic",
       "Ecuador","Estonia","Finland","France","Germany","Greece","Guatemala","Honduras","China, Hong Kong SAR",
       "Hungary","Iceland","India","Indonesia" ,"Ireland","Israel","Italy","Jamaica","Japan",
       "Kazakhstan","Kenya","Republic of Korea","Kyrgyzstan","Latvia", "Lithuania","North Macedonia",
       "Madagascar","Malawi","Malaysia",
       "Mauritius","Mexico","Republic of Moldova","Morocco","Netherlands","New Zealand","Nigeria","Norway","Panama",
       "Paraguay","Peru","Philippines","Poland","Portugal","Romania","Russian Federation","Sierra Leone" ,
       "Singapore","Slovakia","Slovenia","Spain","Sri Lanka","Sweden","Switzerland","Syrian Arab Republic",
       "Taiwan", "Tajikistan","Thailand","Turkey","Ukraine","United Kingdom","Uruguay","United States","Venezuela (Bolivarian Republic of)",
       "Zambia","Zimbabwe")

developed<-c("Australia","Austria","Belgium","Canada","Denmark","Finland","France","Germany","Greece","China, Hong Kong SAR",
             "Iceland","Ireland","Israel","Italy","Japan","Republic of Korea","Netherlands","New Zealand","Norway","Portugal","Singapore","Spain","Sweden",
             "Switzerland","Taiwan","United Kingdom","United States")


year=c(1990,1995,2000,2005,2010,2015,2019)

ii=which((pwt100$country %in% bhz) & (pwt100$year %in% year))
df<-pwt100[ii,]


ii=which(df$country %in% developed)
dfa<-df[ii,]
dfb<-df[-ii,]


####### developed and developing countries ########

year=c(1990,1995,2000,2005,2010,2015,2019)
res.agg=matrix(NA,nrow=length(year),ncol=8)
res.simple=matrix(NA,nrow=length(year),ncol=8)


for (i in 1:length(year)) {
  
  cat (i,"\n")
  
  if (i<length(year)){
    
    i1=which(df$year==year[i])
    df1=df[i1,]
    x1=cbind(df1$emp,df1$cn)
    y1=matrix(df1$rgdpo,ncol=1)
    
    i2=which(df$year==year[i+1])
    df2=df[i2,]
    x2=cbind(df2$emp,df2$cn)
    y2=matrix(df2$rgdpo,ncol=1)
    
    
  } else {
    
    i1=which(df$year==year[1])
    df1=df[i1,]
    x1=cbind(df1$emp,df1$cn)
    y1=matrix(df1$rgdpo,ncol=1)
    i2=which(df$year==year[i])
    df2=df[i2,]
    x2=cbind(df2$emp,df2$cn)
    y2=matrix(df2$rgdpo,ncol=1)
    
  }
  
  tt=coverage.agg(x1=x1,y1=y1,x2=x2,y2=y2)
  
  res.agg[i,1:2]=exp(tt$estimate[c(1,3)])
  res.agg[i,3:4]=exp(tt$bounds1[1,])
  res.agg[i,5:6]=exp(tt$bounds1[2,])
  res.agg[i,7:8]=exp(tt$bounds1[3,])
  
  ###############
  ss=coverage.simple.mean(x1=x1,y1=y1,x2=x2,y2=y2)
  
  res.simple[i,1:2]=exp(ss$estimate[c(1,3)])
  res.simple[i,3:4]=exp(ss$bounds1[1,])
  res.simple[i,5:6]=exp(ss$bounds1[2,])
  res.simple[i,7:8]=exp(ss$bounds1[3,])
  
}


round(res.simple,4)
round(res.agg,4)


#### Latex
year1=c(1990,1995,2000,2005,2010,2015,1990)
year2=c(1995,2000,2005,2010,2015,2019,2019)

tex1=formatC(year1,width=6,digits=0,format="f")
tex2=formatC(year2,width=6,digits=0,format="f")



tex=paste(tex1,"&",tex2)
for (k in 1:ncol(res.simple)) {
  tex = paste(tex,"&",formatC(res.simple[,k],width=6,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")

write(tex,file="./Output/coverage-pwt-simple-all.tex")



tex=paste(tex1,"&",tex2)
for (k in 1:ncol(res.agg)) {
  tex = paste(tex,"&",formatC(res.agg[,k],width=6,digits = 4,format = "f"))
}
tex = paste(tex,"\\\\")

write(tex,file="./Output/coverage-pwt-agg-all.tex")

