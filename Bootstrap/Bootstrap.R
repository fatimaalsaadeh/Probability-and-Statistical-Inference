#central limit theorem upper and lower bounds function
#n: size of the population vector, 
#mean,sd: sd,mean of population vector, 
#return list of lower and upper bound 
CLT<-function(bootvec,mean0, sd0, alpha, n)  {
  CLTL<-mean0-((sd0/sqrt(n))*qt(1-(alpha/2),n-1))
  CLTU<-mean0-((sd0/sqrt(n))*qt((alpha/2),n-1))
  list(CLTL=CLTL, CLTU=CLTU)
}
#bootstrap function
#v1: population vector, 
#nboot: number of times for resampling, 
#alpha: to define the ratio of upper/lower bounds
#return bb: bootsrap bias vector, bsde: bootstrap sd error,
#bBCI: bootstrap confidence interval vector , 
#bNCI: normal distribution confidence interval
#bNCIBYN:normal distribution confidence interval with sd divided by n
Bootsrap<-function(v1, nboot=1000, alpha=0.1) {
  #extract the length mean and standard diviation from the population vector
  n<-length(v1)
  mean0<-mean(v1)
  sd0<-sqrt(var(v1))
  
  #initiate the boot and bias vectors
  bootvec<-NULL
  bootbiasvec<-NULL
  
  #resampling with replacement nboot times 
  for (i in 1:nboot) {
    #create the sample vector from the population vector with replacement 
    samplevec = sample(v1,replace=T)
    mean1=mean(samplevec)
    sd1=sqrt(var(samplevec))
    #create the boot vector by calculating the coefficient of variation  
    bootvec<-c(bootvec,(mean1-mean0)/(sd1/sqrt(n)))
    #create the bootstrap bias vector by calculating the mean difference
    #between original populatio vector and sample vector 
    bootbiasvec<-c(bootbiasvec,mean1-mean0)
  }
  #extract the  mean and standard diviation and standard error from the boot vector
  bootmean = mean(bootvec)
  bootsd= sqrt(var(bootvec))
  bootsde= sqrt(var(bootvec))/sqrt(n)
  #calculate the lower and upper pivotal bound for the bootstrap vector
  uq = quantile(bootvec,1-(alpha/2))
  LB = mean0-(sd0/sqrt(n))*uq
  
  lq = quantile(bootvec,alpha/2)
  UB = mean0-(sd0/sqrt(n))*lq
  
  #calculate the lower and upper bound for the normal distribution
  NLB=mean0-(sd0/sqrt(n))*qnorm(1-(alpha/2))
  NUB=mean0+(sd0/sqrt(n))*qnorm(1-(alpha/2))
  
  
  #calculate the lower and upper percentile bound
  PLB=mean0-quantile(bootbiasvec,1-(alpha/2))
  PUB=mean0-quantile(bootbiasvec,(alpha/2))
  
  #calculate the lower and upper bound for the normal distribution with sd divided by n
  NLBYN=mean0-((sd0/n)*qnorm(1-(alpha/2)))
  NUBYN=mean0+((sd0/n)*qnorm(1-(alpha/2)))

  #calculate the pivotal bootsrap confidence interval
  BCI<-c(LB, UB)
  
  #calculate the normal distribution confidence interval
  NBCI<-c(NLB, NUB)
  
  #calculate the percentile confidence interval
  PBCI<-c(PLB, PUB)
  #calculate the normal distribution with sd divided by n confidence interval
  NCIBYN<-c(NLBYN, NUBYN)
  #calculate the bootstrap bias
  biasmean=mean(bootbiasvec)
  #return list with the calculated parameters
  list(bb=biasmean, bsde=bootsde, bBCI=BCI, bNCI=NBCI, bPBCI=PBCI, bNCIBYN=NCIBYN)
}

#jackknife function
#v1: population vector, 
#alpha: to define the ratio of upper/lower bounds
#return jb: jackknife bias, jsde: jackknife sd error,
#jCI: jackknife confidence interval vector , 
#jNCI: normal distribution confidence interval
Jackknife<-function(v1, alpha=0.1, statfunc=sd)
{
  #extract the length, mean and standard diviation from the population vector
  n<-length(v1)
  mean0<-statfunc(v1)
  sd0<-sd(v1)
  
  #jackknife vector of the n sampled points
  jackvec<-NULL
  
  #computing the jack sample vector.
  for(i in 1:n){
    m1<-statfunc(v1[-i])
    jackvec<-c(jackvec, n*(mean0)-(n-1)*m1)
  }
  
  #extract the  mean and standard diviation and standard error from the boot vector
  jackvecmean = mean(jackvec)
  jackvecsd= sd(jackvec)
  jackvecsde= sqrt(var(jackvec))/sqrt(n)
  
  #jackknife bias 
  jackbias<-mean(jackvec)-mean0
  
  #calculate the lower and upper quantile for the jackknife vector
  uq = qt(1-(alpha/2),n-1)
  LB=mean0-(sd0/sqrt(n))*uq
  
  lq = qt(alpha/2, n-1)
  UB=mean0-(sd0/sqrt(n))*lq

  #calculate the lower and upper quantile for the normal distribution
  NLB=mean0-(sd0/sqrt(n))*qnorm(1-(alpha/2))
  NUB=mean0+(sd0/sqrt(n))*qnorm(1-(alpha/2))
  
  #calculate the jackknife confidence interval
  JCI<-c(LB, UB)
  
  #calculate the normal distribution confidence interval
  NCI<-c(NLB, NUB)
  
  #return list with the calculated parameters
  list(jb=jackbias, jsde=jackvecsde, jCI=JCI, jNCI=NCI)
}

Sim.func<-function(meanv=3, n=30, nsim=1000, alpha=0.1)
{
  #create coverage indicator vectors for bootstrap,normal,
  #central limit theorem,jackknife,normal with sd divided by n
  cbootvec<-NULL
  cbootpercentilevec<-NULL
  cnormvec<-NULL
  ccltvec<-NULL
  cjackvec<-NULL
  cnormbynvec<-NULL
  #calculate real mean
  mulnorm<-(exp(meanv+(1/2)))
  
  #run simulation with nsim times
  for(i in 1:nsim){
    #create a sample lognormal distribution vector
    samplevec<-rlnorm(n,meanv)
    
    #extract mean and sd from the sample vector
    samplemean=mean(samplevec)
    samplesd=sqrt(var(samplevec))
    #run bootsrap resampling
    bootlist<-Bootsrap(samplevec)
    
    #extract the bootstrap confidence interval from the returned bootstrap list
    bootconf<-bootlist$bBCI
    
    #extract the bootstrap percentil confidence interval from the returned bootstrap list
    cbootpercentileconf<-bootlist$bPBCI
    
    #extract the normal distribution confidence interval from the returned bootstrap list
    normconf<-bootlist$bNCI
    
    #extract the normal distribution divided by n confidence interval
    normconfbyn<-bootlist$bNCIBYN
    
    #get the central limit theorem confidence interval for the sample
    
    cltlist<-CLT(samplevec,samplemean,samplesd, alpha, length(samplevec))

    #run jackknife resampling
    jackknifelist<-Jackknife(samplevec, alpha)
    
    #extract the jackknife confidence interval from the returned jackknife list
    jackknifeconf<-jackknifelist$jCI
    
    #calculate the coverage vectors for bootsrap,normal dist.,normal by n, clt and jackknife 
    cbootvec<-c(cbootvec,(bootconf[1]<mulnorm)*(bootconf[2]>mulnorm))
    cnormvec<-c(cnormvec,(normconf[1]<mulnorm)*(normconf[2]>mulnorm))
    cnormbynvec<-c(cnormbynvec,(normconfbyn[1]<mulnorm)*(normconfbyn[2]>mulnorm))
    ccltvec<-c(ccltvec,(cltlist[1]<mulnorm)*(cltlist[2]>mulnorm))
    cjackvec<-c(cjackvec,(jackknifeconf[1]<mulnorm)*(jackknifeconf[2]>mulnorm))


    cbootpercentilevec<-c(cbootpercentilevec,(cbootpercentileconf[1]<mulnorm)*(cbootpercentileconf[2]>mulnorm))
  }
  
  hist(samplevec)
  #calculate and output coverage probability estimates
  list(bootcoverage=(sum(cbootvec)/nsim),
       bootpercentilcoverage=(sum(cbootpercentilevec)/nsim),
       normcoverage=(sum(cnormvec)/nsim),
       normbyncoverage=(sum(cnormbynvec)/nsim),
       cltcoverage=(sum(ccltvec)/nsim),
       jackcoverage=(sum(cjackvec)/nsim)
       )

}



vector1=c(60,75,80,85,90)
#print(Jackknife(vector1))
#print(Bootsrap(vector1))
#print(my.bootstrapci.ml(vector1))
#print(CLT(vector1,mean(vector1),sd(vector1),.05,length(vector1)))
#sampling on 10, 30, 100 population sample size with alpha=0,05
print(Sim.func(n=10, alpha=0.05))
print(Sim.func(alpha=.05))
print(Sim.func( alpha=0.05))


