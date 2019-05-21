#################################################################################
##
## Author:  Nat Goodman
## Created: 19-02-19
##          from dodata.R created 19-02-18
##          from ovrfx.R created 19-02-03 
##          from siglo.R 19-01-01
##          from repwr/R/repwr.R created 17-10-05 
##           and repwr/R/sim.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate data for misig documents
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Top Level Data Generation Function ----
dodat=function() {
  param(datfun);
  datfun();
}
## ---- Simulation Functions ----
vrnorm=Vectorize(rnorm,"mean");
## random-d simulation
dosim_rand=function(n,m,d) {
  param(verbose);
  sapply(n,function(n) {
    if (verbose) print(paste(sep=' ','>>> dosim_rand:',nvq(n)));
    group0=replicate(m,rnorm(n,mean=0));
    group1=vrnorm(n,mean=d);
    mean0=colMeans(group0);
    mean1=colMeans(group1);
    d.raw=mean1-mean0;
    sd0=apply(group0,2,sd);
    sd1=apply(group1,2,sd);
    sd=pooled_sd(sd0,sd1);
    d.sdz=d.raw/sd;
    sim=data.frame(n,d.pop=d,d.sdz,sd,d.raw,mean0,mean1,sd0,sd1,
                   row.names=NULL,stringsAsFactors=F);
    save_sim_rand(sim,n);
  });
  return(T);
}
## fixed-d simulation
dosim_fixd=function(n,m,d) {
  param(verbose);
  cases=expand.grid(n=n,d=d);
  if (nrow(cases)>0)
    apply(cases,1,function(row) {
      n=row[1]; d=row[2];
      if (verbose) print(paste(sep=' ','>>> dosim_fixd:',nvq(n,d)));
      group0=replicate(m,rnorm(n,mean=0));
      group1=replicate(m,rnorm(n,mean=d));
      mean0=colMeans(group0);
      mean1=colMeans(group1);
      d.raw=mean1-mean0;
      sd0=apply(group0,2,sd);
      sd1=apply(group1,2,sd);
      sd=pooled_sd(sd0,sd1);
      d.sdz=d.raw/sd;
      sim=data.frame(n,d.pop=d,d.sdz,sd,d.raw,mean0,mean1,sd0,sd1,
                     row.names=NULL,stringsAsFactors=F);
      save_sim_fixd(sim,n,d);
    });
  return(T);
}
## het-d simulation
dosim_hetd=function(n,m,d,sd) {
  param(verbose);
  cases=expand.grid(n=n,d.het=d,sd.het=sd);
  if (nrow(cases)>0)
    apply(cases,1,function(row) {
      n=row[1]; d.het=row[2]; sd.het=row[3];
      if (verbose) print(paste(sep=' ','>>> dosim_hetd:',nvq(n,d.het,sd.het)));
      d.pop=rnorm(m,mean=d.het,sd=sd.het);
      group0=replicate(m,rnorm(n,mean=0));
      group1=vrnorm(n,mean=d.pop);
      mean0=colMeans(group0);
      mean1=colMeans(group1);
      d.raw=mean1-mean0;
      sd0=apply(group0,2,sd);
      sd1=apply(group1,2,sd);
      sd=pooled_sd(sd0,sd1);
      d.sdz=d.raw/sd;
      sim=data.frame(n,d.het,sd.het,d.pop,d.sdz,sd,d.raw,mean0,mean1,sd0,sd1,
                     row.names=NULL,stringsAsFactors=F);      
      save_sim_hetd(sim,n,d.het,sd.het);
    });
  return(T);
}

## --- Generate Analysis Data ---
#####
## mean significant effect size
## empirical from fixd simulation
domeand_fixd=function(n,d0,base='meand.fixd') {
  cases=expand.grid(n,d0);
  meand=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d0=case[2];
    if (verbose) print(paste(sep=' ','>>> domeand_fixd',nvq(n),nvq(d0)));
    sim=get_sim_fixd(n=n,d=d0);
    d.crit=d_crit(n=n);
    meand=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    data.frame(n,d0,meand,over=meand/d0,row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## empirical from rand simulation
domeand_rand=function(n,base='meand.rand') {
  meand=do.call(rbind,lapply(n,function(n) {
    if (verbose) print(paste(sep=' ','>>> domeand_rand',nvq(n)));
    sim=get_sim_rand(n=n);
    d.crit=d_crit(n=n);
    meand=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    data.frame(n,meand,over=meand/sim$d.pop,row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## empirical from hetd simulation
domeand_hetd=function(n,d.het,sd.het,base='meand.hetd') {
  cases=expand.grid(n,d.het,sd.het);
  meand=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d.het=case[2]; sd.het=case[3]
    if (verbose) print(paste(sep=' ','>>> domeand_hetd',nvq(n,d.het,sd.het)));
    sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
    d.crit=d_htcrit(n=n,sd.het=sd.het);
    meand=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    ## also compute with conventional pvalues
    d.crit=d_crit(n=n);
    meand.tval=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    data.frame(n,d.het,sd.het,meand,over=meand/d.het,meand.tval,over.tval=meand.tval/d.het,
               row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## theoretical from d2t sampling distribution
domeand_d2t=function(n,d0,base='meand.d2t') {
  param(verbose);
  cases=expand.grid(n,d0);
  meand=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d0=case[2];
    if (verbose) print(paste(sep=' ','>>> domeand_d2t',nvq(n),nvq(d0)));
    ## integral below equivalent to commmented out computation
    ## d=seq(d_crit(n),dmax,len=dlen);
    ## wt=d_d2t(n=n,d=d,d0=d0);
    ## meand=weighted.mean(d,wt);
    ##
    d.crit=d_crit(n);
    meand=integrate(function(d) d_d2t(n=n,d=d,d0=d0)*d,lower=d.crit,upper=Inf)$value / integrate(function(d) d_d2t(n=n,d=d,d0=d0),lower=d.crit,upper=Inf)$value;
    data.frame(n,d0,meand,over=meand/d0,row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## domeand_d2ht=function(n,d.het,sd.het,dmax=10,dlen=1e4) {
domeand_d2ht=function(n,d.het,sd.het,base='meand.d2ht') {
  param(verbose);
  if (verbose) print(paste(sep=' ','>>> domeand_d2ht'));
  ## NG 19-04-15: simply call meand_d2ht now that it exists
  cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
  meand=with(cases,meand_d2ht(n,d.het,sd.het));
  ## also compute with conventional pvalues
  meand.tval=with(cases,meand_d2ht(n,d.het,sd.het,d.crit=d_crit(n)));
  meand=data.frame(cases,meand,meand.tval);
  save_data(meand,base=base);
  invisible(meand);
}
## power
## TODO: decide which stubs make sense and implement them
## empirical from fixd simulation
dopower_fixd=function(n,d0,base='power.fixd') { }
## empirical from rand simulation
dopower_rand=function(n,base='power.rand') { }
## empirical from hetd simulation
dopower_hetd=function(n,d.het,sd.het,base='power.hetd') {
  cases=expand.grid(n,d.het,sd.het);
  power=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d.het=case[2]; sd.het=case[3]; 
    if (verbose) print(paste(sep=' ','>>> dopower_hetd',nvq(n,d.het,sd.het)));
    sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
    d.crit=d_htcrit(n=n,sd.het=sd.het);
    power=length(which(sim$d.sdz>=d.crit))/nrow(sim);
    ## also compute with conventional pvalues
    d.crit=d_crit(n=n);
    power.tval=length(which(sim$d.sdz>=d.crit))/nrow(sim);
    ## to compute over-estimate
    power.d2t=power_d2t(n,d.het);
    data.frame(n,d.het,sd.het,power,over=power/power.d2t,power.tval,over.tval=power.tval/power.d2t,
               row.names=NULL);
  }));
  save_data(power,base=base);
  invisible(power);
}
## theoretical from d2t sampling distribution
dopower_d2t=function(n,d0,base='power.d2t') { }
## theoretical from d2ht sampling distribution
dopower_d2ht=function(n,d.het,sd.het,base='power.d2ht') {
  if (param(verbose)) print('>>> dopower_d2ht');
  cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
  power=with(cases,power_d2ht(n,d.het,sd.het));
  ## also compute with conventional pvalues
  power.tval=with(cases,power_d2ht(n,d.het,sd.het,d.crit=d_crit(n)));
  ## to compute over-estimate
  power.d2t=with(cases,power_d2t(n,d.het));
  power=data.frame(cases,power,over=power/power.d2t,power.tval,over.tval=power.tval/power.d2t);
  save_data(power,base=base);
  invisible(power);
}

## pval
## TODO: decide which other functions make sense and implement them
## empirical from hetd simulation
dopval_hetd=function(n,sd.het,sig.level=param(sig.level),base='pval.hetd') {
  cases=expand.grid(n,sd.het);
  pval=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; sd.het=case[2]; 
    if (verbose) print(paste(sep=' ','>>> dopval_hetd',nvq(n,sd.het)));
    sim=get_sim_hetd(n=n,d=0,sd=sd.het);
    d.crit=d_crit(n=n);
    pval=length(which(abs(sim$d.sdz)>=d.crit))/nrow(sim);
    over=pval/sig.level;
    data.frame(n,sd.het,pval,over,row.names=NULL);
  }));
  save_data(pval,base=base);
  invisible(pval);
}
## theoretical from d2ht sampling distribution
dopval_d2ht=function(n,sd.het,sig.level=param(sig.level),base='pval.d2ht') {
  if (param(verbose)) print('>>> dopval_d2ht');
  cases=expand.grid(n=n,sd.het=sd.het);
  pval=with(cases,data.frame(n,sd.het,pval=d2htpval(n=n,sd.het=sd.het,d=d_crit(n))));
  pval$over=pval$pval/sig.level;
  save_data(pval,base=base);
  invisible(pval);
}
## confidence interval
## TODO: decide which other functions make sense and implement them
doci_d2ht=function(n,sd.het,d,conf.level=param(conf.level),base='ci.d2ht') {
  if (param(verbose)) print('>>> ci_d2ht');
  cases=expand.grid(n=n,sd.het=sd.het,d=d);
  ci=with(cases,t(ci_d2ht(n=n,sd.het=sd.het,d=d,conf.level=conf.level)));
  ci.d2t=with(cases,t(ci_d2t(n,d)));
  over=(ci[,2]-ci[,1])/(ci.d2t[,2]-ci.d2t[,1]);
  ## ci=data.frame(cases,ci.lo=ci[,1],ci.hi=ci[,2],ci.d2t.lo=ci.d2t[,1],ci.d2t.hi=ci.d2t[,2],over);
  ci=data.frame(cases,ci.lo=ci[,1],ci.hi=ci[,2],over);
  save_data(ci,base=base);
  invisible(ci);
}
