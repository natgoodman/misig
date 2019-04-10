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
    pval=d2pval(n,d.sdz);
    sim=data.frame(n,d.pop=d,d.sdz,sd,pval,d.raw,mean0,mean1,sd0,sd1,
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
      pval=d2pval(n,d.sdz);
      sim=data.frame(n,d.pop=d,d.sdz,sd,pval,d.raw,mean0,mean1,sd0,sd1,
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
      pval=d2pval(n,d.sdz);
      sim=data.frame(n,d.het,sd.het,d.pop,d.sdz,sd,pval,d.raw,mean0,mean1,sd0,sd1,
                     row.names=NULL,stringsAsFactors=F);
      
      save_sim_hetd(sim,n,d.het,sd.het);
    });
  return(T);
}

## --- Generate Analysis Data ---
## compute mean significant effect size
## empirical from fixd simulation
domeand_fixd=function(n,d0) {
  cases=expand.grid(n,d0);
  meand=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d0=case[2];
    if (verbose) print(paste(sep=' ','>>> domeand_fixd',nvq(n),nvq(d0)));
    sim=get_sim_fixd(n=n,d=d0);
    d.crit=d_crit(n=n);
    meand=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    data.frame(n,d0,meand,over=meand/d0,row.names=NULL);
  }));
  save_meand_empi(meand);
  invisible(meand);
}
## empirical from rand simulation
domeand_rand=function(n) {
  meand=do.call(rbind,lapply(n,function(n) {
    if (verbose) print(paste(sep=' ','>>> domeand_rand',nvq(n)));
    sim=get_sim_rand(n=n);
    d.crit=d_crit(n=n);
    meand=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    data.frame(n,meand,over=meand/d0,row.names=NULL);
  }));
  save_meand_empi(meand);
  invisible(meand);
}
## empirical from hetd simulation
domeand_hetd=function(n,d.het,sd.het) {
  cases=expand.grid(n,d.het,sd.het);
  meand=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d.het=case[2]; sd.het=case[3]
    if (verbose) print(paste(sep=' ','>>> domeand_hetd',nvq(n),nvq(d.het),nvq(sd.het)));
    sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
    d.crit=d_htcrit(n=n,sd.het=sd.het);
    meand=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    ## also compute with conventional pvalues
    d.crit=d_crit(n=n);
    meand.tval=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    data.frame(n,d.het,sd.het,meand,over=meand/d.het,meand.tval,over.tval=meand.tval/d.het,
               row.names=NULL);
  }));
  save_meand_empi(meand);
  invisible(meand);
}
## theoretical from d2t sampling distribution
## domeand_d2t=function(n,d0,dmax=10,dlen=1e4) {
domeand_d2t=function(n,d0) {
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
  save_meand_theo(meand);
  invisible(meand);
}
## theoretical from d2ht sampling distribution
## assumes interp already computed!
## domeand_d2ht=function(n,d.het,sd.het,dmax=10,dlen=1e4) {
domeand_d2ht=function(n,d.het,sd.het) {
  param(verbose);
  cases=expand.grid(n,d.het,sd.het);
  meand=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d.het=case[2]; sd.het=case[3]
    if (verbose) print(paste(sep=' ','>>> domeand_d2ht',nvq(n),nvq(d.het),nvq(sd.het)));
    ## integral below equivalent to commmented out computation
    ## d=seq(d_htcrit(n,d.het,sd.het),dmax,len=dlen);
    ## wt=d_d2ht(n=n,d=d,d.het=d.het,sd.het=sd.het)
    ## meand=weighted.mean(d,wt);
    d.crit=d_htcrit(n=n,sd.het=sd.het);
    meand=
      integrate(function(d)
        d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d)*d,lower=d.crit,upper=Inf)$value / integrate(function(d) d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d),lower=d.crit,upper=Inf)$value;
    ## also compute with conventional pvalues
    d.crit=d_crit(n);
    meand.tval=
      integrate(function(d)
        d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d)*d,lower=d.crit,upper=Inf)$value / integrate(function(d) d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d),lower=d.crit,upper=Inf)$value;
    data.frame(n,d.het,sd.het,meand,over=meand/d.het,meand.tval,over.tval=meand.tval/d.het,
               row.names=NULL);
  }));
  save_meand_theo(meand);
  invisible(meand);
}
