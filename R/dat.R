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
dosim_rand=function(n,m,d,verbose=param(verbose)) {
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
dosim_fixd=function(n,m,d,verbose=param(verbose)) {
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
dosim_hetd=function(n,m,d,sd,verbose=param(verbose)) {
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
########## meand
## mean significant effect size
## empirical from fixd simulation
domeand_fixd=function(n,d0,base='meand.fixd',verbose=param(verbose)) {
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
domeand_rand=function(n,base='meand.rand',verbose=param(verbose)) {
  meand=do.call(rbind,lapply(n,function(n) {
    if (verbose) print(paste(sep=' ','>>> domeand_rand',nvq(n)));
    sim=get_sim_rand(n=n);
    d.crit=d_crit(n=n);
    sim=subset(sim,subset=d.sdz>=d.crit);
    meand=mean(sim$d.sdz);
    over=mean(sim$d.sdz/sim$d.pop);
    data.frame(n,meand,over,row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## empirical from hetd simulation
domeand_hetd=function(n,d.het,sd.het,base='meand.hetd',verbose=param(verbose)) {
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
    data.frame(n,d.het,sd.het,meand,meand.tval,over=meand/d.het,over.tval=meand.tval/d.het,
               row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## theoretical from d2t sampling distribution
domeand_d2t=function(n,d0,base='meand.d2t',verbose=param(verbose)) {
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
domeand_d2ht=function(n,d.het,sd.het,base='meand.d2ht',verbose=param(verbose)) {
  if (verbose) print(paste(sep=' ','>>> domeand_d2ht'));
  ## NG 19-04-15: simply call meand_d2ht now that it exists
  cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
  meand=with(cases,meand_d2ht(n,d.het,sd.het));
  ## also compute with conventional pvalues
  meand.tval=with(cases,meand_d2ht(n,d.het,sd.het,d.crit=d_crit(n)));
  meand=data.frame(cases,meand,meand.tval);
  meand=cbind(meand,with(meand,data.frame(over=meand/d.het,over.tval=meand.tval/d.het)));
  save_data(meand,base=base);
  invisible(meand);
}
########## power
## empirical from fixd simulation
dopower_fixd=function(n,d0,base='power.fixd',verbose=param(verbose)) {
  cases=expand.grid(n=n,d0=d0);
  power=do.call(rbind,apply(cases,1,function(case) {
    n=case['n']; d0=case['d0'];
    if (verbose) print(paste(sep=' ','>>> dopower_fixd',nvq(n,d0)));
    sim=get_sim_fixd(n=n,d=d0);
    d.crit=d_crit(n=n);
    power=length(which(sim$d.sdz>=d.crit))/nrow(sim);
    power.d2t=power_d2t(n,d0);
    data.frame(case,power,power.d2t);
  }));
  save_data(power,base=base);
  invisible(power);
}
## empirical from rand simulation - makes no sense - each instance has different d.pop
dopower_rand=function(n,base='power.rand',verbose=param(verbose)) { }
## empirical from hetd simulation
dopower_hetd=function(n,d.het,sd.het,base='power.hetd',verbose=param(verbose)) {
  cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
  power=do.call(rbind,apply(cases,1,function(case) {
    n=case['n']; d.het=case['d.het']; sd.het=case['sd.het'];
    if (verbose) print(paste(sep=' ','>>> dopower_hetd',nvq(n,d.het,sd.het)));
    sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
    ## compute with ht pvals
    d.crit=d_htcrit(n=n,sd.het=sd.het);
    power=length(which(sim$d.sdz>=d.crit))/nrow(sim);
    ## also compute with conventional pvalues
    d.crit=d_crit(n=n);
    power.tval=length(which(sim$d.sdz>=d.crit))/nrow(sim);
    ## to compute over-estimate
    power.d2t=power_d2t(n,d.het);
    data.frame(case,power,power.tval,power.d2t,
               over=power/power.d2t,over.tval=power.tval/power.d2t);
  }));
  save_data(power,base=base);
  invisible(power);
}
## theoretical from d2t sampling distribution. not real useful but why not??
dopower_d2t=function(n,d0,base='power.d2t',verbose=param(verbose)) {
  if (verbose) print('>>> dopower_d2t');
  cases=expand.grid(n=n,d0=d0);
  power=with(cases,power_d2t(n,d0));
  power=data.frame(cases,power);
  save_data(power,base=base);
  invisible(power);
}
## theoretical from d2ht sampling distribution
dopower_d2ht=function(n,d.het,sd.het,base='power.d2ht',verbose=param(verbose)) {
  if (verbose) print('>>> dopower_d2ht');
  cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
  power=with(cases,power_d2ht(n,d.het,sd.het));
  ## also compute with conventional pvalues
  power.tval=with(cases,power_d2ht(n,d.het,sd.het,d.crit=d_crit(n)));
  ## to compute over-estimate
  power.d2t=with(cases,power_d2t(n,d.het));
  power=data.frame(cases,power,power.tval,power.d2t,
                   over=power/power.d2t,over.tval=power.tval/power.d2t);
  save_data(power,base=base);
  invisible(power);
}
########## pval
## empirical from fixd simulation
dopval_fixd=function(n,sig.level=param(sig.level),base='pval.fixd',verbose=param(verbose)) {
  pval=sapply(n,function(n) {
    if (verbose) print(paste(sep=' ','>>> dopval_fixd',nvq(n)));
    sim=get_sim_fixd(n=n,d=0);
    d.crit=d_crit(n=n);
    length(which(abs(sim$d.sdz)>=d.crit))/nrow(sim);
  });
  pval=data.frame(n=n,pval,over=pval/sig.level);
  save_data(pval,base=base);
  invisible(pval);
}
## empirical from rand simulation - makes no sense - need many instances with d.pop=0
dopval_rand=function(n,base='pval.rand',verbose=param(verbose)) { }
## empirical from hetd simulation
dopval_hetd=function(n,sd.het,sig.level=param(sig.level),base='pval.hetd',verbose=param(verbose)) {
  cases=expand.grid(n=n,sd.het=sd.het);
  pval=do.call(rbind,apply(cases,1,function(case) {
    n=case['n']; sd.het=case['sd.het'];
    if (verbose) print(paste(sep=' ','>>> dopval_hetd',nvq(n,sd.het)));
    sim=get_sim_hetd(n=n,d=0,sd=sd.het);
    ## compute with ht pvals
    d.crit=d_htcrit(n=n,sd.het=sd.het);
    pval=length(which(abs(sim$d.sdz)>=d.crit))/nrow(sim);
    ## also compute with conventional pvalues
    d.crit=d_crit(n=n);
    pval.tval=length(which(abs(sim$d.sdz)>=d.crit))/nrow(sim);
    data.frame(case,pval,pval.tval,over=pval/sig.level,over.tval=pval.tval/sig.level);
  }));
  save_data(pval,base=base);
  invisible(pval);
}
## theoretical from d2t sampling distribution. not real useful but why not??
dopval_d2t=function(n,sig.level=param(sig.level),base='pval.d2t',verbose=param(verbose)) {
  if (verbose) print('>>> dopval_d2t');
  pval=d2pval(n,d=d_crit(n));
  pval=data.frame(n=n,pval,over=pval/sig.level);
  save_data(pval,base=base);
  invisible(pval);
}
## theoretical from d2ht sampling distribution
dopval_d2ht=function(n,sd.het,sig.level=param(sig.level),base='pval.d2ht',verbose=param(verbose)) {
  if (verbose) print('>>> dopval_d2ht');
  cases=expand.grid(n=n,sd.het=sd.het);
  pval=with(cases,d2htpval(n=n,sd.het=sd.het,d=d_htcrit(n=n,sd.het=sd.het)));
  pval.tval=with(cases,d2htpval(n=n,sd.het=sd.het,d=d_crit(n)));
  pval=data.frame(cases,pval,pval.tval,over=pval/sig.level,over.tval=pval.tval/sig.level);
  save_data(pval,base=base);
  invisible(pval);
}

########## confidence interval
## empiricals test coverage
## empirical from fixd simulation
doci_fixd=
  function(n,d0,conf.level=param(conf.level),m.ci=param(m.ci),base='ci.fixd',
           verbose=param(verbose)) {
    cases=expand.grid(n=n,d.pop=d0);
    ci=do.call(rbind,apply(cases,1,function(case) {
      n=case['n']; d.pop=case['d.pop'];
      if (verbose) print(paste(sep=' ','>>> doci_fixd',nvq(n,d.pop)));
      sim=get_sim_fixd(n=n,d=d.pop);
      ## downsample if necessary
      if (!is.null(m.ci)&&m.ci<nrow(sim)) sim=sim[sample.int(nrow(sim),m.ci),];
      ci=with(sim,t(ci_d2t(n,d.sdz)));
      sim$lo=ci[,1]; sim$hi=ci[,2];
      sim$cover=with(sim,between(d.pop,lo,hi));
      sim$pval=with(sim,d2pval(n,d.sdz));
      cover=length(which(sim$cover))/nrow(sim);
      sim.sig=subset(sim,subset=pval<=0.05);
      cover.sig=length(which(sim.sig$cover))/nrow(sim.sig);
      data.frame(case,cover,cover.sig,over=cover/conf.level,over.sig=cover.sig/conf.level);
    }));
    save_data(ci,base=base);
    invisible(ci);
  }
## empirical from rand simulation
doci_rand=function(n,m.ci=param(m.ci),base='ci.rand',verbose=param(verbose)) {
  ci=do.call(rbind,lapply(n,function(n) {
    if (verbose) print(paste(sep=' ','>>> doci_rand',nvq(n)));
    sim=get_sim_rand(n=n);
    ## downsample if necessary
    if (!is.null(m.ci)&&m.ci<nrow(sim)) sim=sim[sample.int(nrow(sim),m.ci),];
    ci=with(sim,t(ci_d2t(n,d.sdz)));
    sim$lo=ci[,1]; sim$hi=ci[,2];
    sim$cover=with(sim,between(d.pop,lo,hi));
    sim$pval=with(sim,d2pval(n,d.sdz));
    cover=length(which(sim$cover))/nrow(sim);
    sim.sig=subset(sim,subset=pval<=0.05);
    cover.sig=length(which(sim.sig$cover))/nrow(sim.sig);
    data.frame(n,cover,cover.sig,over=cover/conf.level,over.sig=cover.sig/conf.level);
  }));
  save_data(ci,base=base);
  invisible(ci);
}
## empirical from hetd simulation
doci_hetd=
  function(n,d.het,sd.het,conf.level=param(conf.level),m.ci=param(m.ci),base='ci.hetd',
           verbose=param(verbose)) {
    cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
    ci=do.call(rbind,apply(cases,1,function(case) {
      n=case['n']; d.het=case['d.het']; sd.het=case['sd.het'];
      if (verbose) print(paste(sep=' ','>>> doci_hetd',nvq(n,d.het,sd.het)));
      sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
      ## downsample if necessary
      if (!is.null(m.ci)&&m.ci<nrow(sim)) sim=sim[sample.int(nrow(sim),m.ci),];
      ci=with(sim,t(ci_d2ht(n,sd.het,d=d.sdz)));
      sim$lo=ci[,1]; sim$hi=ci[,2];
      sim$cover=with(sim,between(d.het,lo,hi));
      cover=length(which(sim$cover))/nrow(sim);
      tval=with(sim,t(ci_d2t(n,d=d.sdz)));
      sim$tval.lo=tval[,1]; sim$tval.hi=tval[,2];
      sim$cover.tval=with(sim,between(d.het,tval.lo,tval.hi));
      cover.tval=length(which(sim$cover.tval))/nrow(sim);

      sim$pval=with(sim,d2htpval(n,sd.het,d.sdz));
      sim.sig=subset(sim,subset=pval<=0.05);
      cover.sig=length(which(sim.sig$cover))/nrow(sim.sig);

      sim$tval=with(sim,d2pval(n,d.sdz));
      sim.tsig=subset(sim,subset=tval<=0.05);
      cover.tsig=length(which(sim.tsig$cover.tval))/nrow(sim.tsig);
      data.frame(case,cover,cover.tval,cover.sig,cover.tsig,
                 over=cover/conf.level,over.tval=cover.tval/conf.level,
                 over.sig=cover.sig/conf.level,over.tsig=cover.tsig/conf.level);
    }));
    save_data(ci,base=base);
    invisible(ci);
  }
## theoretical from d2t sampling distribution. not real useful but why not??
doci_d2t=function(n,d0,conf.level=param(conf.level),base='ci.d2t',verbose=param(verbose)) {
  if (verbose) print('>>> doci_d2t');
  cases=expand.grid(n=n,d0=d0);
  ci=with(cases,t(ci_d2t(n=n,d=d0,conf.level=conf.level)));
  lo=ci[,1]; hi=ci[,2]; ci=hi-lo;
  ci=data.frame(cases,lo,hi,ci);
  save_data(ci,base=base);
  invisible(ci);
}
## theoretical from d2ht sampling distribution
doci_d2ht=function(n,d,sd.het,conf.level=param(conf.level),base='ci.d2ht',verbose=param(verbose)) {
  if (verbose) print('>>> ci_d2ht');
  cases=expand.grid(n=n,sd.het=sd.het,d=d);
  ci=with(cases,t(ci_d2ht(n=n,sd.het=sd.het,d=d,conf.level=conf.level)));
  ci.d2t=with(cases,t(ci_d2t(n,d)));
  lo=ci[,1]; hi=ci[,2]; ci=hi-lo;
  tval.lo=ci.d2t[,1]; tval.hi=ci.d2t[,2]; tval.ci=tval.hi-tval.lo;
  over=ci/tval.ci;
  ci=data.frame(cases,lo,hi,ci,tval.lo,tval.hi,tval.ci,over);
  save_data(ci,base=base);
  invisible(ci);
}
