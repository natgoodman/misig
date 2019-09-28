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
dosim_rand=function(n,m,d.gen,d.args) {
  param(m1,verbose);
  what='rand';
  sapply(n,function(n) {
    file=filename_sim(what,n);
    if (file.exists(file)&&(is.na(save)||!save)) {
      if (verbose) print(paste('>>> dosim_rand:',basename(file),'exists. skipping'));
      return();
    }
    more=m; i=1;
    while(more>0) {
      m1=min(m1,more);
      if (verbose) print(paste(sep=' ','>>> dosim_rand:',nvq(i,n)));
      group0=replicate(m1,rnorm(n,mean=0));
      d=do.call(d.gen,c(n=m1,d.args));
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
      save_sim_tmp(sim,i,n);
      more=more-m1; i=i+1;
    }
    ## consolidate subfiles into one
    sim=cat_sim(file=file);
  });
  return();
}
## fixed-d simulation
dosim_fixd=function(n,m,d) {
  param(m1,verbose);
  what='fixd';
  cases=expand.grid(n=n,d=d);
  if (nrow(cases)>0)
    withrows(cases,case,{
      file=filename_sim(what,n,d);
      if (file.exists(file)&&(is.na(save)||!save)) {
        if (verbose) print(paste('>>> dosim_fixd:',basename(file),'exists. skipping'));
        return();
      }
      more=m; i=1;
      while(more>0) {
        m1=min(m1,more);
        if (verbose) print(paste(sep=' ','>>> dosim_fixd:',nvq(i,n,d)));
        group0=replicate(m1,rnorm(n,mean=0));
        group1=replicate(m1,rnorm(n,mean=d));
        mean0=colMeans(group0);
        mean1=colMeans(group1);
        d.raw=mean1-mean0;
        sd0=apply(group0,2,sd);
        sd1=apply(group1,2,sd);
        sd=pooled_sd(sd0,sd1);
        d.sdz=d.raw/sd;
        sim=data.frame(n,d.pop=d,d.sdz,sd,d.raw,mean0,mean1,sd0,sd1,
                       row.names=NULL,stringsAsFactors=F);
        save_sim_tmp(sim,i,n,d);
        more=more-m1; i=i+1;
      }
      ## consolidate subfiles into one
      sim=cat_sim(file=file);
    });
  return();
}
## het-d simulation
dosim_hetd=function(n,m,d,sd) {
  param(m1,verbose);
  what='hetd';
  cases=expand.grid(n=n,d=d,sd=sd);
  if (nrow(cases)>0)
    withrows(cases,case,{
      file=filename_sim(what,n,d,sd);
      if (file.exists(file)&&(is.na(save)||!save)) {
        if (verbose) print(paste('>>> dosim_hetd:',basename(file),'exists. skipping'));
        return();
      }
      more=m; i=1;
      while(more>0) {
        m1=min(m1,more);
        if (verbose) print(paste(sep=' ','>>> dosim_hetd:',nvq(i,n,d,sd)));
        d.pop=rnorm(m1,mean=d,sd=sd);
        group0=replicate(m1,rnorm(n,mean=0));
        group1=vrnorm(n,mean=d.pop);
        mean0=colMeans(group0);
        mean1=colMeans(group1);
        d.raw=mean1-mean0;
        sd0=apply(group0,2,sd);
        sd1=apply(group1,2,sd);
        sd2=pooled_sd(sd0,sd1);
        d.sdz=d.raw/sd2;
        sim=data.frame(n,d.het=d,sd.het=sd,d.pop,d.sdz,sd=sd2,d.raw,mean0,mean1,sd0,sd1,
                       row.names=NULL,stringsAsFactors=F);      
        save_sim_tmp(sim,i,n,d,sd);
        more=more-m1; i=i+1;
      }
      ## consolidate subfiles into one
      sim=cat_sim(file=file);
    });
  return();
}

## --- Generate Analysis Data ---
########## meand
## mean significant effect size
#### fix model
## empirical from fixd simulation
domeand_fixd=function(n,d.pop,base='meand.fixd',verbose=param(verbose)) {
  cases=expand.grid(n=n,d.pop=d.pop);
  meand=do.call(rbind,withrows(cases,case,{
    if (verbose) print(paste(sep=' ','>>> domeand_fixd',nvq(n),nvq(d.pop)));
    sim=get_sim_fixd(n=n,d=d.pop);
    meand.all=mean(sim$d.sdz);
    d.crit=d_crit(n=n);
    meand=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    data.frame(n,d.pop,meand,over=meand/d.pop,meand.all,over.all=meand.all/d.pop,row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## empirical from rand simulation
domeand_rand=function(n,base='meand.rand',verbose=param(verbose)) {
  meand=do.call(rbind,lapply(n,function(n) {
    if (verbose) print(paste(sep=' ','>>> domeand_rand',nvq(n)));
    sim=get_sim_rand(n=n);
    meand.all=mean(sim$d.sdz);
    over.all=mean(sim$d.sdz/sim$d.pop);
    d.crit=d_crit(n=n);
    sim=subset(sim,subset=d.sdz>=d.crit);
    meand=mean(sim$d.sdz);
    over=mean(sim$d.sdz/sim$d.pop);
    data.frame(n,meand,over,meand.all,over.all,row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## theoretical from d2t sampling distribution
domeand_d2t=function(n,d.pop,base='meand.d2t',verbose=param(verbose)) {
  cases=expand.grid(n=n,d.pop=d.pop);
  meand=do.call(rbind,withrows(cases,case,{
    if (verbose) print(paste(sep=' ','>>> domeand_d2t',nvq(n),nvq(d.pop)));
    ## integral below equivalent to commmented out computation
    ## d=seq(d_crit(n),dmax,len=dlen);
    ## wt=d_d2t(n=n,d=d,d.pop=d.pop);
    ## meand=weighted.mean(d,wt);
    ##
    meand.all=integrate(function(d) d_d2t(n=n,d=d,d0=d.pop)*d,lower=-Inf,upper=Inf)$value;
    d.crit=d_crit(n);
    meand=integrate(function(d) d_d2t(n=n,d=d,d0=d.pop)*d,lower=d.crit,upper=Inf)$value / integrate(function(d) d_d2t(n=n,d=d,d0=d.pop),lower=d.crit,upper=Inf)$value;
    data.frame(n,d.pop,meand,over=meand/d.pop,meand.all,over.all=meand.all/d.pop,row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
#### het model
## empirical from hetd simulation
domeand_hetd=function(n,d.het,sd.het,base='meand.hetd',verbose=param(verbose)) {
  cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);  
  meand=do.call(rbind,withrows(cases,case,{
    if (verbose) print(paste(sep=' ','>>> domeand_hetd',nvq(n,d.het,sd.het)));
    sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
    meand.all=mean(sim$d.sdz);
    over.all=mean(sim$d.sdz/sim$d.pop);
    d.crit=d_htcrit(n=n,sd.het=sd.het);
    meand=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    ## also compute with conventional pvalues
    d.crit=d_crit(n=n);
    meand.tval=mean(sim$d.sdz[sim$d.sdz>=d.crit]);
    data.frame(n,d.het,sd.het,meand,meand.tval,meand.all,
               over=meand/d.het,over.tval=meand.tval/d.het,over.all=meand.all/d.het,
               row.names=NULL);
  }));
  save_data(meand,base=base);
  invisible(meand);
}
## domeand_d2ht=function(n,d.het,sd.het,dmax=10,dlen=1e4) {
domeand_d2ht=function(n,d.het,sd.het,base='meand.d2ht',verbose=param(verbose)) {
  if (verbose) print(paste(sep=' ','>>> domeand_d2ht'));
  ## NG 19-04-15: simply call meand_d2ht now that it exists
  cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
  meand.all=unlist(withrows(cases,case,{
    integrate(function(d) d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d)*d,lower=-Inf,upper=Inf)$value;
  }));
  meand=with(cases,meand_d2ht(n,d.het,sd.het));
  ## also compute with conventional pvalues
  meand.tval=with(cases,meand_d2ht(n,d.het,sd.het,d.crit=d_crit(n)));
  meand=data.frame(cases,meand,meand.tval,meand.all,
                   over=meand/d.het,over.tval=meand.tval/d.het,over.all=meand.all/d.het);
  save_data(meand,base=base);
  invisible(meand);
}
########## power
## power to detect d.pop.
#### fix model
## only useful to check concordance of simulation w/ theory
## empirical from fixd simulation
dopower_fixd=function(n,d.pop,base='power.fixd',verbose=param(verbose)) {
  cases=expand.grid(n=n,d.pop=d.pop);
  power=do.call(rbind,withrows(cases,case,{
    if (verbose) print(paste(sep=' ','>>> dopower_fixd',nvq(n,d.pop)));
    sim=get_sim_fixd(n=n,d=d.pop);
    d.crit=d_crit(n=n);
    power=length(which(sim$d.sdz>=d.crit))/nrow(sim);
    data.frame(case,power);
  }));
  save_data(power,base=base);
  invisible(power);
}
## theoretical from d2t sampling distribution
dopower_d2t=function(n,d.pop,base='power.d2t',verbose=param(verbose)) {
  if (verbose) print('>>> dopower_d2t');
  cases=expand.grid(n=n,d.pop=d.pop);
  power=with(cases,power_d2t(n,d.pop));
  power=data.frame(cases,power);
  save_data(power,base=base);
  invisible(power);
}
#### het model
## empirical from hetd simulation
dopower_hetd=function(n,d.het,sd.het,base='power.hetd',verbose=param(verbose)) {
  cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
  power=do.call(rbind,withrows(cases,case,{
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
#### fix model
## only useful to check concordance of simulation w/ theory
## empirical from fixd simulation
dopval_fixd=function(n,sig.level=param(sig.dat),base='pval.fixd',verbose=param(verbose)) {
  cases=expand.grid(n=n,sig.level=sig.level);
  pval=do.call(rbind,withrows(cases,case,{
    if (verbose) print(paste(sep=' ','>>> dopval_fixd',nvq(n,sig.level)));
    sim=get_sim_fixd(n=n,d=0);
    d.crit=d_crit(n=n,sig.level=sig.level);
    pval=length(which(abs(sim$d.sdz)>=d.crit))/nrow(sim);
    data.frame(case,pval,over=pval/sig.level);
  }));
  save_data(pval,base=base);
  invisible(pval);
}
## theoretical from d2t sampling distribution.
## not real useful... included for stylistic consistency
dopval_d2t=function(n,sig.level=param(sig.dat),base='pval.d2t',verbose=param(verbose)) {
  if (verbose) print('>>> dopval_d2t');
  cases=expand.grid(n=n,sig.level=sig.level);
  pval=with(cases,d2pval(n,d=d_crit(n,sig.level=sig.level)));
  pval=data.frame(cases,pval,over=pval/cases$sig.level);
  save_data(pval,base=base);
  invisible(pval);
}
#### het model
## empirical from hetd simulation
dopval_hetd=function(n,sd.het,sig.level=param(sig.dat),base='pval.hetd',verbose=param(verbose)) {
  cases=expand.grid(n=n,sig.level=sig.level,sd.het=sd.het);
  pval=do.call(rbind,withrows(cases,case,{
    if (verbose) print(paste(sep=' ','>>> dopval_hetd',nvq(n,sig.level,sd.het)));
    sim=get_sim_hetd(n=n,d=0,sd=sd.het);
    ## compute with ht pvals
    d.crit=d_htcrit(n=n,sd.het=sd.het,sig.level=sig.level);
    pval=length(which(abs(sim$d.sdz)>=d.crit))/nrow(sim);
    ## also compute with conventional pvalues
    d.crit=d_crit(n=n,sig.level=sig.level);
    pval.tval=length(which(abs(sim$d.sdz)>=d.crit))/nrow(sim);
    data.frame(case,pval,pval.tval,over=pval/sig.level,over.tval=pval.tval/sig.level);
  }));
  save_data(pval,base=base);
  invisible(pval);
}
## theoretical from d2ht sampling distribution
dopval_d2ht=function(n,sd.het,sig.level=param(sig.dat),base='pval.d2ht',verbose=param(verbose)) {
  if (verbose) print('>>> dopval_d2ht');
  cases=expand.grid(n=n,sig.level=sig.level,sd.het=sd.het);
  pval=with(cases,{
    pval=d2htpval(n=n,sd.het=sd.het,d=d_htcrit(n=n,sd.het=sd.het,sig.level=sig.level));
    pval.tval=d2htpval(n=n,sd.het=sd.het,d=d_crit(n,sig.level=sig.level));
    over=pval/sig.level;
    over.tval=pval.tval/sig.level;
    data.frame(pval,pval.tval,over,over.tval);
  });
  pval=data.frame(cases,pval);
  save_data(pval,base=base);
  invisible(pval);
}
########## confidence interval
## TODO 19-07-22: boy was I dumb.yes, it is possible to compute ci's from simulation
##   see new ci_sim function in stats.R
## theoretical tests ci width inflation when sig
## empiricals test coverage
#### fix model
## empirical from fixd simulation
doci_fixd=
  function(n,d.pop,conf.level=param(conf.level),m.ci=param(m.ci),base='ci.fixd',
           verbose=param(verbose)) {
    cases=expand.grid(n=n,d.pop=d.pop);
    ci=do.call(rbind,withrows(cases,case,{
      if (verbose) print(paste(sep=' ','>>> doci_fixd',nvq(n,d.pop)));
      sim=get_sim_fixd(n=n,d=d.pop);
      ## downsample if necessary
      if (!is.null(m.ci)&&m.ci<nrow(sim)) sim=sim[sample.int(nrow(sim),m.ci),];
      ci=with(sim,t(ci_d2t(n,d.sdz)));
      lo=ci[,1]; hi=ci[,2];
      sim$ci=hi-lo;
      meanci=mean(sim$ci);
      sim$cover=with(sim,between(d.pop,lo,hi));
      pval=with(sim,d2pval(n,d.sdz));
      cover=length(which(sim$cover))/nrow(sim);
      sim.sig=sim[pval<=0.05,];
      cover.sig=length(which(sim.sig$cover))/nrow(sim.sig);
      meanci.sig=mean(sim.sig$ci);
      data.frame(case,cover,cover.sig,
                 over=cover/conf.level,over.sig=cover.sig/conf.level,
                 meanci,meanci.sig);
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
      lo=ci[,1]; hi=ci[,2];
      sim$ci=hi-lo;
      meanci=mean(sim$ci);
      sim$cover=with(sim,between(d.pop,lo,hi));
      pval=with(sim,d2pval(n,d.sdz));
      cover=length(which(sim$cover))/nrow(sim);
      sim.sig=sim[pval<=0.05,];
      cover.sig=length(which(sim.sig$cover))/nrow(sim.sig);
      meanci.sig=mean(sim.sig$ci);
      data.frame(n,cover,cover.sig,
                 over=cover/conf.level,over.sig=cover.sig/conf.level,
                 meanci,meanci.sig);
    }));
    save_data(ci,base=base);
    invisible(ci);
  }
## theoretical from d2t sampling distribution. not real useful but why not??
doci_d2t=function(n,d.pop,conf.level=param(conf.level),base='ci.d2t',verbose=param(verbose)) {
  if (verbose) print('>>> doci_d2t');
  cases=expand.grid(n=n,d.pop=d.pop);
  ci=with(cases,t(ci_d2t(n=n,d=d.pop,conf.level=conf.level)));
  lo=ci[,1]; hi=ci[,2]; ci=hi-lo;
  ci=data.frame(cases,lo,hi,ci);
  save_data(ci,base=base);
  invisible(ci);
}
#### het model
## empirical from hetd simulation
doci_hetd=
  function(n,d.het,sd.het,conf.level=param(conf.level),m.ci=param(m.ci),base='ci.hetd',
           verbose=param(verbose)) {
    cases=expand.grid(n=n,d.het=d.het,sd.het=sd.het);
    ci=do.call(rbind,withrows(cases,case,{
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
## theoretical from d2ht sampling distribution
doci_d2ht=
  function(n,d.het,sd.het,conf.level=param(conf.level),base='ci.d2ht',verbose=param(verbose)) {
    if (verbose) print('>>> ci_d2ht');
    cases=expand.grid(n=n,sd.het=sd.het,d.het=d.het);
    ci=with(cases,t(ci_d2ht(n=n,sd.het=sd.het,d=d.het,conf.level=conf.level)));
    ci.d2t=with(cases,t(ci_d2t(n,d=d.het)));
    lo=ci[,1]; hi=ci[,2]; ci=hi-lo;
    lo.tval=ci.d2t[,1]; hi.tval=ci.d2t[,2]; ci.tval=hi.tval-lo.tval;
    over=ci/ci.tval;
    ci=data.frame(cases,lo,hi,ci,lo.tval,hi.tval,ci.tval,over);
    save_data(ci,base=base);
    invisible(ci);
  }
