#################################################################################
##
## Author:  Nat Goodman
## Created: 19-02-03 
##          from siglo.R 19-01-01
##          from repwr/R/repwr.R created 17-10-05 
##           and repwr/R/sim.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Explore significance bias when power is low
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
source('R/util.R');
source('R/datman.R');
source('R/doc.R');
## source('R/doc_readme.R');
source('R/doc_ovrfx.R');
source('R/init.R');
source('R/plot.R');
## source('R/sim.R');
source('R/stats.R');

## ---- run ----
## run the program
## parameters defined in init
run=function(need.init=T,doc='ovrfx',...) {
  if (need.init) wrap_fun(init,...);
  need.init=F;
  dosim();                        # do the simulation
  domeand();
  wrap_fun(dodoc,init_doc,...);   # generate figures, tables for doc
}
## ---- Simulation Functions ----
## do the simulation
dosim=function() {
  ## do simulations! each value of n x d
  param(n.rand,m.rand,d.rand);
  dosim_rand(n.rand,m.rand,d.rand);
  param(n.fixd,m.fixd,d.fixd);
  dosim_fixd(n.fixd,m.fixd,d.fixd);
  invisible();
}
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
  cases=expand.grid(n,d);
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
## --- Generate Analysis Data ---
## compute mean significant effect size
## theoretical calculation NOT USED in ovrfx

domeand=function() {
  param(n.fixd,d.fixd);
  domeand_empi(n.fixd,d.fixd);
  ## param(n.meand,d0.meand,verbose);
  ## domeand_theo(n.meand,d0.meand);
  }
domeand_empi=function(n,d0,breaks=100) {
  cases=expand.grid(n,d0);
  meand=do.call(rbind,apply(cases,1,function(row) {
    n=row[1]; d0=row[2];
    sim=get_sim_fixd(n=n,d=d0);
    hist=hist(sim$d.sdz,breaks=breaks,plot=F);
    sig=which(hist$mids>=d_crit(n));
    meand=with(hist,weighted.mean(mids[sig],counts[sig]));
    data.frame(n,d0,meand,over=meand/d0,row.names=NULL);
  }))
  save_meand_empi(meand);
  invisible(meand);
}
domeand_theo=function(n.d0,dmax=10,dlen=1e4) {
  param(verbose);
  d0=round(d0.meand,digits=2);     # to avoid problems with imprecise decimals
  cases=expand.grid(n.meand,d0);
  meand=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d0=case[2]
    if (verbose) print(paste(sep=' ','>>> domeand',nvq(n),nvq(d0)));
    d=seq(d_crit(n),dmax,len=dlen);
    ## wt=d_d2t(n=n,d=d,d0=d0);
    ## NG 19-02-10: it seems I messed up d_d2t.
    ##   needs to be scaled by d2t for reason I don't yet understand...
    ## temporary hack until I fix stats:d_d2t
    ## BUT doesn't affect answer since it's just linear scaling of weights...
    wt=d_d2t(n=n,d=d,d0=d0)*sqrt(n/2);
    meand=weighted.mean(d,wt);
    data.frame(n,d0,meand,over=meand/d0,row.names=NULL);
  }));
  save_meand_theo(meand);
  invisible(meand);
}
