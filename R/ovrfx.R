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
  param(n);
  sapply(n,dosim1);
  invisible();
}
vrnorm=Vectorize(rnorm,"mean");
dosim1=function(n) {
  param(m,d,verbose);
   ## use saved sim if exists and args permit
  sim=get_sim(n,must.exist=F);
  if (!is.null(sim)) return(invisible(sim));
  ## no saved simulation or args say not to use it. run simulation
  if (verbose) print(paste(sep=' ','>>> dosim1:',nvq(n)));
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
  ## optionally save sim and add to in-memory list. usually don't keep - no point
  save_sim(sim,n);
  invisible(sim);
}

## --- Generate Analysis Data ---
## compute mean significant effect size
domeand=function(dmax=10,dlen=1e4) {
  ## use saved meand if exists and args permit
  meand=get_data(meand,must.exist=F);
  if (!is.null(meand)) return(invisible(meand));
  ## no saved meand or args say not to use it. calculate it
  param(n.meand,d0.meand,verbose);
  d0=round(d0.meand,digits=2);     # to avoid problems with imprecise decimals
  cases=expand.grid(n.meand,d0);
  meand=do.call(rbind,apply(cases,1,function(case) {
    n=case[1]; d0=case[2]
    if (verbose) print(paste(sep=' ','>>> domeand',nvq(n),nvq(d0)));
    d=seq(d_crit(n),dmax,len=dlen);
    wt=d_d2t(n=n,d=d,d0=d0);
    meand=weighted.mean(d,wt);
    data.frame(n,d0,meand,over=meand/d0,row.names=NULL);
  }));
  ## optionally save meand and add to in-memory list. usually do - needed for dodetl
  save_data(meand);
  invisible(meand);
}
