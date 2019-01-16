#################################################################################
##
## Author:  Nat Goodman
## Created: 19-01-01
##          from repwr/R/sim.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate simulated data for effit
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
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
