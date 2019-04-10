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
## Generate data for ovrht document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Data Generation for ovrht ----
dat_ovrht=function() {
  param(n.hetd,m.hetd,d.hetd,sd.hetd);
  dosim_hetd(n.hetd,m.hetd,d.hetd,sd.hetd);
  ## generate pval and ci tables
  n=seq(20,400,by=10);
  sd.het=c(0,0.05,0.1,0.2,0.4);
  cases=expand.grid(n=n,sd.het=sd.het);
  pval=with(cases,data.frame(n,sd.het,pval=d2htpval(n=n,sd.het=sd.het,d=d_crit(n))));
  pval$infl=round(digits=2,pval$pval/0.05)
  pval$pval=round(digits=2,pval$pval)
  save_data(pval);
  if (verbose) print('>>> ci');
  ci=with(cases,data.frame(n,sd.het,ci=ci_d2ht(n=n,sd.het=sd.het,d=0)[2,]));
  ci$infl=round(digits=2,ci$ci/ci_d2t(n=cases$n,d=0)[2,]);
  ci$ci=round(digits=2,ci$ci);
  save_data(ci);
  invisible();
}
