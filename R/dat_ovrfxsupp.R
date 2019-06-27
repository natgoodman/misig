#################################################################################
##
## Author:  Nat Goodman
## Created: 19-06-27
##          from dat_ovrhtsupp.R 19-04-09
##          from dat_ovrht.R created 19-02-19
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
## ---- Data Generation for ovrhtsupp ----
dat_ovrhtsupp=function() {
  param(n.hetd,m.hetd,d.hetd,sd.hetd);
  dosim_hetd(n.hetd,m.hetd,d.hetd,sd.hetd);
  ## mean significant effect size
  domeand_hetd(n.hetd,d.hetd,sd.hetd);
  domeand_d2ht(n.hetd,d.hetd,sd.hetd);
  ## power
  dopower_hetd(n.hetd,d.hetd,sd.hetd);
  dopower_d2ht(n.hetd,d.hetd,sd.hetd);
  ## pval
  dopval_hetd(n.hetd,sd.hetd);
  dopval_d2ht(n.hetd,sd.hetd);
  ## ci
  doci_hetd(n.hetd,d.hetd,sd.hetd);
  doci_d2ht(n.hetd,d.hetd,sd.hetd);

  invisible();
}
