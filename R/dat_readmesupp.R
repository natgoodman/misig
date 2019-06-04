#################################################################################
##
## Author:  Nat Goodman
## Created: 19-05-09 
##          from dat_readme.R created 19-02-19
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate data for readme supplement document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Data Generation for README supplement ----
dat_readmesupp=function(...) {
  ## do all types of simulations
  param(n.fixd,m.fixd,d.fixd);
  dosim_fixd(n.fixd,m.fixd,d.fixd);
  param(n.rand,m.rand,d.rand);
  dosim_rand(n.rand,m.rand,d.rand);
  param(n.hetd,m.hetd,d.hetd,sd.hetd);
  dosim_hetd(n.hetd,m.hetd,d.hetd,sd.hetd);

  ## run all dat.R data generation functions
  ## compute four statistics - meand, power, pval, ci
  ##   on five data sources - fixd, rand, hetd, d2t, d2ht
  #### meand
  domeand_fixd(n.fixd,d.fixd)
  domeand_rand(n.rand);
  domeand_hetd(n.hetd,d.hetd,sd.hetd);
  domeand_d2t(n.fixd,d.fixd);
  domeand_d2ht(n.hetd,d.hetd,sd.hetd);
  #### power
  dopower_fixd(n.fixd,d.fixd);
  ## dopower_rand(n); # makes no sense - each instance has different d.pop
  dopower_hetd(n.hetd,d.hetd,sd.hetd);
  dopower_d2t(n.fixd,d.fixd);
  dopower_d2ht(n.hetd,d.hetd,sd.hetd);
  #### pval
  dopval_fixd(n.fixd);
  ## dopval_rand(n); # makes no sense - need many instances with d.pop=0
  dopval_hetd(n.hetd,sd.hetd);
  dopval_d2t(n.fixd);
  dopval_d2ht(n.hetd,sd.hetd);
  #### ci
  doci_fixd(n.fixd,d.fixd);
  doci_rand(n.rand);
  doci_hetd(n.hetd,d.hetd,sd.hetd);
  doci_d2t(n.fixd,d.fixd);
  doci_d2ht(n.hetd,d.hetd,sd.hetd);
  
  invisible();
}

