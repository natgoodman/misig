#################################################################################
##
## Author:  Nat Goodman
## Created: 19-02-19
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate data for readme document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Data Generation for README ----
dat_readme=function(...) {
  ## do all types of simulations
  param(n.rand,m.rand,d.gen,d.args);
  dosim_rand(n.rand,m.rand,d.gen,d.args);
  param(n.fixd,m.fixd,d.fixd);
  dosim_fixd(n.fixd,m.fixd,d.fixd);
  param(n.hetd,m.hetd,d.hetd,sd.hetd);
  dosim_hetd(n.hetd,m.hetd,d.hetd,sd.hetd);

  ## fixed effect scenario - mean significant effect size, aka meand
  ## simulation (fixd, rand)
  domeand_fixd(n.fixd,d.fixd);
  domeand_rand(n.rand);
  ## theoretical 
  domeand_d2t(n.fixd,d.fixd);

  ## het effect scenario - pval, meand
  ## simulation (hetd)
  param(sig.dat);
  dopval_hetd(n.hetd,sd.hetd,sig.dat);
  ## doci_hetd(n.hetd,d.hetd,sd.hetd);
  domeand_hetd(n.hetd,d.hetd,sd.hetd);
  ## theoretical
  dopval_d2ht(n.hetd,sd.hetd,sig.dat);
  domeand_d2ht(n.hetd,d.hetd,sd.hetd);

  invisible();
}

