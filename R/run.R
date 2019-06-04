#################################################################################
##
## Author:  Nat Goodman
## Created: 19-02-18
##          from ovrfx.R created 19-02-03 
##          from siglo.R created 19-01-01
##          from repwr/R/repwr.R created 17-10-05 
##           and repwr/R/sim.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Top level program for misig documents
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
source('R/source.R');                   # source the others

## ---- run ----
## run the program
## parameters defined in init
run=function(need.init=T,...) {
  if (need.init) wrap_fun(init,...);
  need.init=F;
  wrap_fun(dodat,...);             # generate data: eg, simulations, meand
  wrap_fun(dodoc,init_doc,...);     # generate figures, tables for doc
}
