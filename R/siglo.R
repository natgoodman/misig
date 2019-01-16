#################################################################################
##
## Author:  Nat Goodman
## Created: 19-01-01
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
source('R/datman.R');
 source('R/doc.R');
## source('R/doc_readme.R');
source('R/doc_siglo.R');
source('R/init.R');
source('R/plot.R');
source('R/sim.R');
source('R/stats.R');
source('R/util.R');

## ---- run ----
## run the program
## parameters defined in init
run=function(need.init=T,...) {
  if (need.init) wrap_fun(init,...);
  need.init=F;
  dosim();                        # do the simulation
  wrap_fun(dodoc,init_doc,...);   # generate figures, tables for doc
}
