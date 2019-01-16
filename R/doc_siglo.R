#################################################################################
##
## Author:  Nat Goodman
## Created: 19-01-10
##          from repwr/R/doc_resig.R created 18-09-05
## Includes code from repwr/R/docfun_resig.R created 18-10-25
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate figures and tables for siglo blog post
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate Figures and Tables for siglo Blog Post ---
## no sections. only 2 figures
## n.fig is sample size for which figures plotted
doc_siglo=function(n.fig=20,sect=parent(sect,NULL)) {
  d.crit=d_crit(n.fig); d.pop=0.3;
  x='d.sdz'; y='d.pop';
  figblk_start();
  sim=get_sim(n=n.fig);
  title=title_siglo(n=n.fig,'P-values improve as observed effect size grows more extreme');
  dofig(plotdvsd,'big_picture',sim=sim,x=x,y=y,d.crit=d.crit,d.pop=d.pop,xlim=c(-2,2),
        title=title,cex.main=1);
  zoom_in=subset(sim,subset=near(d.sdz,d.crit,0.05));
  title=title_siglo(n=n.fig,'Sharp boundary between nonsignificant and significant p-values');
  dofig(plotdvsd,'zoom_in',sim=zoom_in,x=x,y=y,d.crit=d.crit,d.pop=d.pop,vlab=F,
        title=title,cex.main=1)
  figblk_end();
  ## support statements in text:
  support=data.frame(
    ##   d_crit(n=20)=0.64
    d.crit=d_crit(n=20),
    ##   power(n=20,d=0.3)=15%
    pwr.03.20=power.t.test(n=20,d=0.3)$power,
    ##   d2pval(n=20,d=0.3)=0.35
    pval.03.20=d2pval(n=20,d=0.3),
    ##   need n=87 for d=0.3 to be significant
    n.crit=uniroot(function(n) d2pval(n=n,d=0.3)-0.05,interval=c(20,100))$root,
    ##   power(d=0.3,power=.95)=290
    n.03.95=power.t.test(d=0.3,power=.95)$n,
    ##   power(d=0.5,power=.95)=105
    n.05.95=power.t.test(d=0.5,power=.95)$n);
  dotbl(support);
  invisible();
}
## generate title for doc_siglo
title_siglo=function(n,desc=NULL) {
  fig=paste(sep='','Figure ',figlabel());
  paste(collapse="\n",c(fig,desc,paste_nv(n)));
}
