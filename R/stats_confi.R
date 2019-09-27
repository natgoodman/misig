#################################################################################
##
## Author:  Nat Goodman
## Created: 19-08-04
##          from doc_confi.R created 19-07-16 & bayes_confi.R created 19-07-21
##            with additional content from stats.R
##          doc_confi.R created from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Specialized stats functions for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
########################################
## some of these functions may get moved to stats.R when suitably generalized
## empirical distribution from sim
## NG 19-07-22: VERY INCOMPLETE! just realized I needed this. sigh
##   quick web search suggests this is more complicated than it looks...
## data is vector 
d_sim=function(data) stop('d_sim not yet implemented')
p_sim=function(data,d,...) ecdf(data,...)(d);
q_sim=function(data,p,...) quantile(data,probs=p,...)
r_sim=function(m,data) as.vector(q_sim(data,runif(m,0,1)))

## this one is equivalent to the one below
## ci_sim_=function(data,d,conf.level) {
##   p0=(1-conf.level)/2; p1=1-p0;
##   p_sim=ecdf(data);
##   lo=suppressWarnings(uniroot(function(d0) p_sim(d0)-p0,interval=c(-10,10))$root);
##   hi=suppressWarnings(uniroot(function(d0) p_sim(d0)-p1,interval=c(-10,10))$root);
##   c(lo,hi);
## }
## ci_sim=Vectorize(ci_sim_,vectorize.args=c('d','conf.level'));

## this one adapted from ci_norm. suffix 'q' because uses q_sim
ci_simq=function(data,simplify=T,conf.level=param(conf.level)) {
  p0=(1-conf.level)/2; p1=1-p0;
  lo=quantile(data,probs=p0,type=1);
  hi=quantile(data,probs=p1,type=1);
  ci=rbind(lo,hi);
  colnames(ci)=conf.level;
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
## this one is from first principles
ci_sim=function(data,d0,simplify=T,conf.level=param(conf.level)) {
  p0=1-conf.level; p1=1-p0;
  data.lt=data[data<=d0];
  data.gt=data[data>=d0];
  lo=quantile(data.lt,probs=p0,type=1);
  hi=quantile(data.gt,probs=p1,type=1);
  ci=rbind(lo,hi);
  colnames(ci)=conf.level;
  if (simplify&ncol(ci)==1) ci=as.vector(ci);
  ci;
}
