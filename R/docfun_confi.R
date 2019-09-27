#################################################################################
##
## Author:  Nat Goodman
## Created: 19-08-04
##          from doc_confi.R created 19-07-16
##          from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Specialized data generation functions for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
distr.all=cq(std,sim,simq,d2tpost,normpost,meta,bayes,d2t,norm,prior);
## compute ci stats (lo,hi,coverage) from sim
doci_confi=
  function(sim,n,d0,distr=cq(std,sim,simq,d2tpost,normpost,meta,bayes),mean.prior,sd.prior,prior,
           conf.level=param(conf.level)) {
    if (missing(distr)) distr='std'
    else distr=match.arg(distr,several.ok=TRUE);
    if (('meta' %in% distr)&&(missing(mean.prior)||missing(sd.prior)))
      stop("distr contains 'meta' but 'mean.prior' or 'sd.prior' arguments missing");
    if (('bayes' %in% distr)&&missing(prior))
      stop("distr contains 'bayes' but prior argument missing");
    sim.byn=split(sim,sim$n);
    ci=do.call(rbind,lapply(n,function(n) {
      sim=sim.byn[[as.character(n)]];
      d.pop=sim$d.pop;
      m=length(d.pop);
      if ('meta' %in% distr) meta=meta_normd2t(mean.prior,sd.prior,d0,n);
      if ('bayes' %in% distr) init_bayes(n=n,d0=d0,prior=prior);
      ci=do.call(rbind,lapply(distr,function(distr) {
        ci=switch(
          distr,
          std=doci1(ci_d2t,list(n=n,d=d0),d.pop,conf.level),
          sim=doci1(ci_sim,list(data=d.pop,d0=d0),d.pop,conf.level),
          simq=doci1(ci_simq,list(data=d.pop),d.pop,conf.level),
          d2tpost=doci1(ci_d2tpost,list(n=n,d0=d0),d.pop,conf.level),
          normpost=doci1(ci_norm,list(mean=d0,sd=sqrt(2/n)),d.pop,conf.level),
          meta=doci1(ci_norm,list(mean=meta$mean,sd=meta$sd),d.pop,conf.level),
          bayes=doci1(ci_bayes,list(),d.pop,conf.level),
          stop(paste0('Unknown distribution: ',distr,'. should have been caught earlier!')));
        colnames(ci)=cq(lo,hi,cvr);
        data.frame(distr,ci,stringsAsFactors=F);
      }));
      data.frame(n=n,conf.level,ci);
    }));
    ci;
  }
doci1=function(ci_fun,ci.args,d.pop,conf.level) {
  ci.args=c(ci.args,list(simplify=F,conf.level=conf.level));
  ci=do.call(ci_fun,ci.args);
  lo=ci[1,]; hi=ci[2,];
  cvr=do.call(rbind,lapply(d.pop,function(d.pop) between(d.pop,lo,hi)));
  cvr=colSums(cvr)/length(d.pop);
  cbind(lo,hi,cvr);
}
## compute qq stats (quantile) from sim and distributions
## prior is density, as usual. qprior is quantile
## TODO: if prior useful, generate quantile function from density
doqq_confi=
  function(sim,n,d0,distr=cq(sim,d2tpost,normpost,meta,bayes,prior),mean.prior,sd.prior,
           prior,qprior,rprior,p=seq(0.01,0.99,by=0.01)) {
    if (missing(distr)) distr=cq(sim,d2post,normpost)
    else distr=match.arg(distr,several.ok=TRUE);
    if (('meta' %in% distr)&&(missing(mean.prior)||missing(sd.prior)))
      stop("distr contains 'meta' but 'mean.prior' or 'sd.prior' arguments missing");
    if (('bayes' %in% distr)&&missing(prior))
      stop("distr contains 'bayes' but prior argument missing");
    if (('prior' %in% distr)&&(missing(qprior)||missing(rprior)))
      stop("distr contains 'prior' but 'qprior' or 'rprior' argument missing");
    sim.byn=split(sim,sim$n);
    qq=do.call(rbind,lapply(n,function(n) {
      sim=sim.byn[[as.character(n)]];
      d.pop=sim$d.pop;
      m=length(d.pop);
      if ('meta' %in% distr) meta=meta_normd2t(mean.prior,sd.prior,d0,n);
      if ('bayes' %in% distr) init_bayes(n=n,d0=d0,prior=prior);
      qq=do.call(rbind,lapply(distr,function(distr) {
        qq=switch(
          distr,
          sim=quantile(d.pop,probs=p),
          d2tpost=q_d2t(n=n,p=p,d0=d0),
          normpost=qnorm(p=p,mean=d0,sd=sqrt(2/n)),
          meta=qnorm(p=p,mean=meta$mean,sd=meta$sd),
          bayes=q_bayes(p),
          prior=qprior(p),
          stop(paste0('Unknown distribution: ',distr,'. should have been caught earlier!')));
        rq=switch(
          distr,
          sim=quantile(r_sim(m,d.pop),probs=p),
          d2tpost=quantile(r_d2t(m,n=n,d0=d0),probs=p),
          normpost=quantile(rnorm(m,mean=d0,sd=sqrt(2/n)),probs=p),
          meta=quantile(rnorm(m,mean=meta$mean,sd=meta$sd),probs=p),
          bayes=quantile(r_bayes(m),probs=p),
          prior=quantile(rprior(m),probs=p),
          stop(paste0('Unknown distribution: ',distr,'. should have been caught earlier!')));
        data.frame(distr,qq,rq,stringsAsFactors=F,row.names=NULL);
      }));
      data.frame(n=n,p,qq,row.names=NULL);
    }));
    qq;
  }
