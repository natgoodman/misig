#################################################################################
##
## Author:  Nat Goodman
## Created: 19-08-04
##          from doc_confi.R created 19-07-16
##          from confi.R created 19-07-04
##
## Copyright (C) 2019 Nat Goodman.
## 
## Specialized plot functions for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## wrapper for plotm to show coverage from simulated data
## ci is data.frame from doci_confi
## n, distr optional - if present, limit ci to these
## lty, lwd defaults are for single distribution
plotm_cover=
  function(ci,n,distr,title=NULL,smooth='spline',
           col=n2col(),lty='solid',lwd=2,
           legend='topleft',xlab='confidence level',ylab='coverage',
           col.diag='grey',lty.diag='dotted',lwd.diag=1.5,...) {
    ## limit ci to n, distr if given
    if (!missing(n)) ci=ci[ci$n%in%n,];
    if (!missing(distr)) ci=ci[ci$distr%in%distr,];
    if (nrow(ci)==0)
      stop(paste(
        collapse=' ',c("Nothing to plot: 'ci' empty",
                       if(!missing(n)||!missing(distr)) "after filtering on 'n' and 'distr'")));
    ## extract params from ci
    n=unique(ci$n);
    conf.level=unique(ci$conf.level);
    distr=unique(ci$distr);
    ## create matrix of y-values - each is a line
    cvr=subset(ci,select=c(n,conf.level,distr,cvr));
    cvr.by=split(cvr,paste(sep='.',cvr$n,cvr$distr));
    y=do.call(cbind,lapply(cvr.by,function(cvr) cvr$cvr));
    ## split colnames to extract params
    cname=strsplit(colnames(y),'.',fixed=T);
    ## construct line properties and legend labels and properties
    if (length(distr)>1) {
      if (missing(lty)) lty=distr2lty(distr);
      if (missing(lwd)) lwd=distr2lwd(distr);
      y2col=sapply(cname,function(cname) col[cname[1]]);
      y2lty=sapply(cname,function(cname) lty[cname[2]]);
      y2lwd=sapply(cname,function(cname) lwd[cname[2]]);
      legend.args=list(labels=c('colors depict n',paste_nv('n'),
                                'line types depict distribution',distr2label(distr)),
                       col=c('black',unique(y2col),'black',rep('grey20',length(distr))),
                       lty=c(NA,rep('solid',length(n)),NA,y2lty[distr]),
                       lwd=c(NA,rep(2,length(n)),NA,y2lwd[distr]),
                       seg.len=4);
    } else {
      y2col=sapply(cname,function(cname) col[cname[1]]);
      y2lty=lty; y2lwd=lwd;
      legend.args=list(labels=paste_nv('n'),col=col,lty=lty,lwd=2);
    }
    plotm(x=conf.level,y=y,title=title,col=y2col,lty=y2lty,lwd=y2lwd,smooth=smooth,
          legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,...);
    abline(a=0,b=1,col=col.diag,lty=lty.diag,lwd=lwd.diag);
    grid();
  }
## wrapper for plotm to plot concurves from theory and simulated data
## ci is data.frame from doci_confi
## n, distr optional - if present, limit ci to these
## lty, lwd defaults are for single distribution
plotm_concurve=
  function(ci,n,distr,title=NULL,smooth='spline',
           col=n2col(),lty='solid',lwd=2,
           legend='bottomleft',xlab='predicted d.pop interval',ylab='confidence level',...) {
    ## limit ci to n, distr if given
    if (!missing(n)) ci=ci[ci$n%in%n,];
    if (!missing(distr)) ci=ci[ci$distr%in%distr,];
    if (nrow(ci)==0)
      stop(paste(
        collapse=' ',c("Nothing to plot: 'ci' empty",
                       if(!missing(n)||!missing(distr)) "after filtering on 'n' and 'distr'")));
    ## extract params from ci
    n=unique(ci$n);
    conf.level=unique(ci$conf.level);
    distr=unique(ci$distr);
    ## create matrix of x-values - each is a line
    concrv=subset(ci,select=c(n,conf.level,distr,lo,hi));
    concrv.by=split(concrv,paste(sep='.',concrv$n,concrv$distr));
    x=do.call(cbind,lapply(concrv.by,function(concrv) concrv[,cq(lo,hi)]));
    ## split colnames to extract params
    cname=strsplit(colnames(x),'.',fixed=T);
    ## construct line properties and legend labels and properties
    if (length(distr)>1) {
      if (missing(lty)) lty=distr2lty(distr);
      if (missing(lwd)) lwd=distr2lwd(distr);
      x2col=sapply(cname,function(cname) col[cname[1]]);
      x2lty=sapply(cname,function(cname) lty[cname[2]]);
      x2lwd=sapply(cname,function(cname) lwd[cname[2]]);
      legend.args=list(labels=c('colors depict n',paste_nv('n'),
                                'line types depict distribution',distr2label(distr)),
                       col=c('black',unique(x2col),'black',rep('grey20',length(distr))),
                       lty=c(NA,rep('solid',length(n)),NA,x2lty[distr]),
                       lwd=c(NA,rep(2,length(n)),NA,x2lwd[distr]),
                       seg.len=4);
    } else {
      x2col=sapply(cname,function(cname) col[cname[1]]);
      x2lty=lty; x2lwd=lwd;
      legend.args=list(labels=paste_nv('n'),col=col,lty=lty,lwd=2);
    }
    plotm(x=x,y=conf.level,title=title,col=x2col,lty=x2lty,lwd=x2lwd,smooth=smooth,
          legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,...);
    grid();
  }
## plot hist of d.pop overlaid with various distributions
## TODO: decide whether to keep both norm and normpost
plothist_dpop=
  function(sim,n,d0,distr=cq(d2t,norm,meta,bayes,prior),meta,prior,
           title=NULL,cex.title='auto',legend='topright',
           xlab='d.pop',ylab='probability density',xlim=NULL,ylim=NULL,
           col.hist='grey90',border.hist='grey80',breaks='Sturges',
           col.distr=NULL,lwd.distr=2,lty.distr='solid',
           vhlty='dashed',vhcol='grey50',vhlwd=1,vlab=TRUE,vhdigits=2,
           ...){
    if (missing(distr)) distr=NULL else distr=match.arg(distr,several.ok=TRUE);
    if (('meta' %in% distr)&&missing(meta))
      stop("distr contains 'meta' but meta argument missing");
    if (('bayes' %in% distr)&&missing(prior))
      stop("distr contains 'bayes' but prior argument missing");
    sim=sim[sim$n==n,];
    ld=length(distr);
    col=if (missing(col.distr))
          RColorBrewer::brewer.pal(max(3,ld),'Set1')[1:ld] else rep(col.distr,len=ld);
    lwd=rep(lwd.distr,len=ld);
    lty=rep(lty.distr,len=ld);
    distr2col=setNames(col,distr);
    distr2lty=setNames(lty,distr);
    distr2lwd=setNames(lwd,distr);
    d.pop=sim$d.pop;
    init_bayes(n=n,d0=d0,prior=prior);
    hist.obj=hist(d.pop,breaks=breaks,plot=FALSE);
    if (is.null(xlim)) xlim=range(hist.obj$breaks); # default per R refman
    if (is.null(ylim)&&!is.null(distr)) {
      ## compute max y value across hist and all distributions
      ymax=max(hist.obj$density,do.call(max,lapply(distr,function(distr) {
        switch(distr,
               d2t=d_d2t(n=n,d0=d0,d=d0),
               norm=dnorm(d0,mean=d0,sd=sqrt(2/n)),
               meta=dnorm(meta$mean,mean=meta$mean,sd=meta$sd),
               bayes=optimize(d_bayes,c(-10,10),maximum=T)$objective,
               prior=optimize(prior,c(-10,10),maximum=T)$objective,
               ## none=0,
               stop(paste0('Unknown distribution: ',distr,'. should have been caught earlier!')));
      })));
      ylim=c(0,ymax);
    }
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    plot(hist.obj,freq=F,main=title,cex.main=cex.title,
         col=col.hist,border=border.hist,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,...);
    ## mean.data=mean(d.pop);
    ## vhline(vline=mean.data,vlab=vlab,vhdigits=vhdigits,lty=vhlty,col=vhcol,lwd=vhlwd);
    ## NG 19-09-03: REAL HACK. plot median_bayes. should be done under switch posiibly 
    ##              with option to do medians for other distributions
    vhline(vline=median_bayes(),vlab=vlab,vhdigits=vhdigits,lty=vhlty,col=vhcol,lwd=vhlwd);
    if (!is.null(distr)) {
      ## add the distributions to the plot
      sapply(distr,function(distr) {
        col=distr2col[distr];
        lty=distr2lty[distr];
        lwd=distr2lwd[distr];
        switch(distr,
               d2t=curve(d_d2t(n=n,d0=d0,d=x),col=col,lty=lty,lwd=lwd,add=T),
               norm=curve(dnorm(x,mean=d0,sd=sqrt(2/n)),col=col,lty=lty,lwd=lwd,add=T),
               meta=curve(dnorm(x,mean=meta$mean,sd=meta$sd),col=col,lty=lty,lwd=lwd,add=T),
               bayes=curve(d_bayes,col=col,lty=lty,lwd=lwd,add=T),
               prior=curve(prior,col=col,lty=lty,lwd=lwd,add=T),
               ## none=0,                  # do nothing
               stop(paste0('Unknown distribution: ',distr,'. should have been caught earlier!')));
      });
    }
    grid();
    ## add legend
    if (!is.null(distr)) {
      labels=distr2label(distr);
      col=distr2col[distr];
      lwd=distr2lwd[distr];
      lty=distr2lty[distr];
      legend(legend,legend=labels,title='distribution',col=col,lwd=lwd,lty=lty,cex=0.8,bty='n');
    }
  }
## TODO: decide on names for this vs plotqq_dpop
## wrapper for plotm to show quantiles
## qq is data.frame from doqq_confi
## n, distr optional - if present, limit qq to these
## lty, lwd defaults are for single distribution
plotm_qq=
  function(qq,n,distr,title=NULL,smooth='spline',
           col=n2col(),lty='solid',lwd=2,
           legend='topleft',xlab='quantile level (p)',ylab='d.pop', ...) {
    ## limit qq to n, distr if given
    if (!missing(n)) qq=qq[qq$n%in%n,];
    if (!missing(distr)) qq=qq[qq$distr%in%distr,];
    if (nrow(qq)==0)
      stop(paste(
        collapse=' ',c("Nothing to plot: 'qq' empty",
                       if(!missing(n)||!missing(distr)) "after filtering on 'n' and 'distr'")));
    ## extract params from qq
    n=unique(qq$n);
    p=unique(qq$p);
    distr=unique(qq$distr);
    ## create matrix of y-values - each is a line
    qq.by=split(qq,paste(sep='.',qq$n,qq$distr));
    y=do.call(cbind,lapply(qq.by,function(qq) qq$qq));
    ## split colnames to extract params
    cname=strsplit(colnames(y),'.',fixed=T);
    ## construct line properties and legend labels and properties
    if (length(distr)>1) {
      if (missing(lty)) lty=distr2lty(distr);
      if (missing(lwd)) lwd=distr2lwd(distr);
      y2col=sapply(cname,function(cname) col[cname[1]]);
      y2lty=sapply(cname,function(cname) lty[cname[2]]);
      y2lwd=sapply(cname,function(cname) lwd[cname[2]]);
      legend.args=list(labels=c('colors depict n',paste_nv('n'),
                                'line types depict distribution',distr2label(distr)),
                       col=c('black',unique(y2col),'black',rep('grey20',length(distr))),
                       lty=c(NA,rep('solid',length(n)),NA,y2lty[distr]),
                       lwd=c(NA,rep(2,length(n)),NA,y2lwd[distr]),
                       seg.len=4);
    } else {
      y2col=sapply(cname,function(cname) col[cname[1]]);
      y2lty=lty; y2lwd=lwd;
      legend.args=list(labels=paste_nv('n'),col=col,lty=lty,lwd=2);
    }
    plotm(x=p,y=y,title=title,col=y2col,lty=y2lty,lwd=y2lwd,smooth=smooth,
          legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,...);
    grid();
  }

## TODO: decide on names for this vs plotm_qq
## qq plot of (typically) d.pop vs. various distributions
## qq is data.frame from doqq_confi
## n (single value!), distr optional - if present, limit qq to these
## in any case, qq must contain single value of n
## x tells which column to use for 'x' axis quantiles
plotqq_dpop=
  function(qq,n,d0,distr,x='sim',rand,title=NULL,cex.title='auto',legend='topleft',
           xlab=paste('quantiles from',distr2label(x)),ylab='quantiles from distribution(s)',
           cex=0.5,pch=19,xlim=NULL,ylim=NULL,
           col.distr='black',col.rand='grey',col.diag='red',lty.diag='dotted',lwd.diag=1,...) {
    ## limit qq to n, distr if given
    if (!missing(n)) qq=qq[qq$n%in%n,];
    if (!missing(rand)) rq=qq[qq$distr%in%rand,] else rq=rq[0,];
    if (!missing(distr)) qq=qq[qq$distr%in%distr,];
    if (nrow(qq)==0&&nrow(rq)==0)
      stop(paste(
        collapse=' ',c("Nothing to plot: 'qq' and 'rq' empty",
                       if(!missing(n)||!missing(distr)||!missing(rand))
                         "after filtering on 'n', 'distr', 'rand'")));
    ## extract params from qq
    n=unique(qq$n);
    p=unique(qq$p);
    distr=unique(qq$distr);
    rand=unique(rq$distr);
    ## further error checking. make sure 'n' is single value and 'x' in distr
    if (length(n)>1) stop(paste("'n' (from 'qq' after filtering) must be single value, not",
                                paste(collapse=', ',n)))
    if (x%notin%distr) stop(paste("'x' distribution",x,"not in qq"));
    distr=setdiff(distr,x);
    ld=length(distr);
    lr=length(rand);
    if (ld+lr==0)
      stop("Nothing to plot: 'qq' (after filtering) has one distribution and rand is empty");
    ## create vector of x-values and matrices of y-values and y-rand values
    xx=qq[qq$distr==x,'qq'];
    qq=qq[qq$distr!=x,];
    qq.by=split(qq,qq$distr);
    yy=do.call(cbind,lapply(qq.by,function(qq) qq$qq));
    if (ld>1) yy=yy[,distr];            # reorder if multiple columns. NOT if single column
                                        #   else R turns it into a vector... sigh
    col.yy=if(ld>0)
             setNames(if (missing(col.distr)) RColorBrewer::brewer.pal(max(3,ld),'Set1')[1:ld]
                      else rep(col.distr,len=ld),distr);
    col=sapply(colnames(yy),function(cname) col.yy[cname]);
    labels=distr2label(distr);
    if (lr>0) {
      rq.by=split(rq,rq$distr);
      ry=do.call(cbind,lapply(rq.by,function(rq) rq$rq));
      if (lr>1) ry=ry[,rand];           # reorder if multiple columns. NOT if single column
                                        #   else R turns it into a vector... sigh
      yy=cbind(yy,ry);
      col.ry=setNames(if (missing(col.rand)) colorRampPalette(cq(gray30,gray70))(lr)
                      else rep(col.rand,len=lr),rand);
      col=c(col,sapply(colnames(ry),function(cname) col.ry[cname]))
      labels=c(labels,paste('random',distr2label(rand)));
    }
    legend.args=list(labels=c('colors depict distribution',labels),
                     col=c('black',col));
    if (is.null(xlim)) xlim=range(xx);
    if (is.null(ylim)) ylim=range(yy);
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    matplot(xx,yy,xlim=xlim,ylim=ylim,main=title,cex.main=cex.title,
            xlab=xlab,ylab=ylab,pch=pch,cex=cex,col=col);
    abline(a=0,b=1,col=col.diag,lty=lty.diag,lwd=lwd.diag);
    grid();
    ## add legend
    legend(legend,legend=labels,title='distribution',col=col,pch=pch,cex=0.8,bty='n');
  }

####################
## generate standard line properties & legend labels 
## n to color
n2col=function(n=param(n.confi))
  setNames(RColorBrewer::brewer.pal(max(3,length(n)),'Set1')[1:length(n)],n);
## distr to line type. by default, sim is 'dotted', rest are 'solid','dashed', 'dotdash', ...
distr2lty=
  function(distr,lty.sim='dotted',lty.all=cq(solid,dashed,dotdash,longdash,twodash)) {
  lty=NULL;
  if ('sim' %in% distr) {
    lty=c(sim=lty.sim);
    distr=setdiff(distr,'sim');
  }
  l=length(distr);
  if (l>length(lty.all)) stop('distr longer than line types. is distr right?');
  c(lty,setNames(lty.all[1:l],distr));
}
## distr to line width. by default, sim is 1.5, rest are 2
distr2lwd=function(distr,lwd.sim=1.5,lwd.all=c(2)) {
  lwd=NULL;
  if ('sim' %in% distr) {
    lwd=c(sim=lwd.sim);
    distr=setdiff(distr,'sim');
  }
  lwd.all=rep(lwd.all,len=length(distr));
  c(lwd,setNames(lwd.all,distr));
}
## distr to legend label
## distr.all defined in docfun_confi.R
distr2label=function(distr=distr.all) {
  c(std='standard',
    sim='simulated data',
    simq='simulated data (using q_sim)',
    d2tpost='d2t posterior',
    normpost='normal posterior',
    meta='meta-analytic posterior',
    bayes='bayesian posterior',
    d2t='d2t',norm='normal',prior='bayesian prior')[distr];
}
