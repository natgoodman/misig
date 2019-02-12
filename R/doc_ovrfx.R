#################################################################################
##
## Author:  Nat Goodman
## Created: 19-02-03 
##          from doc_siglo.R created 19-01-10
##          from repwr/R/doc_resig.R created 18-09-05
## Includes code from repwr/R/docfun_resig.R created 18-10-25
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate figures and tables for ovrfx blog post
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate Figures and Tables for ovrfx Blog Post ---
## no sections. only 2 figures
## n.fig is sample size for which figures plotted
doc_ovrfx=function(n.fig1=20,d.fig1=0.3,sect=parent(sect,NULL)) {
  param(n.meand,d0.meand);
  meand=get_meand_empi();
  meand.byd=split(meand,meand$d0);
  ## make splines that interpolate meand
  ##   empirically determined smooth.spline as interp function, and spar=0.3
  spar=0.3;
  meand.fit=lapply(meand.byd,function(meand) with(meand,smooth.spline(n,meand,spar=spar)))
  over.fit=lapply(meand.byd,function(meand) with(meand,smooth.spline(n,over,spar=spar)))
  ## support statements in text or drawn on figures
  ## do it upfront 'cuz some used in figures
  ## d_crit(n=20)=0.64
  d.crit1=d_crit(n.fig1); 
  ## mean observed significant effect size and overestimate for n=20
  ##   d=0.3, mean=0.81, over=2.7x
  ##   d=0.5, mean=0.86, over=1.7x
  ##   d=0.7, mean=0.93, over=1.3x
  meand20=sapply(meand.fit,function(meand.fit) predict(meand.fit,20)$y);
  over20=sapply(over.fit,function(over.fit) predict(over.fit,20)$y);
  ## n that achieves over=1.25x
  ##   d=0.3 n=121
  ##   d=0.5 n=47
  ##   d=0.7 n=26
  nover=sapply(over.fit,function(over.fit) {
    n2over=function(n) predict(over.fit,n)$y;
    uniroot(function(n) n2over(n)-1.25,interval=c(10,200))$root});
  ## save supporting values in table
  names(meand20)=paste(sep='_','meand20',names(meand20));
  names(over20)=paste(sep='_','over20',names(over20));
  names(nover)=paste(sep='_','nover',names(nover));
  support=data.frame(d.crit1,t(meand20),t(over20),t(nover));
  dotbl(support);

  ## draw the figures
  ## figure 1
  title=title_ovrfx('P-values improve as observed effect size grows more extreme',n=n.fig1);
  x='d.sdz'; y='d.pop';
  sim=get_sim_rand(n=n.fig1);
  dofig(plotdvsd,'big_picture',sim=sim,x=x,y=y,vline=c(-d.crit1,d.crit1),xlim=c(-2,2),
        title=title,cex.main=1);
  ## figure 2
  title=title_ovrfx('All significant results are beyond 0.3',d=d.fig1,n=n.fig1);
  sim.zoom=subset(sim,subset=near(d.pop,d.fig1,tol=0.02));
  dofig(plotdvsd,'zoom',sim=sim.zoom,x=x,y=y,vline=c(-d.crit1,d.fig1,d.crit1),hline=d.fig1,hlab=F,
        legend=F,title=title,cex.main=1);
  ## figure 3
  title=title_ovrfx('Histogram of observed effect size',d=d.fig1,n=n.fig1);
  sim=get_sim_fixd(n=n.fig1,d=d.fig1);
  dofig(plothist,'hist',sim=sim,vline=c(-d.crit1,d.fig1,d.crit1),legend.x0=-1.15,xlim=c(-1.2,1.5),
        title=title,cex.main=1);
  ## figure 4
  title=title_ovrfx('Average observed effect size improves as n increases');
  x=seq(min(n.meand),max(n.meand),by=1);
  y=do.call(cbind,lapply(meand.fit,function(meand.fit) predict(meand.fit,x)$y));
  col=setNames(RColorBrewer::brewer.pal(ncol(y),'Set1'),names(meand.byd));
  legend.labels=
    c(paste(sep=' ','mean sig effect for true effect',names(meand.byd)),
      paste(sep=' ','n achieving 1.25x for true effect',names(meand.byd)));
  col=c(col,col);
  lty=c(rep('solid',len=3),rep('dotted',3));
  lwd=2;
  dofig(plotmeand,'meand',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        meand20=meand20,nover=nover,ylim=c(0.3,meand20[3]),
        xlab='sample size',ylab='effect size');
  invisible();
}
## generate title for doc_ovrfx
title_ovrfx=function(desc=NULL,d=NULL,n=NULL) {
  fig=paste(sep='','Figure ',figlabel());
  paste(collapse="\n",
        c(fig,
          paste(sep='. ',desc,
                paste(collapse=', ',c(if(!is.null(d))paste_nv(d),if(!is.null(n))paste_nv(n))))));
}

## plot meand vs n, d (figure 4)
plotmeand=function(x,y,title,col='black',lty='solid',lwd=1,legend.labels,
                 meand20=NULL,nover=NULL,...) {
  param(n.meand,d0.meand);
  ## draw the main plot
  plotm(x=x,y=y,col=col,lty=lty,lwd=lwd,title=title,legend.labels=legend.labels,smooth=F,
        ...);
  ## horizontal grid-like lines for d.pop
  abline(h=d0.meand,col=col,lty='dotted',lwd=1);
  ## horizontal lines for averages with n=20
  ## if (!is.null(meand20)) {
  ##   hline(x=20,y=meand20,col=col,lty='dotted',lwd=0.5,text=round(meand20,digits=2));
  ## }
  ## vertical line for n=20
  vline(x=20,y0=0,y=1,col='grey',lty='dotted',lwd=1,text=20);
  ## horizontal & vertical lines for overestimates
  if (!is.null(nover)) {
    hline(x=nover,y=1.25*d0.meand,col=col,lty='dotted',lwd=2,text=paste(sep='','1.25x',d0.meand));
    vline(x=nover,y=1.25*d0.meand,col=col,lty='dotted',lwd=2,text=round(nover));
  }
}
