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
doc_siglo=function(n.fig1=20,d.fig1=0.3,sect=parent(sect,NULL)) {
  ## support statements in text or drawn on figures
  ## do it upfront 'cuz some used in figures
  ##   d_crit(n=20)=0.64
  d.crit1=d_crit(n.fig1);
  ## compute/get data for plotting critical and mean observed effect sizes
  meand=get_data(meand);
  meand=subset(meand,subset=d0%in%c(0.3,0.5,0.7));
  ## TEMPORARY until I fix domeand
  meand=subset(meand,subset=n>=20);
  ##
  n=sort(unique(meand$n));
  d.crit=d_crit(n);
  ##  mean observed significant effect size
  d.mean3=meand[meand$d0==0.3,];
  d.mean5=meand[meand$d0==0.5,];
  ##  n that achieves d.crit=d.pop
  n.crit3=uniroot(function(n) d2pval(n=n,d=0.3)-param(sig.level),interval=c(10,1000))$root;
  n.crit5=uniroot(function(n) d2pval(n=n,d=0.5)-param(sig.level),interval=c(10,1000))$root;
  ##  need n=122 for d=0.3 to achieve over=1.25, n=45 for d=0.5
  n.over3=n2over(meand,d0=0.3);
  n.over5=n2over(meand,d0=0.5);
  ## save supporting values in table
  support=data.frame(
    ##   d_crit(n=20)=0.64. computed above
    d.crit=d.crit1,
    ##   for n=20, d=0.3, average=0.82, overestimate of 2.7x
    ##   for n=20, d=0.5, average=0.87, overestimate of 1.7x
    d.mean3=d.mean3[d.mean3$n==n.fig1,'meand'],
    d.mean5=d.mean5[d.mean5$n==n.fig1,'meand'],
    d.over3=d.mean3[d.mean3$n==n.fig1,'over'],
    d.over5=d.mean5[d.mean5$n==n.fig1,'over'],
    ##  to achieve 1.25x, for d=0.3, n=122, for d=0.5, n=46
    n.over3=n.over3,
    n.over5=n.over5,
    ##   need n=87 for d=0.3 to be significant, n=32 for d=0.5
    n.crit3=n.crit3,
    n.crit5=n.crit5);
  dotbl(support);

  ## draw the figures
  ## figure 1
  title=title_siglo('P-values improve as observed effect size grows more extreme');
  x='d.sdz'; y='d.pop';
  sim=get_sim(n=n.fig1);
  sim.pos_dsdz=subset(sim,subset=d.sdz>=0);
  dofig(plotdvsd,'big_picture',sim=sim.pos_dsdz,x=x,y=y,d.crit=d.crit1,d.pop=d.fig1,xlim=c(0,2),
        title=title);
  ## figure 2
  title=title_siglo('Sharp boundary between nonsignificant and significant p-values');
  sim.zoom_in=subset(sim,subset=near(d.sdz,d.crit1,0.05));
  dofig(plotdvsd,'zoom_in',sim=sim.zoom_in,x=x,y=y,d.crit=d.crit1,d.pop=d.fig1,legend=F,
        vline=d.crit1,vlab=F,title=title);
  ## figure 3
  title=title_siglo('Critical and average observed effect sizes improve as n increases',n=NULL);
  ## y=data.frame(
  ##   d.crit,d.mean3=d.mean3[,'meand'],d.mean5=d.mean5[,'meand'],d.mean7=d.mean7[,'meand']);
  y=data.frame(
    d.crit,d.mean3=d.mean3[,'meand'],d.mean5=d.mean5[,'meand']);
  col=setNames(RColorBrewer::brewer.pal(ncol(y),'Set1'),cq(d.crit,0.3,0.5))
  legend.labels=
    c('critical effect size',
      'mean significant effect size for d.pop=0.3','mean significant effect size for d.pop=0.5');
  ## dofig(plotm,'dcrit_dmean',x=n,y=y,title=title,cex.title=1,lwd=2,col=col,legend='topright',
  ##       legend.labels=
  ##         c('critical d','mean d for d.pop=0.3','mean d for d.pop=0.5','mean d for d.pop=0.7'),
  ##       xlab='sample size',ylab='effect size',
  ##       vline=round(c(n.fig1,n.crit3,n.crit5)),hline=c(0.3,0.5,0.7));
  dofig(plotdvsn,'dcrit_dmean',x=n,y=y,title=title,lwd=2,col=col,
        legend='topright',legend.labels=legend.labels,
        n.crit=c(n.crit3,n.crit5),n.over=c(n.over3,n.over5),
        xlab='sample size',ylab='effect size');
  invisible();
}
## generate title for doc_siglo
title_siglo=function(desc=NULL,n=parent(n.fig1)) {
  fig=paste(sep='','Figure ',figlabel());
  paste(collapse="\n",c(fig,desc,if(!is.null(n)) paste_nv(n)));
}
## compute n that achieves desired overestimate for given d0 by interpolating meand
n2over=function(meand,d0,over=1.25) {
  md=meand[meand$d0==d0,cq(n,over)];
  ## empirically determined smooth.spline as interp function, and spar=0.1
  n2over=function(n) predict(with(md,smooth.spline(n,over,spar=0.1)),n)$y;
  uniroot(function(n) n2over(n)-1.25,interval=c(10,200))$root;
}
## plot d.crit, d.mean vs n (figure 3)
plotdvsn=function(x,y,title,col='black',lty='solid',lwd=1,legend.labels,
                  d.pop=c(0.3,0.5),n.crit=NULL,n.over=NULL,...) {
  col.pop=col[as.character(d.pop)];
  col.crit=col['d.crit'];
  ## add legend text for extra lines and adjust line properties
  legend.labels=c(legend.labels,
                  'minimum significant sample size',
                  'sample size that achieves 1.25x 0.3',
                  'sample size that achieves 1.25x 0.5');
  col=c(rep(col,len=3),rep(col,len=3));
  lty=c(rep(lty,len=3),rep('dotted',3));
  lwd=c(rep(lwd,len=3),rep(2,3));
  ## draw the main plot
  plotm(x=x,y=y,col=col,lty=lty,lwd=lwd,title=title,legend.labels=legend.labels,...);
  ## horizontal grid-like lines for d.pop
  abline(h=d.pop,col=col.pop,lty='dotted',lwd=1);
  ## horizontal & vertical lines for overestimates
  if (!is.null(n.over)) {
    hline(x=n.over,y=1.25*d.pop,col=col.pop,lty='dotted',lwd=2,text=paste(sep='','1.25x',d.pop));
    vline(x=n.over,y=1.25*d.pop,col=col.pop,lty='dotted',lwd=2,text=round(n.over));
  }
  ## vertical lines for critical crossovers
  if (!is.null(n.crit)) {
    vline(x=n.crit,y=d.pop,col=col.crit,lty='dotted',lwd=2,text=round(n.crit));
  }
}
