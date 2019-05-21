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
## no sections. only 4 figures
## n.fig is sample size for which figures plotted
doc_ovrfx=function(n.fig1=20,d.fig1=0.3,sect=parent(sect,NULL)) {
  param(n.fixd,d.fixd);
  meand.theo=get_data(meand.d2t);
  meand.t.byd=split(meand.theo,meand.theo$d0);
  ## make splines that interpolate meand, over
  ##   empirically determined smooth.spline as interp function, and spar=0.3
  spar=0.3;
  meand.t.fit=lapply(meand.t.byd,function(meand) with(meand,smooth.spline(n,meand,spar=spar)));
  over.t.fit=lapply(meand.t.byd,function(meand) with(meand,smooth.spline(n,over,spar=spar)));
  ## support statements in text or drawn on figures
  ## do it upfront 'cuz some used in figures
  ## d_crit(n=20)=0.64
  d.crit1=d_crit(n.fig1); 
  ## mean observed significant effect size and overestimate for n=20
  ##   d=0.3, mean=0.81, over=2.7x
  ##   d=0.5, mean=0.86, over=1.7x
  ##   d=0.7, mean=0.93, over=1.3x
  meand20=sapply(meand.t.fit,function(fit) predict(fit,20)$y);
  over20=sapply(over.t.fit,function(fit) predict(fit,20)$y);
  ## n that achieves over=1.25x
  ##   d=0.3 n=122
  ##   d=0.5 n=47
  ##   d=0.7 n=26
  nover=sapply(over.t.fit,function(fit) {
    n2over=function(n) predict(fit,n)$y;
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
  title=title_ovrfx('Histogram of observed effect size',d=d.fig1,n=n.fig1);
  sim=get_sim_fixd(n=n.fig1,d=d.fig1);
  ylim=c(0,d_d2t(n=100,d=0.3,d0=0.3));  # set ylim to match figure 3
  dofig(plothist,'hist',sim=sim,vline=c(-d.crit1,d.fig1,d.crit1),title=title,cex.main=1,
        xlim=c(-2,2),ylim=ylim);
  ## dofig(plothist,'hist',sim=sim,vline=c(-d.crit1,d.fig1,d.crit1),title=title,cex.main=1,
  ##       legend.x0=-1.85,xlim=c(-2,2),ylim=c(0, d_d2t(n=100,d=0.3,d0=0.3)));
  ## figure 3
  title=title_ovrfx('Sampling distributions show impact of increasing n and d');
  dofig(plotsmpldist,'smpl_dist',vline=c(-d.crit1,d.fig1,d.crit1),title=title);
  ## figure 4
  title=title_ovrfx('Average observed effect size improves as n increases');
  x=seq(min(n.fixd),max(n.fixd),by=1);
  y=do.call(cbind,lapply(meand.t.fit,function(fit) predict(fit,x)$y));
  col=setNames(RColorBrewer::brewer.pal(ncol(y),'Set1'),names(meand.t.byd));
  legend.labels=
    c(paste(sep=' ','mean sig effect for true effect',names(meand.t.byd)),
      paste(sep=' ','n achieving 1.25x for true effect',names(meand.t.byd)));
  col=rep(col,2);
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
## plot sampling distributions (figure 3)
## mostly hard-coded because no easy way to automate placement of text
plotsmpldist=function(vline=NULL,title,dlim=c(-2,2)) {
  ## do n=100 first because it's way taller than the others
  plotpvsd(n=100,d0=0.3,dlim=dlim,add=F,vline=vline,title=title);
  plotpvsd(n=20,d0=0.3,dlim=dlim,add=T);
  plotpvsd(n=20,d0=0.7,dlim=dlim,add=T);
  dtext(n=20,d0=0.3,y=0.5,side='left');
  dtext(n=20,d0=0.7,y=0.5,side='right');
  dtext(n=100,d0=0.3,y=2.5,side='left');
}
## helper function for placing text on sampling distributions (figure 3)
dtext=function(n,d0,y,label=paste(sep=', ',paste_nv('d',d0),paste_nv('n',n)),side=cq(left,right),
               cex=0.9,col='grey10') {
  side=match.arg(side);
  ## invert d_d2t to find d corresponding to y
  y2d=function(d) d_d2t(n=n,d=d,d0=d0)-y;
  if(side=='left') {
    d=uniroot(y2d,interval=c(-10,d0))$root
    label=paste(sep='',label,' ');
    adj=c(1,0);
  } else {
    d=uniroot(y2d,interval=c(d0,10))$root;
    label=paste(sep='',' ',label);
    adj=c(0,0);
  }
  text(d,y,label,cex=cex,col=col,adj=adj);
}

## plot meand vs n, d (figure 4)
plotmeand=function(x,y,title,col='black',lty='solid',lwd=1,legend.labels,
                 meand20=NULL,nover=NULL,...) {
  param(n.fixd,d.fixd);
  ## draw the main plot
  plotm(x=x,y=y,col=col,lty=lty,lwd=lwd,title=title,legend.labels=legend.labels,smooth=F,
        ...);
  ## horizontal grid-like lines for d.pop
  abline(h=d.fixd,col=col,lty='dotted',lwd=1);
  ## horizontal lines for averages with n=20
  ## if (!is.null(meand20)) {
  ##   hline(x=20,y=meand20,col=col,lty='dotted',lwd=0.5,text=round(meand20,digits=2));
  ## }
  ## vertical line for n=20
  vline(x=20,y0=0,y=1,col='grey',lty='dotted',lwd=1,text=20);
  ## horizontal & vertical lines for overestimates
  if (!is.null(nover)) {
    hline(x=nover,y=1.25*d.fixd,col=col,lty='dotted',lwd=2,text=paste(sep='','1.25x',d.fixd));
    vline(x=nover,y=1.25*d.fixd,col=col,lty='dotted',lwd=2,text=round(nover));
  }
}
