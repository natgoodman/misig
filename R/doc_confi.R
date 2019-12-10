#################################################################################
##
## Author:  Nat Goodman
## Created: 19-07-16
##          from confi.R created 19-07-04
## with additional code copied from confx.R 19-11-21
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate figures and tables for confi document
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
library(wesanderson);                   # for line colors
## --- Generate Figures and Tables for confi Blog Post ---
## no sections.
## TRN=T is supposed to create figures that work for Bob Reed's TRN Word format
##   but DIDN't do the right thing when I made the plots for the post. The scaling
##   messed up the in-plot text, and I didn't figure out how to create the right sized
##   png file for 1 row, 2 cols
##   Instead, we ended up living with regular plots (1 figure per file)
## NG 19-12-09: fixed code for TRN figures. not great but works...
doc_confi=
  function(sect=NULL,TRN=FALSE,n.fig=40,dobs.fig=0.5,dpop.fig=0.05,conf.fig=0.95) {
    if (TRN) figdev_start('samp_distr',nrow=1);
    figblk_start();
    dofig(plot_null,'null',n=n.fig,d.obs=dobs.fig);
    dofig(plot_mult,'mult',n=n.fig,d.obs=dobs.fig);
    figblk_end();
    if (TRN) figdev_start('plaus_confi',nrow=1);
    figblk_start();
    dofig(plot_plaus,'plaus',d.obs=dobs.fig);
    dofig(plot_confi,'confi',d.obs=dobs.fig);
    figblk_end();
    if (TRN) figdev_start('cicvr',nrow=1,ncol=1);
    dofig(plot_cicvr,'cicvr',d.pop=dpop.fig,conf.level=conf.fig);
    figdev_end();
    invisible();
  }
## Figure 1a. hist and sampling distribution under the null
plot_null=
  function(n=40,d.pop=0,d.obs=0.5,
           xlab='d.obs',ylab='probability density',title=NULL,
           col.hist='grey90',border.hist='grey80',legend.x0=-0.75,dinc=1e-3) {
   if (is.null(title)) title=figtitle('Histogram and sampling distribution under the null',n=n);
   sim=get_sim_fixd(n=n,d=d.pop);
    ## draw the histogram
    plothist(sim,col=col.hist,border=border.hist,title=title,xlab=xlab,ylab=ylab,
             legend.x0=legend.x0);
    ## add the sampling distribution
    plotpvsd(n=n,d0=d.pop,vline=d.obs,dinc=dinc,vhcol='grey10',vlab=FALSE,add=T);
  }
## Figure 1b. multiple sampling distributions
plot_mult=
  function(n=40,d.pop=c(seq(-0.1,by=0.1,len=4),seq(1.1,by=-0.1,len=4)),d.obs=0.5,
           title=NULL,xlab='d',ylab='probability density',
           xlim=c(-0.5,1.5),lwd=1,
           dinc=1e-3,d=seq(min(xlim),max(xlim),by=dinc)) {
    if (length(d.pop)<1) return();
    if (is.null(title)) title=figtitle('Sampling distributions for several values of d.pop',n=n);
    ## plot 1st d.pop separately to set graphical params
    d0=d.pop[1];
    ymax=d_d2t(n=n,d0=d0,d=d0)+strheight(d0);
    plotpvsd(n=n,d0=d0,d=d,vline=d.obs,col=d2plauscol(n=n,d0=d0,d=d),lwd=lwd,
             title=title,xlab=xlab,ylab=ylab,legend.label='plaus',legend.xscale=0.05,
             vhcol='grey10',vlab=FALSE,
             xlim=xlim,ylim=c(0,ymax));
    ## plot rest in a lopp
    sapply(tail(d.pop,-1),
           function(d0) plotpvsd(n=n,d0=d0,d=d,col=d2plauscol(n=n,d0=d0,d=d),lwd=lwd,add=T));
    ## add labels to each curve
    sapply(d.pop,function(d0) textd2t(n=40,d0=d0,label=d0));
    return();
  }
## Figure 2a. plaus vs. d.pop for several values of n
plot_plaus=
  function(n=c(40,100,400),d.obs=0.5,howmany=c(3,1),
           xlab='d.pop',ylab='plausibility value (plaus-value)',title=NULL,
           xlim=c(-0.2,1.2),ylim=c(0,1),cex.label=0.5,
           col=setNames(wes_palette("Moonrise2",n=length(n),type='continuous'),n),
           lwd=3,smooth='none',
           dinc=1e-3,d=seq(min(xlim),max(xlim),by=dinc)) {
    if (is.null(title)) title=figtitle('Plausibility-value vs. d.pop',d.obs=d.obs);
    y=do.call(cbind,lapply(n,function(n) y=d2plaus(n=n,d=d.obs,d0=d)));
    ## col=setNames(RColorBrewer::brewer.pal(ncol(y),'Accent'),n);
    legend.labels=paste(sep='=','n',n);
    plotm(x=d,y=y,title=title,lwd=lwd,col=col,smooth=smooth,
          legend='left',legend.labels=legend.labels,
          hline=c(0,0.05,0.5),hlab=c(FALSE,TRUE,TRUE),
          xlab=xlab,ylab=ylab);
    ## draw ci boundaries of interest
    ci.n=head(n,howmany[1]); ci.col=rep(head(col,howmany[1]),each=2);
    ci=as.vector(ci_d2t(n=ci.n,d=d.obs,conf.level=0.95));
    vline(x=ci,y=0.05,col=ci.col,lty='dotted',lwd=2,cex=cex.label,text=round(ci,digits=2));
    ci.n=tail(n,howmany[2]); ci.col=rep(tail(col,howmany[2]),each=2);
    ci=as.vector(ci_d2t(n=ci.n,d=d.obs,conf.level=0.5));
    vline(x=ci,y=0.5,col=ci.col,lty='dotted',lwd=2,cex=cex.label,text=round(ci,digits=2));
  }
## Figure 2b. confidence vs. d.pop for several values of n
plot_confi=
  function(n=c(40,100,400),d.obs=0.5,howmany=c(3,1),
           xlab='d.pop',ylab='confidence',title=NULL,
           xlim=c(-0.2,1.2),ylim=c(0,1),cex.label=0.5,
           col=setNames(wes_palette("Moonrise2",n=length(n),type='continuous'),n),
           lwd=3,smooth='none',
           dinc=1e-3,d=seq(min(xlim),max(xlim),by=dinc)) {
    if (is.null(title)) title=figtitle('Confidence vs. d.pop',d.obs=d.obs);
    y=do.call(cbind,lapply(n,function(n) y=1-d2plaus(n=n,d=d.obs,d0=d)));
    ## col=setNames(RColorBrewer::brewer.pal(ncol(y),'Set1'),n);
    legend.labels=paste(sep='=','n',n);
    plotm(x=d,y=y,title=title,lwd=lwd,col=col,smooth=smooth,
          legend='left',legend.labels=legend.labels,
          hline=c(0,0.5,0.95),hlab=c(FALSE,TRUE,TRUE),
          xlab=xlab,ylab=ylab);
    ## draw ci boundaries of interest
  ci.n=head(n,howmany[1]); ci.col=rep(head(col,howmany[1]),each=2);
    ci=as.vector(ci_d2t(n=ci.n,d=d.obs,conf.level=0.95));
    vline(x=ci,y=0.95,col=ci.col,lty='dotted',lwd=2,cex=cex.label,text=round(ci,digits=2));
    ci.n=tail(n,howmany[2]); ci.col=rep(tail(col,howmany[2]),each=2);
    ci=as.vector(ci_d2t(n=ci.n,d=d.obs,conf.level=0.5));
    vline(x=ci,y=0.5,col=ci.col,lty='dotted',lwd=2,cex=cex.label,text=round(ci,digits=2));
  }
## Figure 3. sampling distribution plus CIs to show coverage property
plot_cicvr=
  function(n=40,d.pop=0.05,conf.level=0.95,
           p=c(seq(0.025,0.975,by=0.025),seq(0.005,0.02,by=0.005),seq(0.98,0.995,by=0.005)),
           ## p=seq(0.01,0.99,by=0.01),
           d.ci=q_d2t(n=n,d0=d.pop,p=p),
           xlab='d.obs',ylab='probability density',title=NULL,vline='d.pop',
           xlim=d.pop+c(-1,1),lwd=2,
           ci.col=plaus2col(c(0.1,0.012)),ci.lty=cq(solid,dotted),
           ## ci.lwd=c(0.5,1),
           ci.lwd=c(1.5,2),
           ci.pch=19,ci.cex=c(0.5,0.75),ci.tol=1e-5,
           dinc=1e-3,d=seq(min(xlim),max(xlim),by=dinc)) {
    if (is.null(title)) title=figtitle('Coverage of confidence intervals',n=n,d.pop=d.pop);
    ci.col=rep(ci.col,len=2); ci.lty=rep(ci.lty,len=2); ci.lwd=rep(ci.lwd,len=2);
    ci.pch=rep(ci.pch,len=2); ci.cex=rep(ci.cex,len=2);
    plotpvsd(n=n,d0=d.pop,d=d,col=d2plauscol(n=n,d0=d.pop,d=d),lwd=lwd,
             title=title,xlab=xlab,ylab=ylab,legend.label='plaus',legend.xscale=0.075,
             xlim=xlim);
    ci=ci_d2t(n=n,d=d.ci,conf.level=conf.level);
    y=d_d2t(n=n,d=d.ci,d0=d.pop);
    col=ifelse(between(d.pop,ci[1,],ci[2,],ci.tol),ci.col[1],ci.col[2]);
    lty=ifelse(between(d.pop,ci[1,],ci[2,],ci.tol),ci.lty[1],ci.lty[2]);
    lwd=ifelse(between(d.pop,ci[1,],ci[2,],ci.tol),ci.lwd[1],ci.lwd[2]);
    pch=ifelse(between(d.pop,ci[1,],ci[2,],ci.tol),ci.pch[1],ci.pch[2]);
    cex=ifelse(between(d.pop,ci[1,],ci[2,],ci.tol),ci.cex[1],ci.cex[2]);
    hline(y=y,x0=ci[1,],x=ci[2,],col=col,lty=lty,lwd=lwd);
    points(d.ci,y,col=col,pch=pch,cex=cex);
    if ('d.pop' %in% vline) vhline(vline=d.pop,col='black',lty='solid');
    if ('d.conf' %in% vline) {
      d.conf=d_conf(n=n,d0=d.pop,conf.level=conf.level);
      vhline(vline=d.conf,col=ci.col[1],lty='dotted');
    }
  }

## d to plaus. adapted from d2pval
d2plaus=function(n,d,d0) {
  p=p_d2t(n=n,d=d,d0=d0);
  ## use data frame to extend all args to same length. ifelseneeds this
  cases=data.frame(d,d0,p);
  with(cases,ifelse(d<=d0,2*p,2*(1-p)));
}
## plaus to color. synonym for pval2col
plaus2col=pval2col
## d to color via plaus
d2plauscol=function(n,d0,d) plaus2col(d2plaus(n,d0,d));

## helper function for placing text on sampling distributions
textd2t=function(n,d0=0,label,where=cq(above,below,left,right),cex=0.7,col='grey10') {
  where=match.arg(where);
  y=d_d2t(n,d0=d0,d=d0);
  cex=cex*par('cex');
  label.bump=strheight(label,cex=cex)*0.5
  if (where=='above') {
    y=y+label.bump;
    adj=c(0.5,0);
  } else if (where=='below') {
    y=y-label.bump;
    adj=c(0.5,1);
  } else if (where=='left') {
    label=paste(sep='',label,' ');
    adj=c(1,0);
  } else {
    label=paste(sep='',' ',label);
    adj=c(0,0);
  }
  text(d0,y,label,cex=cex,col=col,adj=adj);
}

