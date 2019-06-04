#################################################################################
##
## Author:  Nat Goodman
## Created: 19-02-19
##          from doc_ovrfx.R created 19-02-03 
##          from doc_siglo.R created 19-01-10
##          from repwr/R/doc_resig.R created 18-09-05
## Includes code from repwr/R/docfun_resig.R created 18-10-25
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate figures and tables for ovrht blog post
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate Figures and Tables for ovrht Blog Post ---

## no sections. only 4 figures
## n.fig is sample size for which figures plotted
doc_ovrht=function(sect=parent(sect,NULL)) {
  ## support statements in text or drawn on figures
  pval=get_data(pval.d2ht);
  ci=get_data(ci.d2ht);
  support=list(
    ## for n=200, sd.het=0.2: htpval=0.38, ci_d2ht=+-0.44. ci_d2t=+-0.20
    pval1=subset(pval,subset=(n==200&sd.het==0.2)),
    ci1=subset(ci,subset=(n==200&sd.het%in%c(0.2,0))),
    ## for *n=200*, sd.het=0.05: pval inflation 1.59, ci inflation 1.12
    ## for *n=200*, sd.het=0.4: pval inflation 12.68, ci inflation 4.12
    pval2=subset(pval,subset=(n==200&sd.het%in%c(0.05,0.4))),
    ci2=subset(ci,subset=(n==200&sd.het%in%c(0.05,0.4))),
    ## for sd.het=0.05, n=c(20,100,200,400)
    ## pval inflation 1.05,1.29,1.59,2.19; ci inflation 1.01,1.06,1.12,1.22
    pval3=subset(pval,subset=(n%in%c(20,100,200,400)&sd.het==0.05)),
    ci3=subset(ci,subset=(n%in%c(20,100,200,400)&sd.het==0.05)),
    ## for sdhet=0.2, n=c(20,400)
    ## pval inflation 1.90,10.26; ci inflation 1.18,3.00
    pval4=subset(pval,subset=(n%in%c(20,400)&sd.het==0.20)),
    ci4=subset(ci,subset=(n%in%c(20,400)&sd.het==0.20)),
    end=NA);                          # placeholder for last entry
  support=lapply(support,function(x) round(x,digits=2));
  dotbl(support,obj.ok=T);
  ## draw the figures
  ## figure 1
  n.fig=200; d.fig=0; sd.fig=0.2;
  title=figtitle('Het histogram and conventional sampling distribution under the null',
                 sd.het=sd.fig,n=n.fig);
  sim=get_sim_hetd(n=n.fig,d=d.fig,sd=sd.fig);
  dofig(plothist_d2t,'hist_d2t',sim=sim,title=title,d=d.fig,n=n.fig);
  ## figure 2
  title=figtitle('P-value inflation worsens as sd.het and n increase');
  dofig(plotpval_over,'pval_inflation',pval=pval,title=title);
  ## figure 3
  title=figtitle('Confidence interval inflation worsens as sd.het and n increase');
  dofig(plotci_over,'ci_inflation',ci=ci,title=title);
  ## figure 4
  title=figtitle('Conventional and het sampling distributions for two values of n');
  n.fig=c(20,200); d.fig=0; sd.fig=0.2;
  dofig(plotd2t_d2ht,'d2t_d2ht',title=title,n=c(20,200),d=d.fig,sd.het=0.2);
    
  invisible();
}
## figure 1. histogram and d_d2t distribution
plothist_d2t=
  function(sim,title,n,d,xlim=d+c(-1,1),ylim=c(0,d_d2t(n=n,d=0)),
           col.hist='grey90',border.hist='grey80',
           ci.col='black',ci.lty='dotted',ci.lwd=2,
           ...) {
    plothist(sim=sim,col=col.hist,border=border.hist,title=title,
             xlim=xlim,ylim=ylim,ylab='probability density');
    d.crit=d_crit(n);
    plotpvsd(n=n,add=T,dinc=1e-3,vline=c(-d.crit,d.crit),...)
    ## plotpvsd(n=n,add=T,dinc=1e-3,...);
    ## draw ci
    x=ci_d2t(n=n,d=d);
    y=d_d2t(n=n,d=x);
    segments(x[1],y[1],x[2],y[2],col=ci.col,lty=ci.lty,lwd=ci.lwd);
    invisible();
}
## figure 2. p-value inflation
plotpval_over=function(pval,title,lwd=2) {
    pval=subset(pval,subset=(sd.het!=0));
    pval.bysd=split(pval,pval$sd.het);
    n=unique(pval$n);
    y=do.call(cbind,lapply(pval.bysd,function(pval) pval$over.tval));
    col=setNames(RColorBrewer::brewer.pal(ncol(y),'Set1'),names(pval.bysd));
    legend.labels=paste(sep='=','sd.het',names(pval.bysd));
    plotm(x=n,y=y,title=title,lwd=lwd,col=col,smooth='spline',
          legend='topleft',legend.labels=legend.labels,
          xlab='sample size',ylab='inflation (ratio of correct to conventional value)');
}
## figure 3. confidence interval inflation
plotci_over=function(ci,title,lwd=2) {
    ci=subset(ci,subset=(sd.het!=0));
    ci.bysd=split(ci,ci$sd.het);
    n=unique(ci$n);
    y=do.call(cbind,lapply(ci.bysd,function(ci) ci$over));
    col=setNames(RColorBrewer::brewer.pal(ncol(y),'Set1'),names(ci.bysd));
    legend.labels=paste(sep='=','sd.het',names(ci.bysd));
    plotm(x=n,y=y,title=title,lwd=lwd,col=col,smooth='spline',
          legend='topleft',legend.labels=legend.labels,
          xlab='sample size',ylab='inflation (ratio of correct to conventional value)');
}
## figure 4. sampling distributions
## mostly hard-coded because no easy way to automate placement of text
plotd2t_d2ht=
  function(title,n,d,sd.het,xlim=d+c(-1,1),ylim=c(0,d_d2t(n=n[2],d=0)),col.d2ht='grey70') {
    plotpvsd(n=n[1],title=title,xlim=xlim,ylim=ylim,ylab='probability density');
    plotpvsd(n=n[2],add=T);
    plotpvsd(n=n[1],sd.het=sd.het,col=col.d2ht,add=T);
    plotpvsd(n=n[2],sd.het=sd.het,col=col.d2ht,add=T);
    textd2t_d2h(n=n[1],sd.het=0,d=0,label=paste_nv('n',n[1]),where='above');
    textd2t_d2h(n=n[2],sd.het=0,d=0,label=paste_nv('n',n[2]),where='below');
    textd2t_d2h(n=n[1],sd.het=0.2,d=0,label=paste_nv('n',n[1]),where='below');
    textd2t_d2h(n=n[2],sd.het=0.2,d=0,label=paste_nv('n',n[2]),where='above');
}

## helper function for placing text on sampling distributions (figure 4)
textd2t_d2h=function(n,sd.het,d=0,label,where=cq(above,below,left,right),cex=0.8,col='grey10') {
  where=match.arg(where);
  y=d_d2ht(n,sd.het=sd.het,d=0);
  if (where=='above') {
    y=y+strheight(label,cex=cex)*0.5;
    adj=c(0.5,0);
  } else if (where=='below') {
    y=y-strheight(label,cex=cex)*0.5;
    adj=c(0.5,1);
  } else if (where=='left') {
    label=paste(sep='',label,' ');
    adj=c(1,0);
  } else {
    label=paste(sep='',' ',label);
    adj=c(0,0);
  }
  text(d,y,label,cex=cex,col=col,adj=adj);
}
