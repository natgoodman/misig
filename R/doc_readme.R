#################################################################################
##
## Author:  Nat Goodman
## Created: 19-02-18
##          from doc_readme.R created 19-02-03 
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
## --- Generate Figures and Tables for misig README---
doc_readme=function(sect=NULL) {
  sect.all=cq(plotdvsd,plothist,plotpvsd,plotm,table);
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
  sapply(sect,function(sect) {
    ## compute section number. from stackoverflow.com/questions/5577727
    ## sectnum=which(sect==sect.all)[1];
    ## sect_start(sectnum); 
##### plotdvsd
    if (sect=='plotdvsd') {
      figblk_start();
      n=min(param(n.rand));
      d.crit=d_crit(n);
      sim=get_sim_rand(n=n);
      x='d.sdz'; y='d.pop';
      title=title_readme(c('Scatter plot',y,'vs',x),n=n);
      dofig(plotdvsd,NULL,sim=sim,x=x,y=y,vline=c(-d.crit,d.crit),xlim=c(-2,2),
            title=title,cex.main=1);
      sim.zoom=subset(sim,subset=near(d.sdz,d.crit,0.10));
      x='d.sdz'; y='d.pop';
      title=title_readme(c('Scatter plot.',y,'vs',x,'showing sharp boundary at',
                           round(d.crit,digits=2)),n=n);
      dofig(plotdvsd,'zoom',sim=sim.zoom,x=x,y=y,vline=d.crit,vlab=T,legend=F,
            title=title,cex.main=0.9);
    }
##### plothist
    if (sect=='plothist') {
      figblk_start();
      n=range(param(n.fixd));
      d=range(param(d.fixd));
      figname=cq(blue,red);
      sapply(1:2,function(i) {
        n.fig=n[i]; d.fig=d[i];
        d.crit=d_crit(n.fig);
        sim=get_sim_fixd(n=n.fig,d=d.fig);
        title=title_readme('Histogram of observed effect size',d.pop=d.fig,n=n.fig);
        dofig(plothist,figname[i],sim=sim,vline=c(-d.crit,d.fig,d.crit),title=title,cex.main=1,
              breaks=20);
      });
    }
##### plotpvsd 
    if (sect=='plotpvsd') {
      figblk_start();
      n=range(param(n.hetd));
      d=range(param(d.hetd));
      sd=range(param(sd.hetd));
      y=cq(density,cumulative);
      figname=substr(y,1,3);
      sapply(1:2,function(i) {
        n.fig=n[i]; d.fig=d[i]; sd.fig=sd[i]; y.fig=y[i]; 
        d.crit=d_htcrit(n.fig,sd.het=sd.fig);
        desc=paste(sep='','het distribution (',y.fig,')');
        title=title_readme(desc,d.het=d.fig,sd.het=sd.fig,n=n.fig);
        dofig(plotpvsd,figname[i],n=n.fig,d.het=d.fig,sd.het=sd.fig,y=y.fig,
              vline=c(-d.crit,d.crit),xlim=c(-1,2),title=title,cex.main=1);
      });
    }
##### plotm 
    if (sect=='plotm') {
      figblk_start();
      param(n.fixd,d.fixd);
      meand.simu=get_data(meand.fixd);
      meand.theo=get_data(meand.d2t);
      meand.s.byd=split(meand.simu,meand.simu$d0);
      meand.t.byd=split(meand.theo,meand.theo$d0);
      x=n.fixd;
      y=do.call(cbind,c(lapply(meand.s.byd,function(meand) meand$over),
                        lapply(meand.t.byd,function(meand) meand$over)));
      col=setNames(RColorBrewer::brewer.pal(length(d.fixd),'Set1'),d.fixd);
      legend.labels=
        c(paste(sep='=','meand.simu. d.pop',d.fixd),paste(sep='=','meand.theo. d.pop',d.fixd));
      col=rep(col,2);
      lty=c(rep('dotted',len=length(d.fixd)),rep('solid',length(d.fixd)));
      lwd=2;     
      title=title_readme('Line plot of mean significant observed effect size inflation (raw)');
      dofig(plotm,'meand',x=x,y=y,col=col,lty=lty,lwd=lwd,smooth=F,title=title,cex.main=0.9,
            legend.labels=legend.labels,legend='topright',
            xlab='sample size',ylab='inflation (ratio of actual to correct)');
      param(n.hetd,sd.hetd);
      pval.simu=get_data(pval.hetd);
      pval.theo=get_data(pval.d2ht);
      pval.s.bysd=split(pval.simu,pval.simu$sd.het);
      pval.t.bysd=split(pval.theo,pval.theo$sd.het);
      x=n.hetd;
      y=do.call(cbind,c(lapply(pval.s.bysd,function(pval) pval$over),
                        lapply(pval.t.bysd,function(pval) pval$over)));
      col=setNames(RColorBrewer::brewer.pal(length(sd.hetd),'Set1'),sd.hetd);
      legend.labels=
        c(paste(sep='=','pval.simu. sd.hetd',sd.hetd),paste(sep='=','pval.theo. sd.hetd',sd.hetd));
      col=rep(col,2);
      lty=c(rep('dotted',len=length(sd.hetd)),rep('solid',length(sd.hetd)));
      lwd=2;     
      title=title_readme('Line plot of pval inflation (smooth=aspline)');
      dofig(plotm,'pval',x=x,y=y,col=col,lty=lty,lwd=lwd,smooth=T,title=title,cex.main=1,
            legend.labels=legend.labels,legend='topleft',
            xlab='sample size',ylab='inflation (ratio of actual to correct)');  
    }
    if (sect=='table') {
      param(n.fixd,d.fixd);
      support=do.call(rbind,lapply(cq(meand.fixd,meand.d2t),function(what) {
        meand=get_data(list=what);
        meand.byd=split(meand,meand$d0);
        ## make splines that interpolate meand, over
        ##   empirically determined smooth.spline as interp function, and spar=0.3
        spar=0.3;
        meand.fit=lapply(meand.byd,function(meand) with(meand,smooth.spline(n,meand,spar=spar)));
        over.fit=lapply(meand.byd,function(meand) with(meand,smooth.spline(n,over,spar=spar)));
        meand20=sapply(meand.fit,function(fit) predict(fit,20)$y);
        over20=sapply(over.fit,function(fit) predict(fit,20)$y);
        nover=sapply(over.fit,function(fit) {
          n2over=function(n) predict(fit,n)$y;
          uniroot(function(n) n2over(n)-1.25,interval=c(10,200))$root});
        ## save supporting values in table
        names(meand20)=paste(sep='_','md20',names(meand20));
        names(over20)=paste(sep='_','ov20',names(over20));
        names(nover)=paste(sep='_','nov',names(nover));
        data.frame(what,d.crit=d_crit(20),t(meand20),t(over20),t(nover));
      }));
      dotbl(support);
    }
  });
} 
## generate title for doc_readme
## use CAP arg names to reduce conflicts with partial arg matching
title_readme=function(TEXT=NULL,...,SEP=' ') {
  dots=unlist(list(...));
  fig=paste(sep='','Figure ',figlabel());
  if (!is.null(dots)) {
    dots=paste(collapse=', ',sapply(names(dots),function(name) paste(sep='=',name,dots[name])));
    TEXT=paste(collapse='. ',c(TEXT,dots));
  }
  paste(collapse="\n",c(fig,TEXT));
}
