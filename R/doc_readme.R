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
    sectnum=which(sect==sect.all)[1];
    sect_start(sectnum); 
##### plotdvsd
    if (sect=='plotdvsd') {
      ## n.rand should have two elements. make sure, else Rmd will break 
      n=param(n.rand); n=n[c(1,length(n))];
      sapply(n,function(n) {
        figblk_start();
        d.crit=d_crit(n);
        sim=get_sim_rand(n=n);
        x='d.sdz'; y='d.pop';
        title=title_readme(c('Scatter plot.',y,'vs',x),n=n);
        dofig(plotdvsd,'scat_pop_sdz',sim=sim,x=x,y=y,vline=c(-d.crit,d.crit),xlim=c(-2,2),
              title=title,cex.main=1);
        x='d.pop'; y='d.sdz';
        title=title_readme(c('Scatter plot.',y,'vs',x),n=n);
        dofig(plotdvsd,'scat_sdz_pop',sim=sim,x=x,y=y,hline=c(-d.crit,d.crit),ylim=c(-2,2),
              title=title,cex.main=1);
      });
      figblk_start();
      sapply(n,function(n) {
        d.crit=d_crit(n);
        sim=get_sim_rand(n=n);
        sim.zoom=subset(sim,subset=near(d.sdz,d.crit,0.10));
        x='d.sdz'; y='d.pop';
        title=title_readme(c('Scatter plot.',y,'vs',x,'showing sharp boundary at'),
                           round(d.crit,digits=2),n=n);
        dofig(plotdvsd,'zoom',sim=sim.zoom,x=x,y=y,vline=d.crit,vlab=T,legend=F,
              title=title,cex.main=0.9);
      });
    }
##### plothist
    if (sect=='plothist') {
      ## n.fixd, d.fixd too big to plot all. grab 1st and last elements.
      n=param(n.fixd); n=n[c(1,length(n))];
      d=param(d.fixd); d=d[c(1,length(d))];
      sapply(n,function(n) {
        figblk_start();
        d.crit=d_crit(n);
        d.fig=d[1];
        sim=get_sim_fixd(n=n,d=d.fig);
        title=title_readme('Histogram of observed effect size',d=d.fig,n=n);
        dofig(plothist,'hist',sim=sim,vline=c(-d.crit,d.fig,d.crit),title=title,cex.main=1,
              breaks=20);
        d.fig=d[2];
        sim=get_sim_fixd(n=n,d=d.fig);
        title=title_readme('Histogram of observed effect size',d=d.fig,n=n);
        dofig(plothist,'hist',sim=sim,vline=c(-d.crit,d.fig,d.crit),title=title,cex.main=1,
              breaks=20);
      });
    }
##### plotpvsd 
    if (sect=='plotpvsd') {
      ## n.fixd, d.fixd too big to plot all. grab 1st and last elements.
      n=param(n.fixd); n=n[c(1,length(n))];
      d=param(d.fixd); d=d[c(1,length(d))];
      cases=data.frame(n=n,d=d);
      apply(cases,1,function(case) {
        figblk_start();
        n=case[1]; d=case[2];
        d.crit=d_crit(n);
        title=title_readme('NULL distribution (density)',d=d,n=n);
        dofig(plotpvsd,'distn',n=n,d0=0,y='density',vline=c(-d.crit,d.crit),xlim=c(-1,2),
              title=title,cex.main=1);
        title=title_readme('NULL distribution (cumulative)',d=d,n=n);
        dofig(plotpvsd,'distn',n=n,d0=0,y='cumulative',vline=c(-d.crit,d.crit),xlim=c(-1,2),
              title=title,cex.main=1);
        title=title_readme('Sampling distribution (density)',d=d,n=n);
        dofig(plotpvsd,'distn',n=n,d0=d,y='density',vline=c(-d.crit,d.crit),xlim=c(-1,2),
              title=title,cex.main=1);
        title=title_readme('Sampling distribution (cumulative)',d=d,n=n);
        dofig(plotpvsd,'distn',n=n,d0=d,y='cumulative',vline=c(-d.crit,d.crit),xlim=c(-1,2),
              title=title,cex.main=1);
      });
    }
##### plotm 
    if (sect=='plotm') {
      param(n.meand,d.meand);
      meand.empi=get_meand_empi();
      meand.theo=get_meand_theo();
      meand.e.byd=split(meand.empi,meand.empi$d0);
      meand.t.byd=split(meand.theo,meand.theo$d0);
      x=n.meand;
      y=do.call(cbind,c(lapply(meand.e.byd,function(meand) meand$meand),
                        lapply(meand.t.byd,function(meand) meand$meand)));
      col=setNames(RColorBrewer::brewer.pal(length(d.meand),'Set1'),d.meand);
      legend.labels=
        c(paste(sep='=','meand.empi. d.pop',d.meand),paste(sep='=','meand.theo. d.pop',d.meand));
      col=rep(col,2);
      lty=c(rep('dotted',len=length(d.meand)),rep('solid',length(d.meand)));
      lwd=2;     
      title=title_readme('Line plot of mean significant observed effect size (raw)');
      dofig(plotm,'meand',x=x,y=y,col=col,lty=lty,lwd=lwd,smooth=F,title=title,
            legend.labels=legend.labels,legend='topright',
            xlab='sample size',ylab='effect size');
      title=title_readme('Line plot of mean significant observed effect size (smooth=aspline)');
      dofig(plotm,'meand',x=x,y=y,col=col,lty=lty,lwd=lwd,smooth=T,title=title,
            legend.labels=legend.labels,legend='topright',
            xlab='sample size',ylab='effect size');
      title=
        title_readme('Line plot of mean significant observed effect size (smooth=smooth.spline)');
      dofig(plotm,'meand',x=x,y=y,col=col,lty=lty,lwd=lwd,smooth='spline',title=title,
            legend.labels=legend.labels,legend='topright',
            xlab='sample size',ylab='effect size');
      title=title_readme('Line plot of mean significant observed effect size (smooth=loess)');
      suppressWarnings(
        dofig(plotm,'meand',x=x,y=y,col=col,lty=lty,lwd=lwd,smooth='loess',title=title,
            legend.labels=legend.labels,legend='topright',
            xlab='sample size',ylab='effect size'));
    }
    if (sect=='table') {
      param(n.meand,d.meand);
      meand.empi=get_meand_empi();
      meand.theo=get_meand_theo();
      meand.e.byd=split(meand.empi,meand.empi$d0);
      meand.t.byd=split(meand.theo,meand.theo$d0);
      ## make splines that interpolate meand, over
      ##   empirically determined smooth.spline as interp function, and spar=0.3
      spar=0.3;
      meand.e.fit=lapply(meand.e.byd,function(meand) with(meand,smooth.spline(n,meand,spar=spar)));
      meand.t.fit=lapply(meand.t.byd,function(meand) with(meand,smooth.spline(n,meand,spar=spar)));
      over.t.fit=lapply(meand.t.byd,function(meand) with(meand,smooth.spline(n,over,spar=spar)));
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
      names(meand20)=paste(sep='_','md20',names(meand20));
      names(over20)=paste(sep='_','ov20',names(over20));
      names(nover)=paste(sep='_','nov',names(nover));
      support=data.frame(d.crit=d_crit(20),t(meand20),t(over20),t(nover));
      dotbl(support);
    }
  });
}
 
## generate title for doc_readme
title_readme=function(desc=NULL,d=NULL,n=NULL,sep=' ') {
  fig=paste(sep='','Figure ',figlabel());
  paste(collapse="\n",
        c(fig,
          paste(sep='. ',paste(collapse=sep,desc),
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
  param(n.meand,d.meand);
  ## draw the main plot
  plotm(x=x,y=y,col=col,lty=lty,lwd=lwd,title=title,legend.labels=legend.labels,smooth=F,
        ...);
  ## horizontal grid-like lines for d.pop
  abline(h=d.meand,col=col,lty='dotted',lwd=1);
  ## horizontal lines for averages with n=20
  ## if (!is.null(meand20)) {
  ##   hline(x=20,y=meand20,col=col,lty='dotted',lwd=0.5,text=round(meand20,digits=2));
  ## }
  ## vertical line for n=20
  vline(x=20,y0=0,y=1,col='grey',lty='dotted',lwd=1,text=20);
  ## horizontal & vertical lines for overestimates
  if (!is.null(nover)) {
    hline(x=nover,y=1.25*d.meand,col=col,lty='dotted',lwd=2,text=paste(sep='','1.25x',d.meand));
    vline(x=nover,y=1.25*d.meand,col=col,lty='dotted',lwd=2,text=round(nover));
  }
}
