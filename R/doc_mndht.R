#################################################################################
##
## Author:  Nat Goodman
## Created: 19-04-09
##          from doc_ovrht.R created 19-02-19
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
doc_ovrht=function(n.fig1=20,d.fig1=0.3,sd.fig1=0.2,sect=parent(sect,NULL)) {
  meand.empi=get_meand_empi();
  meand.theo=get_meand_theo();

    ##### this block not yet ported

  meand.e.byd=split(meand.empi,meand.empi$d0);
  meand.t.byd=split(meand.theo,meand.theo$d0);
  ## make splines that interpolate meand, over
  ##   empirically determined smooth.spline as interp function, and spar=0.3
  spar=0.3;
  meand.e.fit=lapply(meand.e.byd,function(meand) with(meand,smooth.spline(n,meand,spar=spar)));
  ## over.e.fit=lapply(meand.e.byd,function(meand) with(meand,smooth.spline(n,over,spar=spar)));
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
  meand20=predict(meand.e.fit,20)$y;
  over20=predict(over.e.fit,20)$y;
  ## n that achieves over=1.25x
  ##   d=0.3 n=122
  ##   d=0.5 n=47
  ##   d=0.7 n=26
  n2over=function(n) predict(over.e.fit,n)$y;
  nover=uniroot(function(n) n2over(n)-1.25,interval=c(10,200))$root;
  ## save supporting values in table
  support=data.frame(d.crit1,meand20,over20,nover);
  dotbl(support);
  ##### 

  ## draw the figures
  ## figure 1
  title=title_ovrht('Histogram and sampling distribution under NULL',d.het=NULL,sd.het=sd.fig1,n=n.fig1);
  sim=get_sim_hetd(n=n.fig1,d=0,sd=sd.fig1)
  dofig(plothist_d2ht,'histd2ht_NULL',sim=sim,title=title,n=n.fig1,d.het=0,sd.het=c(0,sd.fig1));
  ## figure 2
  title=title_ovrht('Histogram and sampling distribution non-NULL',
                    d.het=d.fig1,sd.het=sd.fig1,n=n.fig1);
  dofig(plothist_d2ht,'histd2ht_nonNULL',title,n=n.fig1,d.het=d.fig1,sd.het=c(0,sd.fig1));

  ## figure 3
  title=title_ovrht('Sampling distributions show impact of increasing n and d');
  dofig(plotsmpldist,'smpl_dist',vline=c(-d.crit1,d.fig1,d.crit1),title=title);

 ## figure 4
  title=title_ovrfx('Average observed effect size improves as n increases');
  x=seq(min(n.hetd),max(n.hetd),by=1);
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
## TODO: probably don't need d
## generate title for doc_ovrht
title_ovrht=function(desc=NULL,d.het=NULL,sd.het=NULL,n=NULL) {
  fig=paste(sep='','Figure ',figlabel());
  paste(collapse="\n",
        c(fig,
          paste(sep='. ',desc,
                paste(collapse=', ',
                      c(if(!is.null(d.het))paste_nv(d.het),
                        if(!is.null(sd.het))paste_nv(sd.het),
                        if(!is.null(n))paste_nv(n))))));
}


## TODO: NOT USED. wrong sampling distribution!!
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
dtext=
  function(n,d0,y,label=paste(sep=', ',paste_nv('d',d0),paste_nv('n',n)),side=cq(left,right),
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

## figure 1, 2. plot histogram and d_d2t distribution
plothist_d2ht=
  function(sim,title,n,d.het,sd.het,xlim=d.het+c(-1.5,1.5),lwd.sig=2,
           plot.dcrit=T,plot.dcplus=plot.dcrit,plot.dcminus=plot.dcrit,
           dc.col='black',dc.lwd=1,dc.lty='dotted',dc.cex=0.75,dc.text=F) {
    ylim=c(0,d_d2t(n=n,d=0));             # set ylim to fit d2t
    plothist(sim=sim,col='grey90',border='grey80',title=title,cex.main=1,xlim=xlim,ylim=ylim);
    sapply(sd.het,function(sd.het) {
      plotpvsd(n=n,d.het=d.het,sd.het=sd.het,dlim=xlim,lwd.sig=lwd.sig,add=T);
      if (plot.dcplus|plot.dcminus) {
        d.crit=d_htcrit(n,sd.het);
        dc=NULL;
        if (plot.dcplus) dc=c(dc,d.crit);
        if (plot.dcminus) dc=c(dc,-d.crit);
        y=d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=dc);
        text=if(dc.text) round(dc,digits=2) else NULL;
        vline(x=dc,y=y,col=dc.col,lty=dc.lty,lwd=dc.lwd,cex=dc.cex,text=text);
      }});
  invisible();
}

## plot meand vs n (figure 4)
plotmeandn=
  function(n,d.het,sd.het,title,col='black',lty='solid',lwd=1,legend.labels=NULL,...) {
    cases=expand.grid(d.het=d.het,sd.het=sd.het);
    x=seq(min(n),max(n),by=1);
    meand=do.call(cbind,lapply(d.het,function(d.het) {
      print(d.het);
      d.crit=d_htcrit(n=n,sd.het=sd.het);
      meand.het=meand_d2ht(n,d.het,sd.het,d.crit);
      ## do it next with t-pvals
      d.crit=d_crit(n=n);
      meand.t=meand_d2ht(n,d.het,sd.het,d.crit);
      ## finally, t data and pvals
      meand.tt=meand_d2ht(n,d.het,sd.het=0,d.crit);
      cbind(meand.het,meand.t,meand.tt);
    }));
    meand=splinem(n,meand,xout=x);

    overd=do.call(cbind,lapply(d.het[d.het>0],function(d.het) {
      print(d.het);
      d.crit=d_htcrit(n=n,sd.het=sd.het);
      meand.het=meand_d2ht(n,d.het,sd.het,d.crit);
      ## do it next with t-pvals
      d.crit=d_crit(n=n);
      meand.t=meand_d2ht(n,d.het,sd.het,d.crit);
      ## finally, t data and pvals
      meand.tt=meand_d2ht(n,d.het,sd.het=0,d.crit);
      cbind(meand.het,meand.t,meand.tt)/d.het;
    }));
    overd=splinem(n,overd,xout=x);

    ## draw the main plot
  plotm(x=x,y=y,col=col,lty=lty,lwd=lwd,title=title,legend.labels=legend.labels,smooth=F,
        ...);
  ## horizontal grid-like lines for d.pop
  abline(h=d.mean,col=col,lty='dotted',lwd=1);
  ## horizontal lines for averages with n=20
  ## if (!is.null(meand20)) {
  ##   hline(x=20,y=meand20,col=col,lty='dotted',lwd=0.5,text=round(meand20,digits=2));
  ## }
  ## vertical line for n=20
  vline(x=20,y0=0,y=1,col='grey',lty='dotted',lwd=1,text=20);
  ## horizontal & vertical lines for overestimates
  if (!is.null(nover)) {
    hline(x=nover,y=1.25*d.mean,col=col,lty='dotted',lwd=2,text=paste(sep='','1.25x',d.mean));
    vline(x=nover,y=1.25*d.mean,col=col,lty='dotted',lwd=2,text=round(nover));
  }
}
