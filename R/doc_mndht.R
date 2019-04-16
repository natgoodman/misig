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
doc_mndht=function(sect=parent(sect,NULL)) {
  ## NG 19-04-15: use general get_data instead of special get_meand_theo
  ##   first step toward full switchover to save_data/get_data
  ## meand=get_meand_theo();
  get_data(meand,power);
  ## meand inflation - ratio of meand to true d
  meand[,cq(over,over.tval)]=with(meand,cbind(meand/d.het,meand.tval/d.het));
  ## power inflation - ratio of computed power to actual power
  power[,cq(over,over.tval)]=
    with(power,cbind(power_d2t(n,d.het)/power,power_d2t(n,d.het)/power.tval));

  ## meand.byd=split(meand.theo,meand.theo$d0);
  ## support=list(
  ##   ## support statements in text or drawn on figures
  ##   ## do it upfront 'cuz some used in figures
  ##   ## d_crit(n=20)=0.64
  ##   d.crit1=d_crit(n.fig1),
  ##   ## mean observed significant effect size and overestimate for n=20
  ##   ##   d=0.3, mean=0.81, over=2.7x
  ##   ##   d=0.5, mean=0.86, over=1.7x
  ##   ##   d=0.7, mean=0.93, over=1.3x
  ##   meand20=predict(meand.e.fit,20)$y,
  ##   over20=predict(over.e.fit,20)$y,
  ##   ## n that achieves over=1.25x
  ##   ##   d=0.3 n=122
  ##   ##   d=0.5 n=47
  ##   ##   d=0.7 n=26
  ##   n2over=function(n) predict(over.e.fit,n)$y,
  ##   nover=uniroot(function(n) n2over(n)-1.25,interval=c(10,200))$root,
  ##   end=NULL);               # placeholder for last entry
  ## dotbl(support,obj.ok=T);

  ## draw the figures
  ## figure 1
  n.fig=200; d.fig=0.3; sd.fig=0.2;
  title=title_mndht('P-values improve as observed effect size grows more extreme',
                    n=n.fig,d.het=d.fig,sd.het=sd.fig);
  sim=get_sim_hetd(n=n.fig,d=d.fig,sd=sd.fig);
  x='d.sdz'; y='d.pop';
  d.crit=d_crit(n.fig); d.htcrit=d_htcrit(n.fig,sd.fig);
  vl=c(d.fig,-d.crit,d.crit,-d.htcrit,d.htcrit);
  hl=d.fig;
  xlim=c(-1,1.5);
  dofig(plotdvsd,'big_picture',sim=sim,x=x,y=y,vline=vl,hline=hl,xlim=xlim,
        title=title,cex.main=0.75);
  ## figure 2
  title=title_mndht('Histogram of observed effect size',n=n.fig,d.het=d.fig,sd.het=sd.fig1);
  ylim=c(0,d_d2ht(n=n.fig,d.het=d.fig,sd.het=sd.fig,d=d.fig));  # set ylim to match figure 3
  dofig(plothist,'hist',sim=sim,vline=vl,title=title,cex.main=0.75,xlim=xlim,ylim=ylim);
  ## figure 3
  figblk_start();
  n.fig=200; d.fig=c(0.3,0.3,0.3,0.7,0.7); sd.fig=c(0,0.1,0.2,0.1,0.2);
  text.fig=list(label=paste_nv('sd.het',sd.fig),
                y=c(3.5,2.25,0.75,2.25,0.75),side=cq(left,left,left,right,right));
  title=title_mndht('Sampling distributions with varying sd.het and d',n=n.fig);
  dofig(plotsmpldist,'smpl_dist',n=n.fig,d.het=d.fig,sd.het=sd.fig,
        text=text.fig,vline=c(-d.htcrit,d.fig,d.htcrit),title=title,xlim=xlim,ylim=ylim);
  ##
  n.fig=c(200,20,20); d.fig=c(0.3,0.3,0.7); sd.fig=0.2;
  text.fig=list(label=paste_nv('n',n.fig),y=c(1.5,0.5,0.5),side=cq(left,left,right));
  title=title_mndht('Sampling distributions with varying n and d',sd.het=sd.fig);
  dofig(plotsmpldist,'smpl_dist',n=n.fig,d.het=d.fig,sd.het=sd.fig,
        text=text.fig,vline=c(-d.htcrit,d.fig,d.htcrit),title=title,xlim=xlim,ylim=ylim);
  figblk_end();

  ## figure 4 - meand
  figblk_start();
  ## 4a,b - actual meand - not inflation
  sd.fig=0.2;
  meand.byd=split(meand,meand$d.het);
  bsln.byd=lapply(meand.byd,function(meand) subset(meand,subset=sd.het==0,select=c(meand)));
  meand.byd=lapply(meand.byd,
                  function(meand) subset(meand,subset=sd.het==sd.fig,select=c(meand,meand.tval)))
  title=title_mndht('Mean significant effect size vs. n, d',sd.het=sd.fig);
  x=unique(meand$n);
  y=cbind(do.call(cbind,meand.byd),do.call(cbind,bsln.byd));  
  d.het=names(meand.byd);
  l=length(d.het);
  col=RColorBrewer::brewer.pal(l,'Set1');
  col=c(rep(col,each=2),col);
  lty=c(rep(cq(solid,dashed),l),rep('dotted',l));
  lwd=c(rep(2,2*l),rep(1,l));
  legend.labels=
    c(sapply(d.het,
             function(d.het) c(paste(sep='','het data and p-values. d.het=',d.het),
                               paste(sep='','het data, conventional p-values. d.het=',d.het))),
      paste(sep='','baseline. d.het=',d.het));
  dofig(plotm,'meand',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        xlab='sample size',ylab='effect size',smooth=F);
  ##
  d.fig=0.3;
  bsln=subset(meand,subset=(sd.het==0&d.het==d.fig),select=c(meand));
  meand.fig=subset(meand,subset=sd.het!=0&d.het==d.fig);
  meand.bysd=split(meand.fig,meand.fig$sd.het);
  meand.bysd=lapply(meand.bysd,function(meand) meand[,cq(meand,meand.tval)]) 

  title=title_mndht('Mean significant effect size vs. n, sd.het',d.het=d.fig);
  x=unique(meand$n);
  y=cbind(do.call(cbind,meand.bysd),bsln);  
  sd.het=names(meand.bysd);
  l=length(sd.het);
  col=RColorBrewer::brewer.pal(l,'Set1');
  col=c(rep(col,each=2),'black');
  lty=c(rep(cq(solid,dashed),l),'dotted');
  lwd=c(rep(2,2*l),1);
  legend.labels=
    c(sapply(sd.het,
             function(sd.het) c(paste(sep='','het data and p-values. sd.het=',sd.het),
                               paste(sep='','het data, conventional p-values. sd.het=',sd.het))),
      'baseline');
  dofig(plotm,'meand',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        xlab='sample size',ylab='effect size',smooth=F);
  ## 4c,d - meand inflation - ratio of meand to true d
  sd.fig=0.2;
  meand.byd=split(meand,meand$d.het);
  bsln.byd=lapply(meand.byd,function(meand) subset(meand,subset=sd.het==0,select=c(over)));
  over.byd=lapply(meand.byd,
                  function(meand) subset(meand,subset=sd.het==sd.fig,select=c(over,over.tval)))
  title=title_mndht('Significant effect size inflation vs. n, d',sd.het=sd.fig);
  x=unique(meand$n);
  y=cbind(do.call(cbind,over.byd),do.call(cbind,bsln.byd));  
  d.het=names(meand.byd);
  l=length(d.het);
  col=RColorBrewer::brewer.pal(l,'Set1');
  col=c(rep(col,each=2),col);
  lty=c(rep(cq(solid,dashed),l),rep('dotted',l));
  lwd=c(rep(2,2*l),rep(1,l));
  legend.labels=
    c(sapply(d.het,
             function(d.het) c(paste(sep='','het data and p-values. d.het=',d.het),
                               paste(sep='','het data, conventional p-values. d.het=',d.het))),
      paste(sep='','baseline. d.het=',d.het));
  dofig(plotm,'meand',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        xlab='sample size',ylab='inflation',smooth=F);

  ##
  d.fig=0.3;
  bsln=subset(meand,subset=(sd.het==0&d.het==d.fig),select=c(over));
  meand.fig=subset(meand,subset=sd.het!=0&d.het==d.fig);
  meand.bysd=split(meand.fig,meand.fig$sd.het);
  over.bysd=lapply(meand.bysd,function(meand) meand[,cq(over,over.tval)]) 

  title=title_mndht('Significant effect size inflation vs. n, sd.het',d.het=d.fig);
  x=unique(meand$n);
  y=cbind(do.call(cbind,over.bysd),bsln);  
  sd.het=names(meand.bysd);
  l=length(sd.het);
  col=RColorBrewer::brewer.pal(l,'Set1');
  col=c(rep(col,each=2),'black');
  lty=c(rep(cq(solid,dashed),l),'dotted');
  lwd=c(rep(2,2*l),1);
  legend.labels=
    c(sapply(sd.het,
             function(sd.het) c(paste(sep='','het data and p-values. sd.het=',sd.het),
                               paste(sep='','het data, conventional p-values. sd.het=',sd.het))),
      'baseline');
  dofig(plotm,'meand',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        xlab='sample size',ylab='inflation',smooth=F);
  figblk_end();

  ## figure 5 - power
  figblk_start();
  ## 5a,b - actual power - not inflation
  sd.fig=0.2;
  power.byd=split(power,power$d.het);
  bsln.byd=lapply(power.byd,function(power) subset(power,subset=sd.het==0,select=c(power)));
  power.byd=lapply(power.byd,
                  function(power) subset(power,subset=sd.het==sd.fig,select=c(power,power.tval)))
  title=title_mndht('Power vs. n, d',sd.het=sd.fig);
  x=unique(power$n);
  y=cbind(do.call(cbind,power.byd),do.call(cbind,bsln.byd));  
  d.het=names(power.byd);
  l=length(d.het);
  col=RColorBrewer::brewer.pal(l,'Set1');
  col=c(rep(col,each=2),col);
  lty=c(rep(cq(solid,dashed),l),rep('dotted',l));
  lwd=c(rep(2,2*l),rep(1,l));
  legend.labels=
    c(sapply(d.het,
             function(d.het) c(paste(sep='','het data and p-values. d.het=',d.het),
                               paste(sep='','het data, conventional p-values. d.het=',d.het))),
      paste(sep='','baseline. d.het=',d.het));
  dofig(plotm,'power',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        xlab='sample size',ylab='power',smooth=F);
  ##
  d.fig=0.3;
  bsln=subset(power,subset=(sd.het==0&d.het==d.fig),select=c(power));
  power.fig=subset(power,subset=sd.het!=0&d.het==d.fig);
  power.bysd=split(power.fig,power$sd.het);
  power.bysd=lapply(power.bysd,function(power) power[,cq(power,power.tval)]) 

  title=title_mndht('Power vs. n, sd.het',d.het=d.fig);
  x=unique(power$n);
  y=cbind(do.call(cbind,power.bysd),bsln);  
  sd.het=names(power.bysd);
  l=length(sd.het);
  col=RColorBrewer::brewer.pal(l,'Set1');
  col=c(rep(col,each=2),'black');
  lty=c(rep(cq(solid,dashed),l),'dotted');
  lwd=c(rep(2,2*l),1);
  legend.labels=
    c(sapply(sd.het,
             function(sd.het) c(paste(sep='','het data and p-values. sd.het=',sd.het),
                               paste(sep='','het data, conventional p-values. sd.het=',sd.het))),
      'baseline');
  dofig(plotm,'power',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        xlab='sample size',ylab='power',smooth=F);
  ## 5c,d - power inflation - ratio of computed power to actual power
  sd.fig=0.2;
  power.byd=split(power,power$d.het);
  bsln.byd=lapply(power.byd,function(power) subset(power,subset=sd.het==0,select=c(over)));
  over.byd=lapply(power.byd,
                  function(power) subset(power,subset=sd.het==sd.fig,select=c(over,over.tval)))
  title=title_mndht('Power inflation vs. n, d',sd.het=sd.fig);
  x=unique(power$n);
  y=cbind(do.call(cbind,over.byd),do.call(cbind,bsln.byd));  
  d.het=names(over.byd);
  l=length(d.het);
  col=RColorBrewer::brewer.pal(l,'Set1');
  col=c(rep(col,each=2),col);
  lty=c(rep(cq(solid,dashed),l),rep('dotted',l));
  lwd=c(rep(2,2*l),rep(1,l));
  legend.labels=
    c(sapply(d.het,
             function(d.het) c(paste(sep='','het data and p-values. d.het=',d.het),
                               paste(sep='','het data, conventional p-values. d.het=',d.het))),
      paste(sep='','baseline. d.het=',d.het));
  dofig(plotm,'over',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        xlab='sample size',ylab='inflation',smooth=F);
  ##
  d.fig=0.3;
  bsln=subset(power,subset=(sd.het==0&d.het==d.fig),select=c(over));
  power.fig=subset(power,subset=sd.het!=0&d.het==d.fig);
  power.bysd=split(power.fig,power.fig$sd.het);
  over.bysd=lapply(power.bysd,function(power) power[,cq(over,over.tval)]) 

  title=title_mndht('Power inflation vs. n, sd.het',d.het=d.fig);
  x=unique(power$n);
  y=cbind(do.call(cbind,over.bysd),bsln);  
  sd.het=names(over.bysd);
  l=length(sd.het);
  col=RColorBrewer::brewer.pal(l,'Set1');
  col=c(rep(col,each=2),'black');
  lty=c(rep(cq(solid,dashed),l),'dotted');
  lwd=c(rep(2,2*l),1);
  legend.labels=
    c(sapply(sd.het,
             function(sd.het) c(paste(sep='','het data and p-values. sd.het=',sd.het),
                               paste(sep='','het data, conventional p-values. sd.het=',sd.het))),
      'baseline');
  dofig(plotm,'over',x=x,y=y,title=title,cex.main=1,lwd=lwd,lty=lty,col=col,        
        legend='topright',legend.labels=legend.labels,
        xlab='sample size',ylab='inflation',smooth=F);
  figblk_end();

  invisible();
}
## generate title for doc_mndht
title_mndht=function(...,sep=' ',n=NULL,d.het=NULL,sd.het=NULL) {
  fig=paste(sep='','Figure ',figlabel());
  desc=paste(sep=sep,...);
  nv=nvq(sep=', ',ignore=T,n,d.het,sd.het);
  if (!nzchar(nv)) nv=NULL;
  paste(collapse="\n",c(fig,paste(collapse='. ',c(desc,nv))));
}

## TODO: change function name - getting clobbered by later doc_ file
##   push vline stuff into function
## plot sampling distributions (figure 3)
## xlim not presently used
## mostly hard-coded because no easy way to automate placement of text
plotsmpldist=function(n,d.het,sd.het,text=list(),
                      vline=NULL,title,dlim=c(-2,2),xlim,ylim) {
  ## TODO: set up empty plot before doing mapply
  plotpvsd(n=n[1],d.het=d.het[1],sd.het=sd.het[1],dlim=dlim,add=F,vline=vline,
           title=title,xlim=xlim,ylim=ylim);
  mapply(function(n,d.het,sd.het) plotpvsd(n=n,d.het=d.het,sd.het=sd.het,dlim=dlim,add=T),
         n,d.het,sd.het);
  do.call(Vectorize(dtext),c(list(n=n.fig,d.het=d.fig,sd.het=sd.fig),text));
  return();
}
  ## ## do n=100 first because it's way taller than the others
  ##  plotpvsd(n=100,d0=0.3,dlim=dlim,add=F,vline=vline,title=title);
  ##  plotpvsd(n=20,d0=0.3,dlim=dlim,add=T);
  ##  plotpvsd(n=20,d0=0.7,dlim=dlim,add=T);

  ## TODO: place text somehow
  ## dtext(n=20,d0=0.3,y=0.5,side='left');
  ## dtext(n=20,d0=0.7,y=0.5,side='right');
  ## dtext(n=100,d0=0.3,y=2.5,side='left');
## }
## helper function for placing text on sampling distributions (figure 3)
dtext=function(n,d.het,sd.het,label,y,side,cex=0.9,col='grey10') {
  ## invert d_d2ht to find d corresponding to y
  y2d=function(d) d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d)-y;
  if(side=='left') {
    d=uniroot(y2d,interval=c(-10,d.het))$root
    label=paste(sep='',label,' ');
    adj=c(1,0);
  } else {
    d=uniroot(y2d,interval=c(d.het,10))$root;
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
