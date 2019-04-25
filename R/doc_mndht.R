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

  sim=get_sim_hetd(n=n.fig,d=d.fig,sd=sd.fig);
  ## downsample to 1e4 so plotting will be fast
  m=param(m.hetd);
  if (m>1e4) sim=sim[sample.int(m,1e4),];
  ## draw the figures
  ## figure 1
  n.fig=200; d.fig=0.3; sd.fig=0.2;
  title=title_mndht('P-values improve as observed effect size grows more extreme',
                    n=n.fig,d.het=d.fig,sd.het=sd.fig);
  x='d.sdz'; y='d.pop';
  d.crit=d_crit(n.fig); d.htcrit=d_htcrit(n.fig,sd.fig);
  vl=c(d.fig,-d.crit,d.crit,-d.htcrit,d.htcrit);
  hl=d.fig;
  xlim=c(-1,1.5);
  dofig(plotdvsd,'big_picture',sim=sim,x=x,y=y,vline=vl,hline=hl,xlim=xlim,
        title=title,cex.main=0.75);
  ## figure 2
  title=title_mndht('Histogram of observed effect size',n=n.fig,d.het=d.fig,sd.het=sd.fig1);
  ylim=c(0,d_d2t(n=n.fig,d=d.fig,d0=d.fig));  # set ylim to match figure 3
  dofig(plothist,'hist',sim=sim,vline=vl,title=title,cex.main=0.75,xlim=xlim,ylim=ylim);
  ## figure 3
  figblk_start();
  n.fig=200; d.fig=c(0.3,0.3,0.3,0.7,0.7); sd.fig=c(0,0.1,0.2,0.1,0.2);
  text.fig=list(label=paste_nv('sd.het',sd.fig),
                y=c(3.5,2.25,0.75,2.25,0.75),side=cq(left,left,left,right,right));
  title=title_mndht('Sampling distributions with varying sd.het and d.het',n=n.fig);
  dofig(plotsmpldist,'smpldist',n=n.fig,d.het=d.fig,sd.het=sd.fig,
        text=text.fig,vline=c(-d.htcrit,d.fig,d.htcrit),title=title,xlim=xlim,ylim=ylim);
  ##
  n.fig=c(200,20,20); d.fig=c(0.3,0.3,0.7); sd.fig=0.2;
  text.fig=list(label=paste_nv('n',n.fig),y=c(1.5,0.5,0.5),side=cq(left,left,right));
  title=title_mndht('Sampling distributions with varying n and d.het',sd.het=sd.fig);
  dofig(plotsmpldist,'smpldist',n=n.fig,d.het=d.fig,sd.het=sd.fig,
        text=text.fig,vline=c(-d.htcrit,d.fig,d.htcrit),title=title,xlim=xlim,ylim=ylim);
  ## figure 4 - meand
  ## TODO: too many figures for post. more suitable for supp 
  figblk_start();
  ## 4a,b - meand & inflation by d
  d.fig=c(0.3,0.7); sd.fig=0.2;
  dofig(plot_byd,'meand_byd',what='meand',y='meand',d.fig=d.fig,sd.fig=sd.fig);
  dofig(plot_byd,'meandover_byd',what='meand',y='over',d.fig=d.fig,sd.fig=sd.fig);
  ## 4c,d - meand & inflation by sd
  d.fig=0.3; sd.fig=c(0.1,0.2);
  dofig(plot_bysd,'meand_bysd',what='meand',y='meand',d.fig=d.fig,sd.fig=sd.fig);
  dofig(plot_bysd,'meandover_bysd',what='meand',y='over',d.fig=d.fig,sd.fig=sd.fig);
  ## figure 5 - meand
  ## TODO: too many figures for post. more suitable for supp 
  figblk_start();
  ## 5a,b - power & inflation by d
  d.fig=c(0.3,0.7); sd.fig=0.2;
  dofig(plot_byd,'power_byd',what='power',y='power',d.fig=d.fig,sd.fig=sd.fig,legend='right');
  dofig(plot_byd,'powerover_byd',what='power',y='over',d.fig=d.fig,sd.fig=sd.fig,legend='right');
  ## 5c,d - power & inflation by sd
  d.fig=0.3; sd.fig=c(0.1,0.2);
  dofig(plot_bysd,'power_bysd',what='power',y='power',d.fig=d.fig,sd.fig=sd.fig,legend='right');
  dofig(plot_bysd,'powerover_bysd',what='power',y='over',d.fig=d.fig,sd.fig=sd.fig,legend='right');
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
plotsmpldist=function(n,d.het,sd.het,text=list(),
                      vline=NULL,title,dlim=c(-2,2),xlim,ylim) {
  ## TODO: set up empty plot before doing mapply
  plotpvsd(n=n[1],d.het=d.het[1],sd.het=sd.het[1],dlim=dlim,add=F,vline=vline,
           title=title,xlim=xlim,ylim=ylim);
  mapply(function(n,d.het,sd.het) plotpvsd(n=n,d.het=d.het,sd.het=sd.het,dlim=dlim,add=T),
         n,d.het,sd.het);
  do.call(Vectorize(dtext),c(list(n=n.fig,d.het=d.het,sd.het=sd.het),text));
  return();
}
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
## plot meand, power (raw or inflation) split by d.het
plot_byd=
  function(what=cq(meand,power),y=what,d.fig=NULL,sd.fig=0.2,
           title=NULL,xlab='sample size',ylab=NULL,cex.main=0.75,legend='topright',...) {
    what=match.arg(what);
    y=match.arg(y,cq(meand,power,over));
    if (is.null(title)) title=title_byd(what,y,sd.fig);
    if (is.null(ylab)) ylab=ylab_byd(what,y);
    data=get(what,envir=parent.frame(n=1));
    if (!is.null(d.fig)) data=subset(data,subset=(d.het%in%d.fig));
    data.byd=split(data,data$d.het);
    data.byd=lapply(data.byd,function(dat) {
      data=subset(dat,subset=(sd.het==sd.fig),select=c(y,paste0(y,'.tval')));
      bsln=subset(dat,subset=(sd.het==0),select=c(y));
      colnames(bsln)='bsln';
      cbind(data,bsln);
    });
    x=unique(data$n);
    y=do.call(cbind,data.byd);  
    ## title=title_mndht('Mean significant effect size vs. n, d',sd.het=sd.fig);
    d.het=names(data.byd);
    l=length(d.het);
    col=RColorBrewer::brewer.pal(max(l,3),'Set1')[1:l];
    lty=cq(solid,dashed,dotted);
    lwd=c(2,2,1);
    legend.args=list(labels=c('colors (represent d.het)',
                              paste_nv('d.het'),
                              'line types',
                              'het p-values','conventional p-values','baseline'),
                     col=c('black',col,'black',rep('grey20',l)),
                     lty=c(NA,rep('solid',l),NA,lty),
                     lwd=c(NA,rep(lwd,l),NA,lwd));
    col=rep(col,each=3);
    lty=rep(lty,l);
    lwd=rep(lwd,l);
    plotm(x=x,y=y,title=title,cex.main=cex.main,lwd=lwd,lty=lty,col=col,        
          legend=legend,legend.args=legend.args,
          xlab=xlab,ylab=ylab,smooth=F,...); 
  }
title_byd=function(what,y,sd.het) {
  desc=if (what=='meand') 'Mean significant effect size' else 'Power';
  if (y=='over') desc=paste(sep=' ',desc,'inflation');
  title_mndht(desc,'vs. n, d.het',sd.het=sd.het);
}
ylab_byd=function(what,y) {
  ylab=if (what=='meand') 'effect size' else 'power';
  if (y=='over') ylab=paste(sep=' ',ylab,'inflation');
  ylab;
}
## plot meand, power (raw or inflation) split by sd.het
plot_bysd=
  function(what=cq(meand,power),y=what,d.fig=0.3,sd.fig=NULL,
           title=NULL,xlab='sample size',ylab=NULL,cex.main=0.75,legend='topright',...) {
    what=match.arg(what);
    y=match.arg(y,cq(meand,power,over));
    if (is.null(title)) title=title_bysd(what,y,d.fig);
    if (is.null(ylab)) ylab=ylab_bysd(what,y);
    data=get(what,envir=parent.frame(n=1));
    data=subset(data,subset=(d.het==d.fig));
    data.bsln=split(data,data$sd.het==0);
    bsln=data.bsln[['TRUE']][,y];
    data=data.bsln[['FALSE']];

    if (!is.null(sd.fig)) data=subset(data,subset=(sd.het%in%sd.fig));
    data.bysd=split(data,data$sd.het);
    data.bysd=lapply(data.bysd,function(data) subset(data,select=c(y,paste0(y,'.tval'))));
    x=unique(data$n);
    y=cbind(do.call(cbind,data.bysd),bsln);  
    sd.het=names(data.bysd);
    l=length(sd.het);
    col=RColorBrewer::brewer.pal(max(l,3),'Set1')[1:l];
    lty=cq(solid,dashed,dotted);
    legend.args=list(labels=c('colors (represent sd.het)',
                              paste_nv('sd.het'),
                              'line types',
                              'het p-values','conventional p-values','baseline'),
                     col=c('black',col,'black',rep('grey20',3)),
                     lty=c(NA,rep('solid',l),NA,lty),
                     lwd=c(NA,rep(2,l),NA,2,2,1));
    col=c(rep(col,each=2),'black');
    lty=c(rep(cq(solid,dashed),l),'dotted');
    lwd=c(rep(c(2,2),l),1);
    plotm(x=x,y=y,title=title,cex.main=cex.main,lwd=lwd,lty=lty,col=col,        
          legend=legend,legend.args=legend.args,
          xlab=xlab,ylab=ylab,smooth=F,...); 
  }
title_bysd=function(what,y,d.het) {
  desc=if (what=='meand') 'Mean significant effect size' else 'Power';
  if (y=='over') desc=paste(sep=' ',desc,'inflation');
  title_mndht(desc,'vs. n, sd.het',d.het=d.het);
}
ylab_bysd=function(what,y) {
  ylab=if (what=='meand') 'effect size' else 'power';
  if (y=='over') ylab=paste(sep=' ',ylab,'inflation');
  ylab;
}

