#################################################################################
##
## Author:  Nat Goodman
## Created: 19-05-09
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
## --- Generate Figures and Tables for misig README supplement ---
doc_readmesupp=function(sect=NULL) {
  sect.all=cq(plotdvsd,plothist,plotpvsd,plotm,plotcustom,table);
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all);
  sapply(sect,function(sect) {
    ## compute section number. from stackoverflow.com/questions/5577727
    sectnum=which(sect==sect.all)[1];
    sect_start(sectnum); 
##### plotdvsd
    ## each sim type expected to have >= 2 cases
    if (sect=='plotdvsd') {
      ## fixd doesn't make much sense since there's no d.pop spread, but why not?
      figblk_start();
      param(n.fixd,d.fixd);
      ## use small n else all points significant
      cases=data.frame(n=sort(n.fixd)[1:2],d=range(d.fixd),zoom=c(FALSE,FALSE,TRUE,TRUE));
      apply(cases,1,function(case) withcase(case, {
        sim=get_sim_fixd(n=n,d=d);
        figdvsd(sim=sim,stext='fixd',d.crit=d_crit(n),zoom=zoom,n=n,d.pop=d);
      }));
      ## rand
      figblk_start();
      n=range(param(n.rand));
      cases=data.frame(n=range(param(n.rand)),zoom=c(FALSE,FALSE,TRUE,TRUE));
      apply(cases,1,function(case) withcase(case, {
        sim=get_sim_rand(n=n);
        figdvsd(sim=sim,stext='rand',d.crit=d_crit(n),zoom=zoom,n=n);
      }));
      ## hetd
      figblk_start();
      param(n.hetd,d.hetd,sd.hetd);
      cases=data.frame(n=range(n.hetd),d=range(d.hetd),sd=range(sd.hetd[sd.hetd!=0]),
                       zoom=c(FALSE,FALSE,TRUE,TRUE));
      apply(cases,1,function(case) withcase(case, {
        sim=get_sim_hetd(n=n,d=d,sd=sd);
        figdvsd(sim=sim,stext='hetd',d.crit=d_htcrit(n,sd),zoom=zoom,n=n,d.het=d,sd.het=sd);
      }));
    }
##### plothist
    if (sect=='plothist') {
      ## fixd
      figblk_start();
      param(n.fixd,d.fixd);
      ## use small n else all points significant
      cases=data.frame(n=sort(n.fixd)[1:2],d=range(d.fixd),zoom=c(FALSE,FALSE,TRUE,TRUE));
      apply(cases,1,function(case) withcase(case, {
        sim=get_sim_fixd(n=n,d=d);
        fighist(sim=sim,stext='fixd',d.crit=d_crit(n),zoom=zoom,n=n,d.pop=d);
      }));
      ## rand
      figblk_start();
      n=range(param(n.rand));
      cases=data.frame(n=range(param(n.rand)),zoom=c(FALSE,FALSE,TRUE,TRUE));
      apply(cases,1,function(case) withcase(case, {
        sim=get_sim_rand(n=n);
        fighist(sim=sim,stext='rand',d.crit=d_crit(n),zoom=zoom,n=n);
      }));
      ## hetd
      figblk_start();
      param(n.hetd,d.hetd,sd.hetd);
      cases=data.frame(n=range(n.hetd),d=range(d.hetd),sd=range(sd.hetd[sd.hetd!=0]),
                       zoom=c(FALSE,FALSE,TRUE,TRUE));
      apply(cases,1,function(case) withcase(case, {
        sim=get_sim_hetd(n=n,d=d,sd=sd);
        fighist(sim=sim,stext='hetd',d.crit=d_htcrit(n,sd),zoom=zoom,zbreaks=25,
                n=n,d.het=d,sd.het=sd);
      }));
    }
##### plotpvsd 
    if (sect=='plotpvsd') {
      ## fixd
      figblk_start();
      param(n.fixd,d.fixd);
      ## use small n else all points significant
      cases=data.frame(n=sort(n.fixd)[1:2],d=range(d.fixd));
      apply(cases,1,function(case) withcase(case, {
        figpvsd(mtext='fixd',d.crit=d_crit(n),n=n,d0=d);
      }));
      ## rand - can't do it - each instance has different d.pop
      ## hetd
      figblk_start();
      param(n.hetd,d.hetd,sd.hetd);
      cases=data.frame(n=range(n.hetd),d=range(d.hetd),sd=range(sd.hetd[sd.hetd!=0]));
      apply(cases,1,function(case) withcase(case, {
        figpvsd(mtext='hetd',d.crit=d_htcrit(n,sd),n=n,d.het=d,sd.het=sd);
      }));
    }
##### plotm
    if (sect=='plotm') {
      ## fixd
      figblk_start();
      sapply(cq(none,aspline,spline),function(smooth) {
        ## these show concordance of simulated vs theoretical results
        ## meand
        simu=get_data(meand.fixd); theo=get_data(meand.d2t);
        figm_simutheo(simu,theo,stat='meand',mtext='fixd',smooth=smooth);
        ## pval
        simu=get_data(pval.fixd); theo=get_data(pval.d2t);
        figm_simutheo(simu,theo,stat='pval',mtext='fixd',ylim=c(0,0.1),smooth=smooth);
        ## power
        simu=get_data(power.fixd); theo=get_data(power.d2t);
        figm_simutheo(simu,theo,stat='power',mtext='fixd',smooth=smooth);
        ## this shows ci coverage for all vs sig
        ci=get_data(ci.fixd); ci=subset(ci,subset=d.pop>0);
        figm_ci(ci,mtext='fixd',legend='right',smooth=smooth);
      });
      ## hetd
      figblk_start();
      sapply(cq(none,aspline,spline),function(smooth) {
        ## meand
        simu=get_data(meand.hetd); theo=get_data(meand.d2ht);
        simu=subset(simu,subset=d.het==0); theo=subset(theo,subset=d.het==0); 
        figm_simutheo(simu,theo,stat='meand',mtext='hetd',smooth=smooth,by='sd.het',d.het=0);
        ## pval. pval.tval - more interesting
        simu=get_data(pval.hetd); theo=get_data(pval.d2ht);
        figm_simutheo(simu,theo,stat='pval',mtext='hetd',by='sd.het',ylim=c(0,0.1),smooth=smooth);
        ## power
        simu=get_data(power.hetd); theo=get_data(power.d2ht);
        simu=subset(simu,subset=d.het==0.5); theo=subset(theo,subset=d.het==0.5); 
        figm_simutheo(simu,theo,stat='power',mtext='hetd',smooth=smooth,by='sd.het',
                      legend='right',d.het=0.5);
        ## this shows ci coverage for all vs sig
        ci=get_data(ci.hetd); ci=subset(ci,subset=d.het==0.5);
        figm_ci(ci,mtext='hetd',legend='bottomright',smooth=smooth,by='sd.het',d.het=0.5);
      });
    }
    if (sect=='plotcustom') {
      ## custom plots from other docs
      ## from ovrfx
      figblk_start();
      dofig(plotsmpldist,'ovrfx_smpldist');
      dofig(plotmeand,'ovrfx_meand');
      ## from ovrht
      figblk_start();
      dofig(plothist_d2t,'ovrht_histd2t');
      dofig(plotpval_over,'ovrht_pvalover');
      dofig(plotci_over,'ovrht_ciover');
      dofig(plotd2t_d2ht,'ovrht_sampldist');
    }
    if (sect=='table') {
      ## flat table
      param(n.fixd,d.fixd);
      n.tbl=min(n.fixd);              # typically 20, hence the variable names below
      flat=do.call(rbind,lapply(cq(meand.fixd,meand.d2t),function(what) {
        meand=get_data(list=what);
        meand=subset(meand,subset=d0>0);
        meand.byd=split(meand,meand$d0);
        ## make splines that interpolate meand, over
        ##   empirically determined smooth.spline as interp function, and spar=0.3
        spar=0.3;
        meand.fit=lapply(meand.byd,function(meand) with(meand,smooth.spline(n,meand,spar=spar)));
        over.fit=lapply(meand.byd,function(meand) with(meand,smooth.spline(n,over,spar=spar)));
        meand20=sapply(meand.fit,function(fit) predict(fit,n.tbl)$y);
        over20=sapply(over.fit,function(fit) predict(fit,n.tbl)$y);
        nover=sapply(over.fit,function(fit) {
          n2over=function(n) predict(fit,n)$y;
          uniroot(function(n) n2over(n)-1.25,interval=c(10,200))$root});
        ## save supporting values in table
        names(meand20)=paste(sep='_','md20',names(meand20));
        names(over20)=paste(sep='_','ov20',names(over20));
        names(nover)=paste(sep='_','nov',names(nover));
        data.frame(what,d.crit=d_crit(20),t(meand20),t(over20),t(nover));
      }));
      dotbl(flat);
      ## object table
      param(n.hetd,d.hetd,sd.hetd);
      meand=get_data(meand.d2ht);
      pval=get_data(pval.d2ht);
      ci=get_data(ci.d2ht);
      cases=data.frame(n.tbl=range(n.hetd),d.tbl=range(d.hetd),sd.tbl=range(sd.hetd));
      obj1=with(cases[1,],{
        list(
          meand1=subset(meand,subset=(n==n.tbl&sd.het==sd.tbl)),
          pval1=subset(pval,subset=(n==n.tbl&sd.het==sd.tbl)),
          ci1=subset(ci,subset=(n==n.tbl&sd.het==sd.tbl))
          )});
      obj2=with(cases[2,],{
        list(
          meand2=subset(meand,subset=(n==n.tbl&sd.het==sd.tbl)),
          pval2=subset(pval,subset=(n==n.tbl&sd.het==sd.tbl)),
          ci2=subset(ci,subset=(n==n.tbl&sd.het==sd.tbl))
        )});
      obj=c(obj1,obj2);
      obj=lapply(obj,function(x) round(x,digits=2));
      dotbl(obj,obj.ok=T);
    }
  });
}

## wrapper for plotdvsd
figdvsd=function(sim,stext,d.crit,zoom,...) {
    ztext=NULL;                       # assume not zoom
    if (zoom) {                       # change ztext if zoom
      sim=subset(sim,subset=near(d.sdz,d.crit,0.10*d.crit));
      if (nrow(sim)==0) return();   # nothing to plot
      ztext=paste(sep=' ','showing sharp boundary at',round(d.crit,digits=2))
    } 
    x='d.sdz'; y='d.pop';
    title=figtitle(c('Scatter plot',y,'vs',x,ztext),sim=stext,...);
    figname=paste(collapse='_',c('pop_sdz',stext,if(zoom) 'zoom'));
    dofig(plotdvsd,figname,sim=sim,x=x,y=y,vline=c(-d.crit,d.crit),title=title);
    x='d.pop'; y='d.sdz';
    title=figtitle(c('Scatter plot',y,'vs',x,ztext),sim=stext,...);
    figname=paste(collapse='_',c('sdz_pop',stext,if(zoom) 'zoom'));
    dofig(plotdvsd,figname,sim=sim,x=x,y=y,hline=c(-d.crit,d.crit),title=title);
    return();
}
## wrapper for plothist
fighist=function(sim,stext,d.crit,zoom,zbreaks=10,zdigits=3,...) {
  ## set variables for not zoom
  ztext=NULL; breaks=formals(plothist)$breaks; legend=formals(plothist)$legend;
  vhdigits=formals(plothist)$vhdigits;
  if (zoom) {               # change variables if zoom
    sim=subset(sim,subset=near(d.sdz,d.crit,0.10*d.crit));
    if (nrow(sim)==0) return();   # nothing to plot
    ztext=paste(sep=' ','showing sharp boundary at about',round(d.crit,digits=2))
    breaks=zbreaks;
    legend=FALSE;
    vhdigits=zdigits;
  } 
  title=figtitle(c('Histogram of observed effect size',ztext),sim=stext,...);
  figname=paste(collapse='_',c('hist',stext,if(zoom) 'zoom'));
  dofig(plothist,figname,sim=sim,vline=c(-d.crit,d.crit),breaks=breaks,legend=legend,
        title=title,vhdigits=vhdigits);
  return();
}
## wrapper for plotpvsd
figpvsd=function(n,d0=NULL,d.het=NULL,sd.het=NULL,mtext,d.crit,...) {
  title=figtitle('Sampling distribution (density)',model=mtext,
                 n=n,d0=d0,d.het=d.het,sd.het=sd.het,...);
  figname=paste(collapse='_',c('dens',mtext));
  dofig(plotpvsd,figname,n=n,d0=d0,d.het=d.het,sd.het=sd.het,y='dens',
        vline=c(-d.crit,d0,d.het,d.crit),title=title);
  title=figtitle('Sampling distribution (cumulative)',model=mtext,
                 n=n,d0=d0,d.het=d.het,sd.het=sd.het,...);
  figname=paste(collapse='_',c('cum',mtext));
  dofig(plotpvsd,figname,n=n,d0=d0,d.het=d.het,sd.het=sd.het,y='cum',
        vline=c(-d.crit,d0,d.het,d.crit),title=title);
  return();
}
## wrapper for plotm to show concordance of simulated vs theoretical results
## available statistics: meand, pval, power
## for fixd model, do it by d0 (same as ovrfx)
## for hetd model, do it by sd.het for specific d.het (same as ovrht)
figm_simutheo=
  function(simu,theo,stat,mtext,by=cq(d0,d.het,sd.het),smooth='none',
           d0=NULL,d.het=NULL,sd.het=NULL,
           legend='topright',xlab='sample size',ylab=stat,...) {
    n=unique(simu$n);
    if (missing(by)&length(intersect(by,names(simu)))==0) by=NULL;
    if (!is.null(by)) {
      by=match.arg(by);
      simu.by=split(simu,simu[[by]]);
      theo.by=split(theo,theo[[by]]);
      ## grab stat colums from each group
      y=do.call(cbind,c(lapply(simu.by,function(simu) simu[,stat]),
                        lapply(theo.by,function(theo) theo[,stat])));
      legend.labels=
        c(paste0('simu. ',by,'=',names(simu.by)),paste0('theo. ',by,'=',names(theo.by)));
      ncol=length(simu.by);      
    } else {
      ## grab stat column from each source
      y=cbind(simu[,stat],theo[,stat]);
      legend.labels=cq(simu,theo);
      ncol=1;
    }
    if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none';
    col=setNames(RColorBrewer::brewer.pal(max(3,ncol),'Set1'),ncol);
    col=rep(col,2);
    lty=c(rep('dotted',len=ncol),rep('solid',len=ncol));
    lwd=2;
    title=figtitle(c('Line plot of simulated and theoretical',stat),
                   model=mtext,smooth=smooth,d0=d0,d.het=d.het,sd.het=sd.het);
    figname=paste(collapse='_',c(stat,mtext,smooth));
    dofig(plotm,figname,x=n,y=y,col=col,lty=lty,lwd=lwd,smooth=smooth,
          title=title,legend.labels=legend.labels,legend=legend,xlab=xlab,ylab=ylab,...);
    return();
  }
## wrapper for plotm to show coverage of confidence intervals - all & sig
## for fixd model, do it by d.pop
## for hetd model, do it by sd.het for d.het=0.5
figm_ci=function(ci,mtext,by=cq(d.pop,sd.het),smooth='none',
           d.pop=NULL,d.het=NULL,sd.het=NULL,
           legend='topright',xlab='sample size',ylab='coverage',...) {
  n=unique(ci$n);
  by=match.arg(by);
  ci.by=split(ci,ci[[by]]);
  ## grab desired columns from each group
  y=do.call(cbind,c(lapply(ci.by,function(ci) ci[,'cover']),
                    lapply(ci.by,function(ci) ci[,'cover.sig'])));
  legend.labels=
    c(paste0('cover. ',by,'=',names(ci.by)),paste0('cover.sig. ',by,'=',names(ci.by)));
  ncol=length(ci.by); 
  if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none';
  col=setNames(RColorBrewer::brewer.pal(max(3,ncol),'Set1'),ncol);
  col=rep(col,2);
  lty=c(rep('dotted',len=ncol),rep('solid',len=ncol));
  lwd=2;
  title=figtitle('Line plot of ci coverage',
                 model=mtext,smooth=smooth,d.pop=d.pop,d.het=d.het,sd.het=sd.het);
  figname=paste(collapse='_',c('ci',mtext,smooth));
  dofig(plotm,figname,x=n,y=y,col=col,lty=lty,lwd=lwd,smooth=smooth,
        title=title,legend.labels=legend.labels,legend=legend,xlab=xlab,ylab=ylab,...);
  return();
}

## plot sampling distributions (adapted from ovrfx figure 3)
## params hard-coded because no easy way to automate placement of text
plotsmpldist=function(dlim=c(-2,2),...) {
  title=figtitle('Sampling distributions show impact of increasing n and d');
  ## do n=100, d=0.2 first because it's way taller than the others
  n=c(100,20); d0=c(0.2,0.8); y=c(2.5,0.5); side=cq(left,right);
  d.crit=d_crit(n);
  cases=cbind(expand.grid(n=n,d0=d0),expand.grid(y=y,side=side,stringsAsFactors=FALSE));
  sapply(seq_len(nrow(cases)),function(i) with(cases[i,], {
    if (i==1) plotpvsd(n=n,d0=d0,dlim=dlim,add=F,vline=c(-d.crit,d.crit),title=title)
    else plotpvsd(n=n,d0=d0,dlim=dlim,add=T);
    dtext(n=n,d0=d0,y=y,side=side);
  }));
}
## helper function for placing text on sampling distributions (ovrfx figure 3)
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
## plot meand vs n, d (adapted from ovrfx figure 4)
plotmeand=function(col='black',lty='solid',lwd=1,legend.labels,meand20=NULL,nover=NULL,...) {
  title=figtitle('Average observed effect size improves as n increases');
  meand=get_data(meand.d2t);
  meand=subset(meand,subset=d0>0);
  n=unique(meand$n); d0=unique(meand$d0);
  meand.byd=split(meand,meand$d0);
  ## make splines that interpolate meand, over
  ##   empirically determined smooth.spline as interp function, and spar=0.3
  spar=0.3;
  meand.fit=lapply(meand.byd,function(meand) with(meand,smooth.spline(n,meand,spar=spar)));
  over.fit=lapply(meand.byd,function(meand) with(meand,smooth.spline(n,over,spar=spar)));
  nover=sapply(over.fit,function(fit) {
    n2over=function(n) predict(fit,n)$y;
    uniroot(function(n) n2over(n)-1.25,interval=c(10,200))$root});
  x=seq(floor(min(n,nover)),ceiling(max(n,nover)),by=1);
  y=do.call(cbind,lapply(meand.fit,function(fit) predict(fit,x)$y));
  col=setNames(RColorBrewer::brewer.pal(ncol(y),'Set1'),names(meand.byd));
  legend.labels=
    c(paste(sep=' ','mean sig effect for true effect',names(meand.byd)),
      paste(sep=' ','n achieving 1.25x for true effect',names(meand.byd)));
  col=rep(col,2);
  lty=c(rep('solid',len=3),rep('dotted',3));
  lwd=2;
  ## draw the main plot
  plotm(x=x,y=y,col=col,lty=lty,lwd=lwd,title=title,legend.labels=legend.labels,smooth=F,
        ylim=c(min(d0,y),max(d0,y)),xlab='sample size',ylab='effect size');
  ## horizontal grid-like lines for d0
  abline(h=d0,col=col,lty='dotted',lwd=1);
  ## vertical line for n=20
  vline(x=20,y0=0,y=1,col='grey',lty='dotted',lwd=1,text=20);
  ## horizontal & vertical lines for overestimates
  if (!is.null(nover)) {
    hline(x=nover,y=1.25*d0,col=col,lty='dotted',lwd=2,text=paste(sep='','1.25x',d0));
    vline(x=nover,y=1.25*d0,col=col,lty='dotted',lwd=2,text=round(nover));
  }
}
## histogram and d_d2t distribution (adapted from ovrht figure 1)
plothist_d2t=
  function(n=200,d=0,sd=0.2,xlim=d+c(-1,1),ylim=c(0,d_d2t(n=n,d=0)),
           col.hist='grey90',border.hist='grey80',
           ci.col='black',ci.lty='dotted',ci.lwd=2,...){ 
    title=figtitle('Het histogram and conventional sampling distribution under the null',
                   n=n,d.het=d,sd.het=sd);
    sim=get_sim_hetd(n=n,d=d,sd=sd);
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
## p-value inflation (adapted from ovrht figure 2)
plotpval_over=function(lwd=2,...) {
  title=figtitle('P-value inflation worsens as sd.het and n increase');
  pval=get_data(pval.d2ht);
  pval=subset(pval,subset=(sd.het!=0));
  pval.bysd=split(pval,pval$sd.het);
  n=unique(pval$n);
  y=do.call(cbind,lapply(pval.bysd,function(pval) pval$over.tval));
  col=setNames(RColorBrewer::brewer.pal(max(3,ncol(y)),'Set1'),names(pval.bysd));
  legend.labels=paste(sep='=','sd.het',names(pval.bysd));
  plotm(x=n,y=y,title=title,lwd=lwd,col=col,smooth='spline',
        legend='topleft',legend.labels=legend.labels,
        xlab='sample size',ylab='inflation (ratio of correct to conventional value)');
}
## confidence interval inflation (adapted from ovrht figure 3)
plotci_over=function(lwd=2,...) {
  title=figtitle('Confidence interval inflation worsens as sd.het and n increase',d.het=0);
  ci=get_data(ci.d2ht);
  ci=subset(ci,subset=(sd.het!=0&d==0));
  ci.bysd=split(ci,ci$sd.het);
  n=unique(ci$n);
  y=do.call(cbind,lapply(ci.bysd,function(ci) ci$over));
  col=setNames(RColorBrewer::brewer.pal(max(3,ncol(y)),'Set1'),names(ci.bysd));
  legend.labels=paste(sep='=','sd.het',names(ci.bysd));
  plotm(x=n,y=y,title=title,lwd=lwd,col=col,smooth='spline',
        legend='topleft',legend.labels=legend.labels,
        xlab='sample size',ylab='inflation (ratio of correct to conventional value)');
}
## sampling distributions (adapted from ovrht figure 4)
## mostly hard-coded because no easy way to automate placement of text
plotd2t_d2ht=
  function(n=c(20,200),d=0,sd.het=0.2,xlim=d+c(-1,1),ylim=c(0,d_d2t(n=n[2],d=0)),
           col.d2ht='grey70',...) {
    title=figtitle('Conventional and het sampling distributions for two values of n');
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
