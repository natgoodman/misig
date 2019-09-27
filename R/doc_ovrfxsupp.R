#################################################################################
##
## Author:  Nat Goodman
## Created: 19-06-27
##          from doc_ovrhtsupp.R created 19-04-09
##          from doc_ovrht.R created 19-02-19
##          from doc_ovrfx.R created 19-02-03 
##          from doc_siglo.R created 19-01-10
##          from repwr/R/doc_resig.R created 18-09-05
## Includes code from repwr/R/docfun_resig.R created 18-10-25
##
## Copyright (C) 2019 Nat Goodman.
## 
## Generate figures and tables for ovrht supplement
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## TODO: still just copy of overhtsupp!!
## --- Generate Figures and Tables for ovrfx supplement ---

doc_ovrfxsupp=function(sect=NULL) {
  param(figextra);
  ## sect.all=cq(smpldist,plotpvsd,plotm,plotover);
  sect.all=c(paste(sep='.','simutheo',cq(smpldist,meand,power,pval,ci)),
             paste(sep='.','hetfix',cq(meand,power,pval,ci)));
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all,start=FALSE);
  sapply(sect,function(sect) {
    sect_start(sect,sect.all);
##### plothist - really plothist overlaid with d2ht
    if (sect=='simutheo.smpldist') {
      param(n.hetd,d.hetd,sd.hetd);
      ## start with 4 cases, then select 1st, last unless figextra. TODO: are these the best?
      cases=data.frame(n=rep(range(n.hetd),each=2),d.het=d.hetd[1:4],
                       sd.het=range(sd.hetd[sd.hetd!=0]));
      if (!figextra) cases=cases[c(1,4),];
      withrows(cases,case,{
        title=figtitle('Het histogram and sampling distribution',n=n,d.het=d.het,sd.het=sd.het);
        dofig(plothist_d2ht,'hist',n=n,d.het=d.het,sd.het=sd.het,title=title);
      });
      if (figextra)
        withrows(cases,case,{
          title=figtitle('Simulated vs. theoretical quantiles',n=n,d.het=d.het,sd.het=sd.het);
          dofig(plotqq_d2ht,'qq',n=n,d.het=d.het,sd.het=sd.het,title=title);
        });
      ## QQ plot entire dataset
      cases=expand.grid(n=n.hetd,d.het=d.hetd,sd.het=sd.hetd);
      title=figtitle('Simulated vs. theoretical quantiles across entire dataset');
      dofig(plotqq_d2ht,'qq',cases=cases,title=title);
    }
    if (sect=='simutheo.meand') {
      param(d.hetd,sd.hetd);
      simu=get_data(meand.hetd);
      theo=get_data(meand.d2ht);
      ## want 4 cases if extra, else 1 "middle" value
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,1);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,1);
      cases=expand.grid(stat=sall(meand),what=cq(raw,ovr),stringsAsFactors=FALSE);
      withrows(cases,case,{
        if (what=='ovr') {
          ## remove d.het==0 'cuz over is Inf
          simu=subset(simu,subset=d.het!=0);
          theo=subset(theo,subset=d.het!=0); 
        }
        ## do it by sd.het for each d.het
        sapply(d.het,function(d) {
          simu=subset(simu,subset=d.het==d);
          theo=subset(theo,subset=d.het==d); 
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                           'by sd.het'),d.het=d);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(dhet,d_pretty(d))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
          dofig(plotm_simutheo,figname,simu=simu,theo=theo,stat=stat,by='sd.het',title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'),ylim=lim);
        });
        ## do it by d.het for each sd.het
        sapply(sd.het,function(sd) {
          simu=subset(simu,subset=sd.het==sd);
          theo=subset(theo,subset=sd.het==sd); 
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),'by d.het'),
                         sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
          dofig(plotm_simutheo,figname,simu=simu,theo=theo,stat=stat,by='d.het',title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'),ylim=lim);
        })});
      withrows(cases,case,{
        ## dot plot entire dataset
        ## CAUTION: doesn't fit usual pattern of 4 figures per row. dunno if problem
        if (what=='ovr') {
          ## remove d.het==0 'cuz over is Inf
          simu=subset(simu,subset=d.het!=0);
          theo=subset(theo,subset=d.het!=0); 
        }
        title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                         'across entire dataset'));
        figname=paste(collapse='.',c('plotxy',sfig(stat,what)));
        dofig(plotxy_simutheo,figname,simu=simu,theo=theo,stat=stat,title=title);
      });
    }
    if (sect=='simutheo.power') {
      param(d.hetd,sd.hetd);
      simu=get_data(power.hetd);
      theo=get_data(power.d2ht);
      ## want 4 cases if extra, else 1 "middle" value
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,1);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,1);
      cases=expand.grid(stat=sall(power),what=cq(raw,ovr),stringsAsFactors=FALSE);
      withrows(cases,case,{
        ## do it by sd.het for each d.het
        sapply(d.het,function(d) {
          simu=subset(simu,subset=d.het==d);
          theo=subset(theo,subset=d.het==d);
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                           'by sd.het'),d.het=d);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(dhet,d_pretty(d))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
          dofig(plotm_simutheo,figname,simu=simu,theo=theo,stat=stat,by='sd.het',title=title,
                ylab=slab(stat,what,'(ratio of actual to conventional)'),legend='right',ylim=lim);
        });
        ## do it by d.het for each sd.het 
        sapply(sd.het,function(sd) {
          simu=subset(simu,subset=sd.het==sd);
          theo=subset(theo,subset=sd.het==sd);
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),'by d.het'),
                         sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
          dofig(plotm_simutheo,figname,simu=simu,theo=theo,stat=stat,by='d.het',title=title,
                ylab=slab(stat,what,'(ratio of actual to conventional)'),legend='right',ylim=lim);
        })});
      withrows(cases,case,{
        ## dot plot entire dataset
        ## CAUTION: doesn't fit usual pattern of 4 figures per row. dunno if problem
        title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                         'across entire dataset'));
        figname=paste(collapse='.',c('plotxy',sfig(stat,what)));
        dofig(plotxy_simutheo,figname,simu=simu,theo=theo,stat=stat,title=title,
              xlim=c(0,1),ylim=c(0,1));
      });
    }
    if (sect=='simutheo.pval') {
      simu=get_data(pval.hetd);
      theo=get_data(pval.d2ht);
      cases=expand.grid(stat=sall(pval),what=cq(raw,ovr),stringsAsFactors=FALSE);
      withrows(cases,case,{
        ## pval doesn't have d.het, since always calculated under NULL
        title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                         'by sd.het'),d.het=d);
        figname=sfig(stat,what);
        lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
        dofig(plotm_simutheo,figname,simu=simu,theo=theo,stat=stat,by='sd.het',title=title,
              ylab=slab(stat,what,'(ratio of actual to conventional)'),legend='topleft',ylim=lim);
      });
      withrows(cases,case,{
        ## dot-plot entire dataset
        title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                         'across entire dataset'));
        figname=paste(collapse='.',c('plotxy',sfig(stat,what)));
        lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
        dofig(plotxy_simutheo,figname,simu=simu,theo=theo,stat=stat,title=title,xlim=lim,ylim=lim,
              mergecol=cq(n,sd.het));
      });
    }
    if (sect=='simutheo.ci') {
      ## not much can be checked.
      ## 'simulated' data computes coverage (aka capture rate) for theoretical ci's
      param(conf.level);
      simu=get_data(ci.hetd);
      ## want 4 cases if extra, else 2 "middle" values
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,2);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,2);
      stat='cover';
      ## do it by sd.het for each d.het
      sapply(d.het,function(d) {
        simu=subset(simu,subset=d.het==d);
        title=figtitle('CI coverage (aka capture percentage) by sd.het',d.het=d);
        figname=paste_nv(dhet,d_pretty(d));
        lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.1)
        dofig(plotm_simutheoci,figname,simu=simu,stat=stat,by='sd.het',title=title,ylim=lim);
      });
      ## do it by d.het for each sd.het 
      sapply(sd.het,function(sd) {
        simu=subset(simu,subset=sd.het==sd);
        title=figtitle('CI coverage (aka capture percentage) by d.het',sd.het=sd);
        figname=paste_nv(sdhet,sd_pretty(sd));
        lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.1)
        dofig(plotm_simutheoci,figname,simu=simu,stat=stat,by='d.het',title=title,ylim=lim);
      });
    }
    if (sect=='hetfix.meand') {
      param(d.hetd,sd.hetd);
      theo=get_data(meand.d2ht);
      theo=subset(theo,subset=(d.het!=0)); # since over=Inf
      d.hetd=d.hetd[d.hetd!=0];            # since over=Inf
      sd.hetd=sd.hetd[sd.hetd!=0];         # since baseline
      ## want 1 case for normal, all cases (3 for now) if extra
      ## TODO: may need to pad with empty plot... to maintain 4 figs per row
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,1);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,1);
      stat='meand';
      sapply(cq(raw,ovr),function(what) {
        ## do it by sd.het for each d.het
        sapply(d.het,function(d) {
          theo=subset(theo,subset=d.het==d); 
          title=figtitle(c(stit(stat,what),'by sd.het'),d.het=d);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(dhet,d_pretty(d))));
          dofig(plotm_hetfixbysd,figname,theo=theo,stat=swat(stat,what),title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'));
        });
        ## do it by d.het for each sd.het
        sapply(sd.het,function(sd) {
          theo=subset(theo,subset=(sd.het%in%c(0,sd))); 
          title=figtitle(c(stit(stat,what),'by d.het'),sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          dofig(plotm_hetfixbyd,figname,theo=theo,stat=swat(stat,what),title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'));
        })});
    }
    if (sect=='hetfix.power') {
      param(d.hetd,sd.hetd);
      theo=get_data(power.d2ht);
      ## want 1 case for normal, 4 cases if extra
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,1);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,1);
      stat='power';
      sapply(cq(raw,ovr),function(what) {
        ## do it by sd.het for each d.het
        sapply(d.het,function(d) {
          theo=subset(theo,subset=d.het==d); 
          title=figtitle(c(stit(stat,what),'by sd.het'),d.het=d);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(dhet,d_pretty(d))));
          dofig(plotm_hetfixbysd,figname,theo=theo,stat=swat(stat,what),title=title,legend='right',
                ylab=slab(stat,what,'(ratio of actual to conventional)'));
        });
        ## do it by d.het for each sd.het
        sapply(sd.het,function(sd) {
          theo=subset(theo,subset=(sd.het%in%c(0,sd))); 
          title=figtitle(c(stit(stat,what),'by d.het'),sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          dofig(plotm_hetfixbyd,figname,theo=theo,stat=swat(stat,what),title=title,legend='topright',
                ylab=slab(stat,what,'(ratio of actual to conventional)'));
        })});
    }
    if (sect=='hetfix.pval') {
      ## pval doesn't have d.het, since always calculated under NULL
      param(d.hetd,sd.hetd);
      theo=get_data(pval.d2ht);
      ## want 1 case for normal, 4 cases if extra
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,1);
      stat='pval';
      sapply(cq(raw,ovr),function(what) {
        title=figtitle(c(stit(stat,what),'by sd.het'));
        figname=sfig(stat,what);
        dofig(plotm_hetfixbysd,figname,theo=theo,stat=swat(stat,what),title=title,legend='topleft',
              ylab=slab(stat,what,'(ratio of correct to conventional)'));
        });
    }
    if (sect=='hetfix.ci') {
      param(d.hetd,sd.hetd,conf.level);
      theo=get_data(ci.d2ht);
      ## TODO: why does ci use 'd' instead of d.het?
      ##       rename it d.het so code will work
      colnames(theo)[which(colnames(theo)=='d')]='d.het';
      ## want 1 case for normal, 4 cases if extra
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,1);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,1);
      stat='ci';
      what='raw';
      ## do it by sd.het for each d.het
      sapply(d.het,function(d) {
        theo=subset(theo,subset=d.het==d); 
        title=figtitle(c(stit(stat,what),'by sd.het'),d.het=d);
        figname=paste(collapse='.',c(sfig(stat,what),paste_nv(dhet,d_pretty(d))));
        dofig(plotm_hetfixbysd,figname,theo=theo,stat=swat(stat,what),title=title,legend='right',
              ylab=slab(stat,what,'(ratio of actual to conventional)'));
      });
      ## do it by d.het for each sd.het
      sapply(sd.het,function(sd) {
        theo=subset(theo,subset=(sd.het%in%c(0,sd))); 
        title=figtitle(c(stit(stat,what),'by d.het'),sd.het=sd);
        figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
        dofig(plotm_hetfixbyd,figname,theo=theo,stat=swat(stat,what),title=title,legend='topright',
              ylab=slab(stat,what,'(ratio of actual to conventional)'));
      });
      what='ovr';
      ## do it by sd.het for each d.het
      sapply(d.het,function(d) {
        theo=subset(theo,subset=d.het==d); 
        title=figtitle(c(stit(stat,what),'by sd.het'),d.het=d);
        figname=paste(collapse='.',c(sfig(stat,what),paste_nv(dhet,d_pretty(d))));
        dofig(plotm_hetfix,figname,theo=theo,stat='over',by='sd.het',title=title,legend='topleft',
              ylab=slab(stat,what,'(ratio of actual to conventional)'));
      });
      ## do it by d.het for each sd.het
      sapply(sd.het,function(sd) {
        theo=subset(theo,subset=(sd.het==sd)); 
        title=figtitle(c(stit(stat,what),'by d.het'),sd.het=sd);
        figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
        dofig(plotm_hetfix,figname,theo=theo,stat='over',by='d.het',title=title,legend='topleft',
                ylab=slab(stat,what,'(ratio of actual to conventional)'));
      });
    }
    return();
  })}

plothist_d2ht=
  function(n,d.het,sd.het,d.crit=NULL,xlim=d.het+c(-1.5,1.5),ylim=NULL,title=NULL,
           col.hist='grey90',border.hist='grey80',...){
    if (is.null(d.crit)) d.crit=d_htcrit(n,sd.het);
    if (is.null(ylim)) ylim=c(0,d_d2ht(n=n,d.het=d.het,sd.het=sd.het,d=d.het));
    sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
    plothist(sim=sim,col=col.hist,border=border.hist,title=title,
             xlim=xlim,ylim=ylim,ylab='probability density');
    if (is.null(d.crit)) d.crit=d_htcrit(n,sd.het);
    plotpvsd(n=n,d.het=d.het,sd.het=sd.het,add=T,dinc=1e-3,vline=c(-d.crit,d.crit),...)
    invisible();
  }
## if cases set, do entire dataset (ie, all cases) 
plotqq_d2ht=
  function(n,d.het,sd.het,cases,p=seq(0.01,0.99,by=.01),xlim=NULL,ylim=NULL,title=NULL,
           xlab='theoretical quantiles from het sampling distribution',ylab='simulated quantiles',
           cex=1,cex.cases=0.5,pch=19,col='black',
           col.diag='red',lty.diag='dotted',lwd.diag=1,...) {
    if (missing(cases)) {
      x=q_d2ht(n=n,d.het=d.het,sd.het=sd.het,p=p);
      sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
      y=quantile(sim$d.sdz,probs=p);
    } else {
      qq=withrows(cases,case,{
        x=q_d2ht(n=n,d.het=d.het,sd.het=sd.het,p=p);
        sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
        y=quantile(sim$d.sdz,probs=p);
        data.frame(x,y)})g;
      x=qq$x;
      y=qq$y;
      cex=cex.cases;
    }
    ## cor=cor(x,y);
    if (is.null(xlim)) xlim=range(x,y);
    if (is.null(ylim)) ylim=range(x,y);
    plot(x,y,xlim=xlim,ylim=ylim,main=title,cex.main=cex_title(title),
         xlab=xlab,ylab=ylab,pch=pch,cex=cex,col=col);
    abline(a=0,b=1,col=col.diag,lty=lty.diag,lwd=lwd.diag);
    grid();
  }
## wrapper for plotm to show concordance of simulated vs theoretical results
## available statistics: meand, pval, power
plotm_simutheo=
  function(simu,theo,stat,by=cq(d.het,sd.het),smooth='aspline',d.het=NULL,sd.het=NULL,title=NULL,
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
    if (length(n)==0&is.null(y)) plotempty(title=title,xlab=xlab,ylab=ylab)
    else {
      if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none';
      col=setNames(RColorBrewer::brewer.pal(max(3,ncol),'Set1'),ncol);
      col=rep(col,2);
      lty=c(rep('dotted',len=ncol),rep('solid',len=ncol));
      lwd=2;
      plotm(x=n,y=y,col=col,lty=lty,lwd=lwd,smooth=smooth,
            title=title,legend.labels=legend.labels,legend=legend,xlab=xlab,ylab=ylab,...);
    }
    return();
  }
## wrapper for plotm to show concordance of simulated ci with expectation
plotm_simutheoci=
  function(simu,stat,by=cq(d.het,sd.het),smooth='aspline',d.het=NULL,sd.het=NULL,title=NULL,
           legend='topright',xlab='sample size',ylab='CI coverage(aka capture rate)',...) {
    param(conf.level);
    n=unique(simu$n);
    if (missing(by)&length(intersect(by,names(simu)))==0) by=NULL;
    if (!is.null(by)) {
      by=match.arg(by);
      simu.by=split(simu,simu[[by]]);
      ## grab stat colums from each group
      y=do.call(cbind,lapply(simu.by,function(simu) simu[,stat]));
      legend.labels=paste0('simu. ',by,'=',names(simu.by));
      ncol=length(simu.by);      
    } else {
      ## grab stat column from each source
      y=cbind(simu[,stat]);
      legend.labels=cq(simu);
      ncol=1;
    }
    if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none';
    col=setNames(RColorBrewer::brewer.pal(max(3,ncol),'Set1'),ncol);
    lty=rep('dotted',len=ncol);
    lwd=2;
    plotm(x=n,y=y,col=col,lty=lty,lwd=lwd,smooth=smooth,
          hline=conf.level,vhcol='black',vhlty='solid',
          title=title,legend.labels=legend.labels,legend=legend,xlab=xlab,ylab=ylab,...);
    return();
  }
## dot plot to show concordance of simulated and theoretical results
plotxy_simutheo=
  function(simu,theo,stat,xlim=NULL,ylim=NULL,title=NULL,mergecol=cq(n,d.het,sd.het),
           xlab='theoretical values from het sampling distribution',ylab='simulated values',
           cex=1,pch=19,col='black',col.diag='red',lty.diag='dotted',lwd.diag=1,...) {
    xy=merge(simu,theo,by=mergecol,suffixes=cq(.simu,.theo));
    x=xy[[paste(sep='.',stat,'theo')]];
    y=xy[[paste(sep='.',stat,'simu')]];
    ## cor=cor(x,y);
    if (is.null(xlim)) xlim=range(x,y);
    if (is.null(ylim)) ylim=range(x,y);
    plot(x,y,xlim=xlim,ylim=ylim,main=title,cex.main=cex_title(title),
         xlab=xlab,ylab=ylab,pch=pch,cex=cex,col=col);
    abline(a=0,b=1,col=col.diag,lty=lty.diag,lwd=lwd.diag);
    grid();
  }
## wrapper for plotm to show inflation.
## only used for ci
plotm_hetfix=
  function(theo,stat,by=cq(d.het,sd.het),smooth='aspline',d.het=NULL,sd.het=NULL,title=NULL,
           legend='topright',xlab='sample size',
           ylab=paste(sep=' ',stat,'inflation (ratio of actual to correct)'),...) {
    n=unique(simu$n);
    if (missing(by)&length(intersect(by,names(simu)))==0) by=NULL;
    if (!is.null(by)) {
      by=match.arg(by);
      theo.by=split(theo,theo[[by]]);
      ## grab stat colums from each group
      y=do.call(cbind,lapply(theo.by,function(theo) theo[,stat]));
      legend.labels=c(paste0(stat,' ',by,'=',names(theo.by)));
      ncol=length(theo.by);      
    } else {
      ## grab stat column from each source
      y=theo[,stat];
      legend.labels=stat;
      ncol=1;
    }
    if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none';
    col=setNames(RColorBrewer::brewer.pal(max(3,ncol),'Set1'),ncol);
    lty='solid';
    lwd=2;
    plotm(x=n,y=y,col=col,lty=lty,lwd=lwd,smooth=smooth,
          title=title,legend.labels=legend.labels,legend=legend,xlab=xlab,ylab=ylab,...);
    return();
  }
## wrapper for plotm to show stat (raw or inflation) with het and conventional pvals and baseline
## available statistics: meand, pval, power
plotm_hetfixbysd=
  function(theo,stat,smooth='aspline',d.het=NULL,sd.het=NULL,title=NULL,
           legend='topright',xlab='sample size',
           ylab=paste(sep=' ',stat,'inflation (ratio of actual to correct)'),
           lty.bsln='dotted',lwd.bsln=1.5,...) {
    
    theo.bsln=split(theo,theo$sd.het==0);
    stat2=c(stat,paste0(stat,'.tval'));
    bsln=theo.bsln[['TRUE']][,stat];
    theo=theo.bsln[['FALSE']];

    theo.bysd=split(theo,theo$sd.het);
    theo.bysd=lapply(theo.bysd,function(theo) theo[,stat2])
    x=unique(theo$n);
    ## grab stat colums from each group
    y=cbind(do.call(cbind,theo.bysd),bsln);  

    ncol=length(theo.bysd);      
    sd.het=names(theo.bysd);
    
    if (is.logical(smooth)) smooth=if (smooth) 'aspline' else 'none';
    col=setNames(RColorBrewer::brewer.pal(max(3,ncol),'Set1'),ncol);
    lty=c(cq(solid,dashed),lty.bsln);
    legend.args=list(labels=c('colors (represent sd.het)',paste_nv('sd.het'),
                              'line types',
                              'het p-values','conventional p-values','baseline'),
                     col=c('black',col,'black',rep('grey20',3)),
                     lty=c(NA,rep('solid',ncol),NA,lty),
                     lwd=c(NA,rep(2,ncol),NA,2,2,lwd.bsln));
    col=c(rep(col,each=2),'black');
    lty=c(rep(cq(solid,dashed),ncol),lty.bsln);
    lwd=c(rep(c(2,2),ncol),lwd.bsln);
    plotm(x=x,y=y,title=title,lwd=lwd,lty=lty,col=col,smooth=smooth,       
          legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,...); 
  }
plotm_hetfixbyd=
  function(theo,stat,smooth='aspline',d.het=NULL,sd.het=NULL,title=NULL,
           legend='topright',xlab='sample size',
           ylab=paste(sep=' ',stat,'inflation (ratio of actual to correct)'),
           lty.bsln='dotted',lwd.bsln=1.5,...) {
    stat2=c(stat,paste0(stat,'.tval'));
    theo.byd=split(theo,theo$d.het);    
    x=unique(theo$n);
    theo.byd=lapply(theo.byd,function(theo) {
      theo.bsln=split(theo,theo$sd.het==0);
      bsln=theo.bsln[['TRUE']][,stat];
      theo=theo.bsln[['FALSE']][,stat2]
      cbind(theo,bsln);
    });
    ncol=length(theo.byd);      
    d.het=names(theo.byd);
    y=do.call(cbind,theo.byd);  
    col=RColorBrewer::brewer.pal(max(ncol,3),'Set1')[1:l];
    lty=c(cq(solid,dashed),lty.bsln);
    lwd=c(2,2,lwd.bsln);
    legend.args=list(labels=c('colors (represent d.het)',
                              paste_nv('d.het'),
                              'line types',
                              'het p-values','conventional p-values','baseline'),
                     col=c('black',col,'black',rep('grey20',3)),
                     lty=c(NA,rep('solid',ncol),NA,lty),
                     lwd=c(NA,rep(lwd,ncol),NA,lwd));
    col=rep(col,each=3);
    lty=rep(lty,ncol);
    lwd=rep(lwd,ncol);
    plotm(x=x,y=y,title=title,lwd=lwd,lty=lty,col=col,smooth=smooth,        
          legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,...); 
  }


## generate standard vector of stats for base, eg, 'meand'
sall=function(stat) {
  stat=as.character(pryr::subs(stat));
  c(stat,paste0(stat,'.tval'));
}
## for figname. shorten stat 
sfig=function(stat,what=cq(raw,ovr)) {
  what=match.arg(what);
  head=if (stat=='over'|what=='ovr') 'ovr' else NULL;
  tail=if (grepl('tval',stat)) 'tvl' else NULL;
  if (all(is.null(c(head,tail)))) NULL else paste(collapse='.',c(head,tail));
}
## for figure title. surround stat with ', append 'inflation' if what is 'ovr'
stit=function(stat,what=cq(raw,ovr)) {
  what=match.arg(what);
  stat=paste0("'",stat,"'");
  if (what=='raw') stat else c(stat,'inflation')
}
## for ylab. append 'inflation', desc if what is 'ovr'
slab=function(stat,what=cq(raw,ovr),desc=NULL) {
  what=match.arg(what);
  paste(collapse=' ',c(stat,if(what=='ovr') c('inflation',desc)));
}
## for plotting. get column name
swat=function(stat,what=cq(raw,ovr)) {
  what=match.arg(what);
  if (match.arg(what)=='raw') stat else 'over';
}
    
foo=function() {
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

  ## downsample to 1e4 so plotting will be fast
  m=param(m.hetd);
  if (m>1e4) sim=sim[sample.int(m,1e4),];
  ## draw the figures
  ## figure 1
  n.fig=200; d.fig=0.3; sd.fig=0.2;
  sim=get_sim_hetd(n=n.fig,d=d.fig,sd=sd.fig);
  title=title_ovrhtsupp('P-values improve as observed effect size grows more extreme',
                    n=n.fig,d.het=d.fig,sd.het=sd.fig);
  x='d.sdz'; y='d.pop';
  d.crit=d_crit(n.fig); d.htcrit=d_htcrit(n.fig,sd.fig);
  vl=c(d.fig,-d.crit,d.crit,-d.htcrit,d.htcrit);
  hl=d.fig;
  xlim=c(-1,1.5);
  dofig(plotdvsd,'big_picture',sim=sim,x=x,y=y,vline=vl,hline=hl,xlim=xlim,title=title);
  ## figure 2
  title=title_ovrhtsupp('Histogram of observed effect size',n=n.fig,d.het=d.fig,sd.het=sd.fig1);
  ylim=c(0,d_d2t(n=n.fig,d=d.fig,d0=d.fig));  # set ylim to match figure 3
  dofig(plothist,'hist',sim=sim,vline=vl,title=title,xlim=xlim,ylim=ylim);
  ## figure 3
  figblk_start();
  n.fig=200; d.fig=c(0.3,0.3,0.3,0.7,0.7); sd.fig=c(0,0.1,0.2,0.1,0.2);
  text.fig=list(label=paste_nv('sd.het',sd.fig),
                y=c(3.5,2.25,0.75,2.25,0.75),side=cq(left,left,left,right,right));
  title=title_ovrhtsupp('Sampling distributions with varying sd.het and d.het',n=n.fig);
  dofig(plotsmpldist,'smpldist',n=n.fig,d.het=d.fig,sd.het=sd.fig,
        text=text.fig,vline=c(-d.htcrit,d.fig,d.htcrit),title=title,xlim=xlim,ylim=ylim);
  ##
  n.fig=c(200,20,20); d.fig=c(0.3,0.3,0.7); sd.fig=0.2;
  text.fig=list(label=paste_nv('n',n.fig),y=c(1.5,0.5,0.5),side=cq(left,left,right));
  title=title_ovrhtsupp('Sampling distributions with varying n and d.het',sd.het=sd.fig);
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
## generate title for doc_ovrhtsupp
title_ovrhtsupp=function(...,sep=' ',n=NULL,d.het=NULL,sd.het=NULL) {
  fig=paste(sep='','Figure ',figlabel());
  desc=paste(sep=sep,...);
  nv=nvq(sep=', ',IGNORE=T,n,d.het,sd.het);
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
           title=NULL,xlab='sample size',ylab=NULL,legend='topright',...) {
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
    ## title=title_ovrhtsupp('Mean significant effect size vs. n, d',sd.het=sd.fig);
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
    plotm(x=x,y=y,title=title,lwd=lwd,lty=lty,col=col,        
          legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,smooth=F,...); 
  }
title_byd=function(what,y,sd.het) {
  desc=if (what=='meand') 'Mean significant effect size' else 'Power';
  if (y=='over') desc=paste(sep=' ',desc,'inflation');
  title_ovrhtsupp(desc,'vs. n, d.het',sd.het=sd.het);
}
ylab_byd=function(what,y) {
  ylab=if (what=='meand') 'effect size' else 'power';
  if (y=='over') ylab=paste(sep=' ',ylab,'inflation');
  ylab;
}
## plot meand, power (raw or inflation) split by sd.het
plot_bysd=
  function(what=cq(meand,power),y=what,d.fig=0.3,sd.fig=NULL,
           title=NULL,xlab='sample size',ylab=NULL,legend='topright',...) {
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
    plotm(x=x,y=y,title=title,lwd=lwd,lty=lty,col=col,        
          legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,smooth=F,...); 
  }
title_bysd=function(what,y,d.het) {
  desc=if (what=='meand') 'Mean significant effect size' else 'Power';
  if (y=='over') desc=paste(sep=' ',desc,'inflation');
  title_ovrhtsupp(desc,'vs. n, sd.het',d.het=d.het);
}
ylab_bysd=function(what,y) {
  ylab=if (what=='meand') 'effect size' else 'power';
  if (y=='over') ylab=paste(sep=' ',ylab,'inflation');
  ylab;
}

