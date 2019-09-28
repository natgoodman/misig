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
## Generate figures and tables for ovrht supplement
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## --- Generate Figures and Tables for ovrht supplement ---

doc_ovrhtsupp=function(sect=NULL) {
  param(figextra);
  ## sect.all=cq(smpldist,plotpvsd,plotm,plotover);
  sect.all=c(paste(sep='.','cmp',cq(smpldist,meand,power,pval,ci)),
             paste(sep='.','ovr',cq(meand,power,pval,ci)));
  if (is.null(sect)) sect=sect.all else sect=pmatch_choice(sect,sect.all,start=FALSE);
  sapply(sect,function(sect) {
    sect_start(sect,sect.all);
##### plothist - really plothist overlaid with d2ht
    if (sect=='cmp.smpldist') {
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
    if (sect=='cmp.meand') {
      param(d.hetd,sd.hetd);
      simu=get_data(meand.hetd); theo=get_data(meand.d2ht);
      ## want 4 cases if extra, else 1 "middle" value
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,1);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,1);
      cases=expand.grid(stat=sall(meand),what=cq(raw,ovr),stringsAsFactors=FALSE);
      withrows(cases,case,{
        if (what=='ovr') {
          ## remove d.het==0 'cuz over is Inf
          simu=subset(simu,subset=d.het!=0); theo=subset(theo,subset=d.het!=0); 
        }
        ## do it by sd.het for each d.het
        sapply(d.het,function(d) {
          simu=subset(simu,subset=d.het==d); theo=subset(theo,subset=d.het==d); 
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                           'by sd.het'),d.het=d);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(dhet,d_pretty(d))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01);
          dofig(plotm_cmp,figname,simu=simu,theo=theo,stat=stat,by='sd.het',title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'),ylim=lim);
        });
        ## do it by d.het for each sd.het
        sapply(sd.het,function(sd) {
          simu=subset(simu,subset=sd.het==sd); theo=subset(theo,subset=sd.het==sd); 
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),'by d.het'),
                         sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
          dofig(plotm_cmp,figname,simu=simu,theo=theo,stat=stat,by='d.het',title=title,
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
        dofig(plotxy_cmp,figname,simu=simu,theo=theo,stat=stat,title=title);
      });
    }
    if (sect=='cmp.power') {
      param(d.hetd,sd.hetd);
      simu=get_data(power.hetd); theo=get_data(power.d2ht);
      ## want 4 cases if extra, else 1 "middle" value
      d.het=if(figextra) pick(d.hetd,4) else pick(d.hetd,1);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,1);
      cases=expand.grid(stat=sall(power),what=cq(raw,ovr),stringsAsFactors=FALSE);
      withrows(cases,case,{
        ## do it by sd.het for each d.het
        sapply(d.het,function(d) {
          simu=subset(simu,subset=d.het==d); theo=subset(theo,subset=d.het==d);
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                           'by sd.het'),d.het=d);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(dhet,d_pretty(d))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
          dofig(plotm_cmp,figname,simu=simu,theo=theo,stat=stat,by='sd.het',title=title,
                ylab=slab(stat,what,'(ratio of actual to conventional)'),legend='right',ylim=lim);
        });
        ## do it by d.het for each sd.het 
        sapply(sd.het,function(sd) {
          simu=subset(simu,subset=sd.het==sd); theo=subset(theo,subset=sd.het==sd);
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),'by d.het'),
                         sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
          dofig(plotm_cmp,figname,simu=simu,theo=theo,stat=stat,by='d.het',title=title,
                ylab=slab(stat,what,'(ratio of actual to conventional)'),legend='right',ylim=lim);
        })});
      withrows(cases,case,{
        ## dot plot entire dataset
        ## CAUTION: doesn't fit usual pattern of 4 figures per row. dunno if problem
        title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                         'across entire dataset'));
        figname=paste(collapse='.',c('plotxy',sfig(stat,what)));
        dofig(plotxy_cmp,figname,simu=simu,theo=theo,stat=stat,title=title,
              xlim=c(0,1),ylim=c(0,1));
      });
    }
    if (sect=='cmp.pval') {
      ## NOTE: pval doesn't have d.het, since always calculated under NULL
      param(sd.hetd,sig.dat);
      simu=get_data(pval.hetd); theo=get_data(pval.d2ht);
      cases=expand.grid(stat=sall(pval),what=cq(raw,ovr),stringsAsFactors=FALSE);
      withrows(cases,case,{
        ## do it by sd.het for each sig.level
        sapply(sig.dat,function(sig) {
          simu=subset(simu,subset=sig.level==sig); theo=subset(theo,subset=sig.level==sig);
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                           'by sd.het'),sig.level=sig);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sig,sig)));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.001)
          dofig(plotm_cmp,figname,simu=simu,theo=theo,stat=stat,by='sd.het',title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'),legend='topleft',ylim=lim);
        });
        ## do it by sig.level for each sd.het
        sapply(sd.het,function(sd) {
          simu=subset(simu,subset=sd.het==sd); theo=subset(theo,subset=sd.het==sd);
          title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                           'by sig.level'),
                         sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.001)
          dofig(plotm_cmp,figname,simu=simu,theo=theo,stat=stat,by='sig.level',title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'),legend='left',ylim=lim);
        })});
      withrows(cases,case,{
        ## dot-plot entire dataset
        title=figtitle(c('Concordance of simulated and theoretical',stit(stat,what),
                         'across entire dataset'));
        figname=paste(collapse='.',c('plotxy',sfig(stat,what)));
        lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.01)
        dofig(plotxy_cmp,figname,simu=simu,theo=theo,stat=stat,title=title,xlim=lim,ylim=lim,
              mergecol=cq(n,sd.het,sig.level));
      });
    }
    if (sect=='cmp.ci') {
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
        dofig(plotm_cmpci,figname,simu=simu,stat=stat,by='sd.het',title=title,ylim=lim);
      });
      ## do it by d.het for each sd.het 
      sapply(sd.het,function(sd) {
        simu=subset(simu,subset=sd.het==sd);
        title=figtitle('CI coverage (aka capture percentage) by d.het',sd.het=sd);
        figname=paste_nv(sdhet,sd_pretty(sd));
        lim=round_rng(range(simu[[stat]],theo[[stat]]),u=.1)
        dofig(plotm_cmpci,figname,simu=simu,stat=stat,by='d.het',title=title,ylim=lim);
      });
    }
    if (sect=='ovr.meand') {
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
          dofig(plotm_ovrbysd,figname,theo=theo,stat=swat(stat,what),title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'));
        });
        ## do it by d.het for each sd.het
        sapply(sd.het,function(sd) {
          theo=subset(theo,subset=(sd.het%in%c(0,sd))); 
          title=figtitle(c(stit(stat,what),'by d.het'),sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          dofig(plotm_ovrbyd,figname,theo=theo,stat=swat(stat,what),title=title,
                ylab=slab(stat,what,'(ratio of actual to correct)'));
        })});
    }
    if (sect=='ovr.power') {
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
          dofig(plotm_ovrbysd,figname,theo=theo,stat=swat(stat,what),title=title,legend='right',
                ylab=slab(stat,what,'(ratio of actual to conventional)'));
        });
        ## do it by d.het for each sd.het
        sapply(sd.het,function(sd) {
          theo=subset(theo,subset=(sd.het%in%c(0,sd))); 
          title=figtitle(c(stit(stat,what),'by d.het'),sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          dofig(plotm_ovrbyd,figname,theo=theo,stat=swat(stat,what),title=title,legend='topright',
                ylab=slab(stat,what,'(ratio of actual to conventional)'));
        })});
    }
    if (sect=='ovr.pval') {
      ## pval doesn't have d.het, since always calculated under NULL
      param(d.hetd,sd.hetd,sig.dat);
      theo=get_data(pval.d2ht);
      ## want 1 case for normal, 4 cases if extra
      sig.level=if(figextra) pick(sig.dat,4) else pick(sig.dat,1);
      sd.het=if(figextra) pick(sd.hetd,4) else pick(sd.hetd,1);
      stat='pval';
      sapply(cq(raw,ovr),function(what) {
        ## do it by sd.het for each sig.level
        sapply(sig.level,function(sig) {
          theo=subset(theo,subset=sig.level==sig); 
          title=figtitle(c(stit(stat,what),'by sd.het'),sig.level=sig);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(siglevel,sig)));
          dofig(plotm_ovrbysd,figname,theo=theo,stat=swat(stat,what),title=title,legend='topleft',
                ylab=slab(stat,what,'(ratio of actual to correct)'));
        });
        ## do it by sig.level for each sd.het
        sapply(sd.het,function(sd) {
          theo=subset(theo,subset=(sd.het%in%c(0,sd))); 
          title=figtitle(c(stit(stat,what),'by sig.level'),sd.het=sd);
          figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
          dofig(plotm_ovrbysig,figname,theo=theo,stat=swat(stat,what),title=title,
                legend='topleft',ylab=slab(stat,what,'(ratio of actual to correct)'));
        })});
    }
    if (sect=='ovr.ci') {
      param(d.hetd,sd.hetd,conf.level);
      theo=get_data(ci.d2ht);
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
        dofig(plotm_ovrbysd,figname,theo=theo,stat=swat(stat,what),title=title,legend='right',
              ylab=slab(stat,what,'(ratio of actual to conventional)'));
      });
      ## do it by d.het for each sd.het
      sapply(sd.het,function(sd) {
        theo=subset(theo,subset=(sd.het%in%c(0,sd))); 
        title=figtitle(c(stit(stat,what),'by d.het'),sd.het=sd);
        figname=paste(collapse='.',c(sfig(stat,what),paste_nv(sdhet,sd_pretty(sd))));
        dofig(plotm_ovrbyd,figname,theo=theo,stat=swat(stat,what),title=title,legend='topright',
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
      qq=do.call(rbind,withrows(cases,case,{
        x=q_d2ht(n=n,d.het=d.het,sd.het=sd.het,p=p);
        sim=get_sim_hetd(n=n,d=d.het,sd=sd.het);
        y=quantile(sim$d.sdz,probs=p);
        data.frame(x,y)}));
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
plotm_cmp=
  function(simu,theo,stat,by=cq(d.het,sd.het,sig.level),smooth='aspline',
           d.het=NULL,sd.het=NULL,title=NULL,
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
plotm_cmpci=
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
plotxy_cmp=
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
## wrapper for plotm to show stat (raw or inflation) with het and conventional pvals and baseline
## available statistics: meand, pval, power
plotm_ovrbysd=
  function(theo,stat,smooth='aspline',d.het=NULL,sd.het=NULL,title=NULL,
           legend='topright',xlab='sample size',
           ylab=paste(sep=' ',stat,'inflation (ratio of actual to correct)'),
           lty.bsln='dotted',lwd.bsln=1.5,...) {
    stat2=c(stat,paste0(stat,'.tval'));
    theo.bsln=split(theo,theo$sd.het==0);
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
plotm_ovrbyd=
  function(theo,stat,smooth='aspline',d.het=NULL,sd.het=NULL,title=NULL,
           legend='topright',xlab='sample size',
           ylab=paste(sep=' ',stat,'inflation (ratio of actual to correct)'),
           lty.bsln='dotted',lwd.bsln=1.5,...) {
    plotm_ovrby_notsd(by='d.het',,theo=theo,stat=stat,
                      smooth=smooth,title=title,legend=legend,xlab=xlab,ylab=ylab,
                      lty.bsln=lty.bsln,lwd.bsln=lwd.bsln,...);
    ## stat2=c(stat,paste0(stat,'.tval'));
    ## theo.byd=split(theo,theo$d.het);    
    ## x=unique(theo$n);
    ## theo.byd=lapply(theo.byd,function(theo) {
    ##   theo.bsln=split(theo,theo$sd.het==0);
    ##   bsln=theo.bsln[['TRUE']][,stat];
    ##   theo=theo.bsln[['FALSE']][,stat2]
    ##   cbind(theo,bsln);
    ## });
    ## ncol=length(theo.byd);      
    ## d.het=names(theo.byd);
    ## y=do.call(cbind,theo.byd);  
    ## col=RColorBrewer::brewer.pal(max(ncol,3),'Set1')[1:l];
    ## lty=c(cq(solid,dashed),lty.bsln);
    ## lwd=c(2,2,lwd.bsln);
    ## legend.args=list(labels=c('colors (represent d.het)',
    ##                           paste_nv('d.het'),
    ##                           'line types',
    ##                           'het p-values','conventional p-values','baseline'),
    ##                  col=c('black',col,'black',rep('grey20',3)),
    ##                  lty=c(NA,rep('solid',ncol),NA,lty),
    ##                  lwd=c(NA,rep(lwd,ncol),NA,lwd));
    ## col=rep(col,each=3);
    ## lty=rep(lty,ncol);
    ## lwd=rep(lwd,ncol);
    ## plotm(x=x,y=y,title=title,lwd=lwd,lty=lty,col=col,smooth=smooth,        
    ##       legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,...); 
  }
plotm_ovrbysig=
  function(theo,stat,smooth='aspline',title=NULL,
           legend='topright',xlab='sample size',
           ylab=paste(sep=' ',stat,'inflation (ratio of actual to correct)'),
           lty.bsln='dotted',lwd.bsln=1.5,...) {
    plotm_ovrby_notsd(by='sig.level',theo=theo,stat=stat,
                      smooth=smooth,title=title,legend=legend,xlab=xlab,ylab=ylab,
                      lty.bsln=lty.bsln,lwd.bsln=lwd.bsln,...);
    ## stat2=c(stat,paste0(stat,'.tval'));
    ## theo.bysig=split(theo,theo$sig.level);    
    ## x=unique(theo$n);
    ## theo.bysig=lapply(theo.bysig,function(theo) {
    ##   theo.bsln=split(theo,theo$sd.het==0);
    ##   bsln=theo.bsln[['TRUE']][,stat];
    ##   theo=theo.bsln[['FALSE']][,stat2]
    ##   cbind(theo,bsln);
    ## });
    ## ncol=length(theo.bysig);      
    ## sig.level=names(theo.bysig);
    ## y=do.call(cbind,theo.bysig);  
    ## col=RColorBrewer::brewer.pal(max(ncol,3),'Set1')[1:l];
    ## lty=c(cq(solid,dashed),lty.bsln);
    ## lwd=c(2,2,lwd.bsln);
    ## legend.args=list(labels=c('colors (represent sig.level)',
    ##                           paste_nv('sig.level'),
    ##                           'line types',
    ##                           'het p-values','conventional p-values','baseline'),
    ##                  col=c('black',col,'black',rep('grey20',3)),
    ##                  lty=c(NA,rep('solid',ncol),NA,lty),
    ##                  lwd=c(NA,rep(lwd,ncol),NA,lwd));
    ## col=rep(col,each=3);
    ## lty=rep(lty,ncol);
    ## lwd=rep(lwd,ncol);
    ## plotm(x=x,y=y,title=title,lwd=lwd,lty=lty,col=col,smooth=smooth,        
    ##       legend=legend,legend.args=legend.args,xlab=xlab,ylab=ylab,...); 
  }
plotm_ovrby_notsd=
  function(by,theo,stat,smooth='aspline',d.het=NULL,sd.het=NULL,title=NULL,
           legend='topright',xlab='sample size',
           ylab=paste(sep=' ',stat,'inflation (ratio of actual to correct)'),
           lty.bsln='dotted',lwd.bsln=1.5,...) {
    stat2=c(stat,paste0(stat,'.tval'));
    theo.by=split(theo,theo[[by]]);    
    x=unique(theo$n);
    theo.by=lapply(theo.by,function(theo) {
      theo.bsln=split(theo,theo$sd.het==0);
      bsln=theo.bsln[['TRUE']][,stat];
      theo=theo.bsln[['FALSE']][,stat2]
      cbind(theo,bsln);
    });
    ncol=length(theo.by);      
    by.val=names(theo.by);
    y=do.call(cbind,theo.by);  
    col=RColorBrewer::brewer.pal(max(ncol,3),'Set1')[1:l];
    lty=c(cq(solid,dashed),lty.bsln);
    lwd=c(2,2,lwd.bsln);
    legend.args=list(labels=c(paste0('colors (represent ',by,')'),
                              paste(sep='=',by,by.val),
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
