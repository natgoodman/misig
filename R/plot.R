#################################################################################
##
## Author:  Nat Goodman
## Created: 19-01-09
##          uses code from repwr/R/plot.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Plotting code for effit documents
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- Plot Functions ----
## plot d.y vs d.x, typically d.pop vs d.sdz, colored by pval
## sim is data frame of simulation results
## params for generating pval colors
##   n is sample size - default from sim
##   sd.het is sd for noncentral d2ht distribution - default from sim if exists there
##   distribution is d2t or d2ht - default 'd2t' unless sd.het is set
## title is title
## cex.title is what R calls cex.main. if 'auto' or NULL, use auto-scaling
## legend tells whether to draw pval legend
##   legend.xscale, legend.yscale are fractions of plat ares
##   legend.cex is cex for legend text
## x, y are columns of sim containing x, y values; also used for axis labels
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
#
plotdvsd=
  function(sim,
           n=unique(sim$n),sd.het=if(exists('sd.het',sim))unique(sim$sd.het) else NULL,
           distribution=cq(d2t,d2ht),
           title='',cex.title='auto',x='d.sdz',y='d.pop',
           col=NULL,
           legend=TRUE,legend.xscale=1/8,legend.yscale=1/3,legend.x0=NULL,
           legend.label='p-value',legend.cex=0.75,
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',
           vhlwd=1,vlab=TRUE,hlab=TRUE,vhdigits=2,
           xlab=switch(x,
                       d.sdz="observed effect size (Cohen's d)",
                       d.pop="true effect size",
                       NULL),
           ylab=switch(y,
                       d.sdz="observed effect size (Cohen's d)",
                       d.pop="true effect size",
                       NULL),
           ...) {
    if (nrow(sim)==0) {
      warning('sim is empty. nothing to plot');
      return();
    }
    if (is.null(col)) {
      if (missing(distribution)) distribution=if(is.null(sd.het)) 'd2t' else 'd2ht'
      else distribution=match.arg(distribution);
      if (length(n)==1&length(sd.het)<=1)
        col=d2col(n=n,sd.het=sd.het,distribution=distribution,d=sim$d.sdz)
      else {
        if (length(n)!=1) warning('plotdvsd needs unique n to convert d to color');
        if (length(sd.het)>1) warning('plotdvsd needs unique sd.het to convert d to color');
        col='grey';
      }}
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    plot(sim[,x],sim[,y],col=col,main=title,cex.main=cex.title,xlab=xlab,ylab=ylab,pch=19,
         cex=0.5,...);
    grid();
    ## plot extra lines & values if desired. nop if vline, hline NULL
    vhline(vline=vline,hline=hline,vlab=vlab,hlab=hlab,vhdigits=vhdigits,
           lty=vhlty,col=vhcol,lwd=vhlwd);
    ## plot legend if desired
    if (legend)
      pval_legend(x.scale=legend.xscale,y.scale=legend.yscale,x0=legend.x0,label=legend.label,
                  cex=legend.cex);
  }
## plot histogram of (typically) d.sdz colored by pval
## sim is data frame of simulation results
## params for generating pval colors
##   n is sample size - default from sim
##   sd.het is sd for noncentral d2ht distribution - default from sim if exists there
##   distribution is d2t or d2ht - default 'd2t' unless sd.het is set
## title is title
## legend tells whether to draw pval legend
##   legend.x0 is harcoded x-position of legend - use with care!
##   legend.xscale, legend.yscale are fractions of plat areas
##   legend.cex is cex for legend text
## x is column of sim containing values for histogram; also used for axis labels
## breaks, freq are same-named parameters for hist
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
#
plothist=
  function(sim,
           n=unique(sim$n),sd.het=if(exists('sd.het',sim))unique(sim$sd.het) else NULL,
           distribution=cq(d2t,d2ht),
           title='',cex.title='auto',x='d.sdz',breaks=50,freq=FALSE,add=FALSE,
           col=NULL,border='black',
           legend=TRUE,legend.xscale=1/8,legend.yscale=1/3,legend.x0=NULL,
           legend.label='p-value',legend.cex=0.75,
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',
           vhlwd=1,vlab=TRUE,hlab=TRUE,vhdigits=2,
           xlab="observed effect size (Cohen's d)",ylab="density",
           ...) {
   if (nrow(sim)==0) {
      warning('sim is empty. nothing to plot');
      return();
   }
   hist.obj=hist(sim[,x],breaks=breaks,plot=FALSE);
   if (is.null(col)) {
     if (missing(distribution)) distribution=if(is.null(sd.het)) 'd2t' else 'd2ht'
     else distribution=match.arg(distribution);
     if (length(n)==1&length(sd.het)<=1)
       col=d2col(n=n,sd.het=sd.het,distribution=distribution,d=hist.obj$mids)
     else {
       if (length(n)!=1) warning('plothist needs unique n to convert d to color');
       if (length(sd.het)>1) warning('plothist needs unique sd.het to convert d to color');
       col='grey';
     }}
   if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    plot(hist.obj,col=col,border=border,freq=freq,main=title,cex.main=cex.title,
         xlab=xlab,ylab=ylab,add=add,...);
    if (!add) {
      box(); grid();
    } 
    ## plot extra lines & values if desired. nop if vline, hline NULL
   vhline(vline=vline,hline=hline,vlab=vlab,hlab=hlab,vhdigits=vhdigits,
          lty=vhlty,col=vhcol,lwd=vhlwd);
    ## plot legend if desired
   if (legend)
     pval_legend(x.scale=legend.xscale,y.scale=legend.yscale,x0=legend.x0,label=legend.label, cex=legend.cex);
  }
## plot probability distributions vs. d colored by pval
## n is sample size (for converting d to t or pval)
## d, col, lwd - usually missing - can be used when defaults not desired
## lwd.sig, lwd.nonsig - lwd for sig vs. non sig
##   if lwd set, use for both 
## dlim is range of d.sdz
##   dinc is step-size for extending d across range
## d0 is center for noncentral d2t distribution
## d.het is center for noncentral d2ht distribution
## sd.het is sd for noncentral d2ht distribution
## distribution is d2t or d2ht
## y is probability density or cumulative prob. can be values or keyword
## fill.tail tells whether to fill the distribution tail (density only)
##   boolean or one or more of 'upper','lower','both'. TRUE means c('upper','lower')
## add tells whether to add new plot to existing one
## xunit tells what to write along x-axis: d, t, or both
## title is title
## legend tells whether to draw pval legend
##   legend.xscale, legend.yscale are fractions of plat areas
##   legend.cex is cex for legend text
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
plotpvsd=
  function(n,d0=NULL,d.het=NULL,sd.het=NULL,distribution=cq(d2t,d2ht),y=cq(density,cumulative),
           d,col,
           lty.sig='solid',lty.nonsig=lty.sig,
           lwd=NULL,
           lwd.sig=if(!is.null(lwd)) lwd else 4,lwd.nonsig=if(!is.null(lwd)) lwd else lwd.sig/2,
           dlim=c(-2,2),dinc=.005,
           add=FALSE,fill.tail=FALSE,
           title='',cex.title='auto',d.crit=d_crit(n),
           legend=!add,legend.xscale=1/8,legend.yscale=1/3,legend.x0=NULL,
           legend.label='p-value',legend.cex=0.75,
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',
           vhlwd=1,vlab=TRUE,hlab=TRUE,vhdigits=2,
           xlab="observed effect size (Cohen's d)",ylab="probability",
           ...) {
    if (missing(d)) d=seq(min(dlim),max(dlim),by=dinc);
    if (missing(distribution)) {
      distribution=if(is.null(d.het)&is.null(sd.het)) 'd2t' else 'd2ht';
    }
    else distribution=match.arg(distribution);
    if (mode(y)=='character') {
      y=match.arg(y);
      if (y!='density') fill.tail=FALSE;
      if (missing(ylab)) ylab=if(y=='density') 'probability density' else 'cumulative probability';
      if (distribution=='d2t') {
        y=if(y=='density') d_d2t(n,d,d0) else p_d2t(n,d,d0);
      }
      else {
        y=if(y=='density') d_d2ht(n,d.het=d.het,sd.het=sd.het,d=d)
          else p_d2ht(n,d.het=d.het,sd.het=sd.het,d=d);
      }}
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    if (distribution=='d2t') pval=d2pval(n,d) else pval=d2htpval(n,sd.het,d);
    if (missing(col)||is.null(col)) col=pval2col(pval);
    lty=ifelse(pval<=0.05,lty.sig,lty.nonsig);
    lwd=ifelse(pval<=0.05,lwd.sig,lwd.nonsig);
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    l=length(d);
    if (!add) {
      plot(d,y,type='n',xlab=xlab,ylab=ylab,main=title,cex.main=cex.title,...);
      grid();
    }
    if (l==1) points(d,y,col=col,pch=19)
    else {
      x0=d[-l]; y0=y[-l]; x1=d[-1]; y1=y[-1];
      segments(x0,y0,x1,y1,col=col,lty=lty,lwd=lwd)
    }
    ## grid();
    ## plot extra lines & values if desired. nop if vline, hline NULL
    vhline(vline=vline,hline=hline,vlab=vlab,hlab=hlab,vhdigits=vhdigits,
           lty=vhlty,col=vhcol,lwd=vhlwd);
    ## fill tail if desired. already made sure we're doing density
    if (is.logical(fill.tail)&fill.tail) fill.tail=cq(upper,lower);
    if (!is.logical(fill.tail)) {
      if ('both' %in% fill.tail) fill.tail=cq(upper,lower);
      sapply(fill.tail,function(tail) fill_tail(tail,n,d,d0,d.crit));
    }
    ## plot legend if desired
    if (legend)
      pval_legend(x.scale=legend.xscale,y.scale=legend.yscale,x0=legend.x0,label=legend.label,
                  cex=legend.cex);
  }
## plot multiple lines - my adaptation of matplot - adapted from repwr/plotratm
## x is vector or matrix of x values
## y is vector or matrix of y values
##   like matplot, each line is column of x or y
##   unlike matplot, code currently assumes at most one of x,y is matrix - for smoothing to work
## col, lty, lwd are the usual line properties
## title, cex.title are title and cex for title
## xaxt, yaxt control labels on axes - NOT YET IMPLEMENTED
##   's' or NULL means let R do it
##    else list of axis params, eg, at, labels
##           title='',cex.title='auto'
## type is plot type - passed to matplot. 'n' also turns off extra lines, legend, grid
##   TODO: type='p' should cause legend to draw points instead of lines
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
## smooth whether to smooth data to make plot prettier
##   aspline, spline, loess, none, TRUE, FALSE. default is aspline.
##   TRUE means aspline. FALSE means none
##   spar is for spline
##   span is for loess
## smooth.xy tells which axis is domain of smoothing
##   default: 'x' if both x and y are vector-like, else whichever is vector-like
## legend tells whether and where to draw legend
##   TRUE or word like 'right' means draw legend, FALSE or NULL means no legend
## legend.title is legend title
## legend.args are further legend params
plotm=
  function(x,y,col='black',lty='solid',lwd=1,title='',cex.title='auto',type='l',
           xaxt='s',yaxt='s',
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',
           vhlwd=1,vlab=TRUE,hlab=TRUE,vhdigits=2,
           smooth=cq(aspline,spline,loess,linear,none),smooth.xy=cq(x,y),
           spar=NULL,span=0.75,
           legend=if(is.vector(y)) FALSE else 'right',legend.title=NULL,
           legend.labels=if(is.vector(y)) NULL else colnames(y),
           legend.args=list(where=NULL,x=NULL,y=NULL,cex=0.8,
                            title=legend.title,labels=legend.labels,col=col,lty=lty,lwd=lwd),
           ...) {
    if (is.null(x)) stop("Nothing to plot: 'x' is NULL");
    if (is.null(y)) stop("Nothing to plot: 'y' is NULL");
    if (is.vector(x)) x=as.matrix(x)
    else if (length(dim(x))!=2) stop("'x' must be vector or 2-dimensional matrix-like object");
    if (is.vector(y)) y=as.matrix(y)
    else if (length(dim(y))!=2) stop("'y' must be vector or 2-dimensional matrix-like object");
    if (nrow(x)!=nrow(y)) stop("'x' and 'y' have different number of rows");
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    ## NG 19-12-19: if smooth missing, code passes default (cq(..) in arg list) to smooth
    ##   function which chokes on it. Scary this bug wasn't caught sooner!
    ## smooth=if(is.logical(smooth)&&!smooth) 'none' else smooth;
    smooth.dflt='aspline';
    smooth=if(missing(smooth)) smooth.dflt 
           else {if(is.logical(smooth)) {if(smooth) smooth.dflt else 'none'}
                 else match.arg(smooth);}
    if (smooth!='none') {
      ## at most one of x,y can be matrix-like for smoothing to work
      if (dim(x)[2]>1&dim(y)[2]>1)
        stop("At most one of 'x' or 'y' can have multiple columns for smoothing to work");
      if (missing(smooth.xy)) smooth.xy=if(dim(x)[2]==1) 'x' else 'y';
      ## smooth
      if (smooth.xy=='x') {
        x.smooth=seq(min(x),max(x),len=100);
        y=smooth(x,y,xout=x.smooth,method=smooth,spar=spar,span=span);
        x=x.smooth;
      } else {
        y.smooth=seq(min(y),max(y),len=100);
        x=smooth(x=y,y=x,xout=y.smooth,method=smooth,spar=spar,span=span);
        y=y.smooth;
      }
    }
    matplot(x,y,main=title,cex.main=cex.title,col=col,lty=lty,lwd=lwd,type=type,
            xaxt=xaxt,yaxt=yaxt,...);
    if (type!='n') {
      grid();
      ## plot extra lines & values if desired. nop if vline, hline NULL
      vhline(vline=vline,hline=hline,vlab=vlab,hlab=hlab,vhdigits=vhdigits,
             lty=vhlty,col=vhcol,lwd=vhlwd);
      ## draw legend if desired
      if (is.null(legend)) legend=FALSE
      else if (!is.logical(legend)) {
        legend.args$where=legend;
        legend=TRUE;
      }
      if (legend) do.call(plotm_legend,legend.args);
    }
  }
## wrapper for plotm that pulls values from data frame
## x,y,z are column names
plotm_df=function(data,x,y,z=NULL,xlab=x,ylab=y,zlab=z,...) {
  force(xlab); force(ylab); force(zlab);
  if (!is.data.frame(data)) stop ('data must be data frame');
  if ((length(x)!=1)||(x %notin% colnames(data))) stop('x must contain exactly one column name');
  if ((length(y)!=1)||(y %notin% colnames(data))) stop('y must contain exactly one column name');
  if (!is.null(z)&&((length(z)!=1)||(z %notin% colnames(data))))
    stop('z must be NULL or contain exactly one column name');
  xdata=unique(data[,x,drop=F]);
  ydata=if(is.null(z)) data[,y,drop=F]
    else {
      by=split(data,data[,z]);
      do.call(cbind,lapply(by,function(data) data[,y,drop=F]))
    }
  plotm(xdata,ydata,xlab=xlab,ylab=ylab,legend.title=zlab,...);
}
## empty plot - just title & axes
plotempty=
  function(title='',cex.title='auto',xlab='x',ylab='y',xlim=c(0,1),ylim=c(0,1),
           xaxp=c(xlim,1),yaxp=c(ylim,1),...) {
    if (is.null(cex.title)|cex.title=='auto') cex.title=cex_title(title);
    plot(x=NULL,y=NULL,type='n',main=title,cex.main=cex.title,xlab=xlab,ylab=ylab,
         xlim=xlim,ylim=ylim,xaxp=xaxp,yaxp=yaxp,...);
    xylim=par('usr');                   # limits of disply region
    xmid=mean(xylim[1:2]);
    ymid=mean(xylim[3:4]);
    text(xmid,ymid,'PLOT DELIBERATELY LEFT BLANK',adj=c(0.5,0.5));
    invisible();
}

## helper functions to plot horizontal and vertical line segments
vhline=function(vline=NULL,hline=NULL,vlab=TRUE,hlab=TRUE,vhdigits=2,col=NA,cex=0.75,...) {
  xylim=par('usr');
  vline=vline[which(between(vline,xylim[1],xylim[2]))];
  hline=hline[which(between(hline,xylim[3],xylim[4]))];
  abline(v=vline,h=hline,col=col,...);
  ## write vhline values along axes
  vline=vline[vlab];
  if (length(vline)>0)
    mtext(round(vline,vhdigits),side=1,at=vline,col=col,line=0.25,cex=cex*par('cex'));
  hline=hline[hlab];
  if (length(hline)>0)
    mtext(round(hline,vhdigits),side=2,at=hline,col=col,line=0.25,cex=cex*par('cex'));
}
hline=
  function(y,x0=0,x,col='black',lty='solid',lwd=1,cex=0.75,text=NULL,
           label=list(text=text,side=2,at=y,col=col,line=0.25,cex=cex*par('cex'),las=1)) {
    segments(x0=x0,x1=x,y0=y,y1=y,col=col,lty=lty,lwd=lwd);
    if (!is.null(text)) do.call(mtext,label);
  }
vline=
  function(x,y0=0,y,col='black',lty='solid',lwd=1,cex=0.75,text=NULL,
           label=list(text=text,side=1,at=x,col=col,line=0.25,cex=cex*par('cex'),las=1)) {
    segments(x0=x,x1=x,y0=y0,y1=y,col=col,lty=lty,lwd=lwd);
    if (!is.null(text)) do.call(mtext,label);
  }

## draw pval legend. works for big picture figure and probability plots
pval_legend=function(x.scale,y.scale,x0=NULL,cex,label='p-value') {
  param(brk.pval,col.pval,steps.pvcol,sig.level);
  ## plt=par('usr');                       # plot region in user coordinates
  ## names(plt)=cq(left,right,bottom,top);
  xtkl=par('xaxp');                    # x tick locations
  ytkl=par('yaxp');                    # y tick locations
  names(xtkl)=names(ytkl)=cq(lo,hi,num);
  if (is.null(x0)) x0=xtkl['lo'];
  width=x.scale*(xtkl['hi']-x0);
  x1=x0+width;
  y1=ytkl['hi'];
  height=y.scale*(y1-ytkl['lo']);
  y0=y1-height;
  ## image sometimes leaves blank space when y0 is between tick marks. roundoff problem, I think
  ## works better to stretch y0 to next lower tick
  tkl=seq(ytkl['lo'],y1,len=ytkl['num']+1);
  y0=tkl[findInterval(y0,tkl)];
  height=y1-y0;                         # adjust height for new y0
  x=c(x0,x1);
  y=seq(y0,y1,length.out=2*steps.pvcol+1)[1:(2*steps.pvcol)]
  z=t(as.matrix(rev(head(brk.pval,-1))));
  image(x,y,z,add=TRUE,breaks=brk.pval,col=col.pval);
  ## add legend text
  x1=x1+strwidth(' ',cex=cex);
  text(x1,y0,0,adj=c(0,0),cex=cex)
  text(x1,y0+height/2,sig.level,adj=c(0,0.5),cex=cex)
  text(x1,y0+height,1,adj=c(0,1),cex=cex)
  y1=y1+strheight('p-value',cex=cex);
  text(x0+width/2,y1,label,adj=c(0.5,0.5),cex=cex);
}
## draw plotm legend. adapted from repwr/mesr_legend
## labels and legend are synonyms
plotm_legend=
  function(where=NULL,x=NULL,y=NULL,cex=0.8,bty='n',
           title=NULL,title.col='black',
           col='black',lty='solid',lwd=1,labels=NULL,legend=labels,...) {
    if (is.null(legend)) return();      # nothing to draw
    if (is.null(x)) x=where;
    legend(x,y,bty=bty,legend=legend,cex=cex,col=col,lwd=lwd,lty=lty,
          title=title,title.col=title.col,...);
  }
## plot multiple legends. adapted from repwr/ragm_legend
## legends is list of legend.args - arguments to base::legend or plotm_legend
## where, x, y are starting position
## others used as defaults in each legend
multi_legend=
  function(legends,legend='right',where=legend,x=NULL,y=NULL,cex=0.8,bty='n',
           title=NULL,title.col='black',col='black',lty='solid',lwd=1,
           ## include legend-related plotm args (except legend defined ablove)
           legend.title=NULL,legend.labels=NULL,labels=NULL,
           ...) {
    default.args=
      list(cex=cex,bty=bty,title=title,title.col=title.col,col=col,lty=lty,lwd=lwd,
           ## legend-related plotm args
           legend.title=legend.title,legend.labels=legend.labels,labels=labels,legend=NULL,
           ...);
    if (is.null(x)) x=where;
    sapply(legends,function(legend.args) {
      if (is.null(legend.args)) return();
      legend.args$x=x;
      legend.args$y=y;
      legend.args=fill_defaults(default.args,legend.args);
    ## handle legend-related plotm args
      ##   legend.title, legend.labels, labels, legend (when legend or labels not set)
      legend.args=within(legend.args,{
        if (is.null(title)) title=legend.title;
        if (is.null(legend)) {
          legend=if(is.null(labels)) legend.labels else labels;
        }
        rm(legend.title,legend.labels,labels);
      });
      if (is.null(legend.args$legend)) return();
      ## use 'graphics::legend' to avoid collision with legend arg
      where.next=do.call(graphics::legend,legend.args);
      ## <<- assigns to variables in outer scope, ie, function scope
      ##   from stackoverflow.com/questions/13640157. Thanks!
      ## could also use, eg, assign('x',where.next$rect$left,,envir=parent.frame(n=3))
      x<<-where.next$rect$left;
      y<<-where.next$rect$top-where.next$rect$h;
    })
    return();
  }

## fill tail of probability density
## adapted from https://stackoverflow.com/questions/45527077. Thx!
fill_tail=function(tail=cq(upper,lower),n,d,d0,d.crit) {
  tail=match.arg(tail);
  x=if(tail=='upper') x=d[d>d.crit] else x=d[d<(-d.crit)]
  y=d_d2t(n,x,d0);
  ## toss in extra x, y values to close polygon
  x=c(min(x),x,max(x));
  y=c(0,y,0);
  ## looks nicer if y stops just short of curve
  y=sapply(y,function(y) max(0,y-1e-3));
  ## do it!
  polygon(x=x,y=y,col='grey',border=NA);
}
## TODO: y.gap doesn't really do what I want...
fill_area=function(n,d0,d,col='grey',y.gap=1e-3) {
  y=d_d2t(n,d,d0);
  ## toss in extra x, y values to close polygon
  x=c(min(d),d,max(d));
  y=c(0,y,0);
  ## looks nicer if y stops just short of curve
  y=sapply(y,function(y) max(0,y-y.gap));
  ## do it!
  polygon(x=x,y=y,col=col,border=NA);
}

pval2col=function(pval) {
  param(col.pval,brk.pval,min.pvcol);
  col.pval[findInterval(-log10(clamp_pval(pval,min.pvcol)),brk.pval,all.inside=TRUE)];
}
d2col=function(n,sd.het,distribution,d) {
  pval=if(distribution=='d2t') d2pval(n,d) else d2htpval(n,sd.het,d);
  pval2col(pval);
}
clamp_pval=function(pval,min.pvcol) sapply(pval,function(pval) max(min(pval,1),min.pvcol));

## auto-scale title
cex_title=function(title) {
  xyplt=par('plt');                     # dimensions of plot region
  xplt=xyplt[2]-xyplt[1];               # width of plot region
  min(1,xplt/strwidth(title,units='fig'));
}
