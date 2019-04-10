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
## NG 19-01-09. presently just for siglo document
## plot d.y vs d.x, typically d.pop vs d.sdz, colored by pval
## sim is data frame of simulation results
## title is title
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
  function(sim,title='',cex.title=1,x='d.sdz',y='d.pop',
           legend=T,legend.xscale=1/8,legend.yscale=1/3,legend.cex=0.75,
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',vhlwd=1,vlab=T,hlab=T,vhdigits=2,
           xlab=switch(x,
                       d.sdz="observed effect size (Cohen's d)",
                       d.pop="true effect size",
                       NULL),
           ylab=switch(y,
                       d.sdz="observed effect size (Cohen's d)",
                       d.pop="true effect size",
                       NULL),
           ...) {
    plot(sim[,x],sim[,y],col=pval2col(sim[,'pval']),main=title,cex.main=cex.title,
         xlab=xlab,ylab=ylab,pch=19,cex=0.5,...);
    grid();
    ## plot extra lines if desired. nop if vline, hline NULL
    abline(v=vline,h=hline,lty=vhlty,col=vhcol,lwd=vhlwd);
    ## write values along axes
    if (vlab&!is.null(vline))
      mtext(round(vline,vhdigits),side=1,at=vline,col=vhcol,line=0.25,cex=0.75);
    if (hlab&!is.null(hline))
      mtext(round(hline,vhdigits),side=2,at=hline,col=vhcol,line=0.25,cex=0.75);
    ## plot legend if desired
    if (legend) pval_legend();
  }
## plot histogram of (typically) d.sdz colored by pval
## sim is data frame of simulation results
## n is sample size (for converting d to t or pval)
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
  function(sim,n=unique(sim$n),title='',cex.title=1,x='d.sdz',breaks=50,freq=F,add=F,
           col=NULL,border='black',
           legend=T,legend.x0=NULL,legend.xscale=1/8,legend.yscale=1/3,legend.cex=0.75,
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',vhlwd=1,vlab=T,hlab=T,vhdigits=2,
           xlab="observed effect size (Cohen's d)",ylab="density",
           ...) {
    hist.obj=hist(sim[,x],breaks=breaks,plot=F);
    if (is.null(col)) 
      if (length(n)==1) col=d2col(n=n,d=hist.obj$mids)
      else {
        warn('plothist need unique n to convert d to color');
        col='black';
      }
    plot(hist.obj,col=col,border=border,freq=freq,main=title,cex.main=cex.title,
         xlab=xlab,ylab=ylab,add=add,...);
    if (!add) {
      box(); grid();
    } 
    ## plot extra lines if desired. nop if vline, hline NULL
    abline(v=vline,h=hline,lty=vhlty,col=vhcol,lwd=vhlwd);
    ## write values along axes
    if (vlab&!is.null(vline))
      mtext(round(vline,vhdigits),side=1,at=vline,col=vhcol,line=0.25,cex=0.75);
    if (hlab&!is.null(hline))
      mtext(round(hline,vhdigits),side=2,at=hline,col=vhcol,line=0.25,cex=0.75);
    ## plot legend if desired
    if (legend) pval_legend(x0=legend.x0);
  }
## plot probability distributions vs. d colored by pval
## n is sample size (for converting d to t or pval)
## d, col, lwd - usually missing - can be used when defaults not desired
## lwd.sig, lwd.nonsig - lwd for sig vs. non sig
## dlim is range of d.sdz
##   dinc is step-size for extending d across range
## d0 is center for noncentral d2t distribution
## d.het is center for noncentral d2ht distribution
## sd.het is sd for noncentral d2ht distribution
## distribution is d2t or d2ht
## y is probability density or cumulative prob. can be values or keyword
## fill.tail tells whether to fill the distribution tail (density only)
##   boolean or one or more of 'upper','lower','both'. T means c('upper','lower')
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
           d,col,lwd,lwd.sig=4,lwd.nonsig=lwd.sig/2,dlim=c(-2,2),dinc=.005,add=F,fill.tail=F,
           title='',cex.title=1,d.crit=d_crit(n),
           legend=!add,legend.xscale=1/8,legend.yscale=1/3,legend.cex=0.75,
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',vhlwd=1,vlab=T,hlab=T,vhdigits=2,
           xlab="observed effect size (Cohen's d)",ylab="probability",
           ...) {
    if (missing(d)) d=seq(min(dlim),max(dlim),by=dinc);
    if (missing(distribution)) {
      distribution=if(missing(d.het)&missing(sd.het)) 'd2t' else 'd2ht';
    }
    else distribution=match.arg(distribution);
    if (mode(y)=='character') {
      y=match.arg(y);
      if (y!='density') fill.tail=F;
      if (missing(ylab)) ylab=if(y=='density') 'probability density' else 'cumulative probability';
      if (distribution=='d2t') {
        y=if(y=='density') d_d2t(n,d,d0) else p_d2t(n,d,d0);
      }
      else {
        y=if(y=='density') d_d2ht(n,d.het=d.het,sd.het=sd.het,d=d)
          else p_d2ht(n,d.het=d.het,sd.het=sd.het,d=d);
      }
    }
    if (distribution=='d2t') pval=d2pval(n,d) else pval=d2htpval(n,sd.het,d);
    if (missing(col)) col=pval2col(pval);
    if (missing(lwd)) lwd=ifelse(pval<=0.05,lwd.sig,lwd.nonsig);
    l=length(d);
    if (!add) plot(d,y,type='n',xlab=xlab,ylab=ylab,main=title,cex.main=cex.title,...);
    if (l==1) points(d,y,col=col,pch=19)
    else {
      x0=d[-l]; y0=y[-l]; x1=d[-1]; y1=y[-1];
      segments(x0,y0,x1,y1,col=col,lwd=lwd)
    }
    grid();
    ## plot extra lines if desired. nop if vline, hline NULL
    abline(v=vline,h=hline,lty=vhlty,col=vhcol,lwd=vhlwd);
    ## write values along axes
    if (vlab&!is.null(vline))
      mtext(round(vline,vhdigits),side=1,at=vline,col=vhcol,line=0.25,cex=0.75);
    if (hlab&!is.null(hline))
      mtext(round(hline,vhdigits),side=2,at=hline,col=vhcol,line=0.25,cex=0.75);
    ## fill tail if desired. already made sure we're doing density
    if (is.logical(fill.tail)&fill.tail) fill.tail=cq(upper,lower);
    if (!is.logical(fill.tail)) {
      if ('both' %in% fill.tail) fill.tail=cq(upper,lower);
      sapply(fill.tail,function(tail) fill_tail(tail));
    }
    ## plot legend if desired
    if (legend) pval_legend();
  }
## plot multiple lines - my adaptation of matplot - adapted from repwr/plotratm
## x is vector of x values
## y is vector or matrix of y values - like matplot, each line is column of y
## col, lty, lwd are the usual line properties
## title, cex.title are title and cex for title
## xaxt, yaxt control labels on axes - NOT YET IMPLEMENTED
##   's' or NULL means let R do it
##    else list of axis params, eg, at, labels
##           title='',cex.title=1,
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
## smooth whether to smooth data to make plot prettier
##   aspline, spline, loess, none, T, F. default is aspline. T means aspline. F means none
##   spar is for spline
## legend tells whether and where to draw legend
##   T or word like 'right' means draw legend, F or NULL means no legend
## legend.title is legend title
## legend.args are further legend params
plotm=
  function(x,y,col='black',lty='solid',lwd=1,title='',cex.title=1,
           xaxt='s',yaxt='s',
           vline=NULL,hline=NULL,vhlty='dashed',vhcol='grey50',vhlwd=1,vlab=T,hlab=T,vhdigits=2,
           smooth=c(cq(aspline,spline,loess,none),TRUE,FALSE),spar=NULL,
           legend=if(is.vector(y)) F else 'right',legend.title=NULL,
           legend.labels=if(is.vector(y)) NULL else colnames(y),
           legend.args=list(where=NULL,x=NULL,y=NULL,cex=0.8,
                            title=legend.title,labels=legend.labels,col=col,lty=lty,lwd=lwd),
           ...) {
    if (is.null(x)) stop("Nothing to plot: 'x' is NULL");
    if (is.vector(y)) {
      if (length(x)!=length(y)) stop("'x' and 'y' lengths differ");
      y=as.matrix(y);
    }
    else if (length(dim(y))!=2) {
      stop("'y' must be vector or 2-dimensional matrix-like object");
      if (length(x)!=nrow(y)) stop("'x' and 'y' have different number of rows");
    }
    smooth=if(is.logical(smooth)) if(smooth) 'aspline' else 'none' else match.arg(smooth);
    if (smooth!='none') {
      x.smooth=seq(min(x),max(x),len=100);
      if (smooth=='aspline') y=asplinem(x,y,xout=x.smooth,method='improved')
      else if (smooth=='spline') y=splinem(x,y,xout=x.smooth,spar=spar)
      else y=loessm(x,y,xout=x.smooth);
      x=x.smooth;
    }
    matplot(x,y,type='l',main=title,cex.main=cex.title,col=col,lty=lty,lwd=lwd,
            xaxt=xaxt,yaxt=yaxt,...);
    grid();
    ## plot extra lines if desired. nop if vline, hline NULL
    abline(v=vline,h=hline,lty=vhlty,col=vhcol,lwd=vhlwd);
    ## write values along axes
    if (vlab&!is.null(vline))
      mtext(round(vline,vhdigits),side=1,at=vline,col=vhcol,line=0.25,cex=0.75);
    if (hlab&!is.null(hline))
      mtext(round(hline,vhdigits),side=2,at=hline,col=vhcol,line=0.25,cex=0.75);
    ## draw legend if desired
    if (is.null(legend)) legend=F
    else if (!is.logical(legend)) {
      legend.args$where=legend;
      legend=TRUE;
    }
    if (legend) do.call(plotm_legend,legend.args);
  }
## helper functions to plot horizontal and vertical line segments
hline=
  function(y,x0=0,x,col='black',lty='solid',lwd=1,cex=0.75,text=NULL,
           label=list(text=text,side=2,at=y,col=col,line=0.25,cex=cex,las=1)) {
    segments(x0=x0,x1=x,y0=y,y1=y,col=col,lty=lty,lwd=lwd);
    if (!is.null(text)) do.call(mtext,label);
  }
vline=
  function(x,y0=0,y,col='black',lty='solid',lwd=1,cex=0.75,text=NULL,
           label=list(text=text,side=1,at=x,col=col,line=0.25,cex=cex,las=1)) {
    segments(x0=x,x1=x,y0=y0,y1=y,col=col,lty=lty,lwd=lwd);
    if (!is.null(text)) do.call(mtext,label);
  }

## draw pval legend. works for big picture figure and probability plots
pval_legend=
  function(x.scale=parent(legend.xscale,1/8),y.scale=parent(legend.yscale,1/3),x0=NULL,
           cex=parent(legend.cex,0.75)) {
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
  image(x,y,z,add=T,breaks=brk.pval,col=col.pval);
  ## add legend text
  x1=x1+strwidth(' ',cex=cex);
  text(x1,y0,0,adj=c(0,0),cex=cex)
  text(x1,y0+height/2,sig.level,adj=c(0,0.5),cex=cex)
  text(x1,y0+height,1,adj=c(0,1),cex=cex)
  y1=y1+strheight('p-value',cex=cex);
  text(x0+width/2,y1,"p-value",adj=c(0.5,0.5),cex=cex);
  }
## draw plotm legend. adapted from repwr/mess_legend
plotm_legend=
  function(where=NULL,x=NULL,y=NULL,cex=0.8,bty='n',
           title=NULL,title.col='black',
           col='black',lty='solid',lwd=1,labels=NULL) {
    if (is.null(labels)) return();      # nothing to draw
    if (is.null(x)) x=where;
    legend(x,y,bty=bty,legend=labels,cex=cex,col=col,lwd=lwd,lty=lty,
          title=title,title.col=title.col);
  }

## fill tail of probability density
## adapted from https://stackoverflow.com/questions/45527077. Thx!
fill_tail=
  function(tail=cq(upper,lower),n=parent(n),d=parent(d),d0=parent(d0),d.crit=parent(d.crit)) {
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

pval2col=function(pval) {
  param(col.pval,brk.pval,min.pvcol);
  col.pval[findInterval(-log10(clamp_pval(pval,min.pvcol)),brk.pval,all.inside=T)];
}
d2col=function(n,d) pval2col(d2pval(n=n,d=d));
  
clamp_pval=function(pval,min.pvcol) sapply(pval,function(pval) max(min(pval,1),min.pvcol));

