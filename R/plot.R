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
## legend TBD
## x, y are columns of sim containing x, y values; also used for axis labels
## d.crit is critical d.sdz separating nonsignificant and significant values
## d.pop is presumed true effect size
## xlim, ylim are x, y, limits. default is to use R's default
## vline,hline are vectors of x or y positions for extra vertical or horizontal lines
## vhlty, vhcol, vhlwd are lty, col, lwd for these extra lines
## vlab, hlab contol writing vline, hline values along axes
## vhdigits is number of digits for these values
#
plotdvsd=
  function(sim,title='',x='d.sdz',y='d.pop',d.crit=d_crit(20),d.pop=0.3,legend=T,
           vline=c(d.crit,d.pop),hline=c(d.pop),
           vhlty='dashed',vhcol='grey50',vhlwd=1,vlab=T,hlab=T,vhdigits=2,
           xlab="observed effect size (Cohen's d)",ylab="true effect size",
           ...) {
    plot(sim[,x],sim[,y],col=pval2col(sim[,'pval']),main=title,xlab=xlab,ylab=ylab,
         pch=19,cex=0.5,...);
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
## draw pval legend. specialized for big picture figure
## TODO: generalize
pval_legend=function(width=0.5,height=2,cex=0.75) {
  param(brk.pval,col.pval,steps.pvcol,sig.level);
  plt=par('usr');                       # plot region in user coordinates
  names(plt)=cq(x.left,x.right,y.bottom,y.top);
  x0=trunc(plt['x.left']);
  y0=trunc(plt['y.top'])-height;
  x=c(x0,x0+width);
  y=seq(y0,y0+height,length.out=2*steps.pvcol+1)[1:(2*steps.pvcol)]
  z=t(as.matrix(rev(head(brk.pval,-1))));
  image(x,y,z,add=T,breaks=brk.pval,col=col.pval);
  ## add legend text
  x1=x0+width+strwidth(' ',cex=cex);
  text(x1,y0,0,adj=c(0,0),cex=cex)
  text(x1,y0+height/2,sig.level,adj=c(0,0.5),cex=cex)
  text(x1,y0+height,1,adj=c(0,1),cex=cex)
  y1=max(y)+strheight('p-value',cex=cex);
  text(x0+width/2,y1,"p-value",adj=c(0.5,0.5),cex=cex);
}

pval2col=function(pval) {
  param(col.pval,brk.pval,min.pvcol);
  col.pval[findInterval(-log10(clamp_pval(pval,min.pvcol)),brk.pval,all.inside=T)];
}
clamp_pval=function(pval,min.pvcol) sapply(pval,function(pval) max(min(pval,1),min.pvcol));
