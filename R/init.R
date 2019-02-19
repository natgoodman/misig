#################################################################################
##
## Author:  Nat Goodman
## Created: 19-01-01
##          from repwr/R/init.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Initialization code for effit
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- init ----
## initialization.
## process parameters and store in param environment.
## create output directories if necessary.
doc.all=cq(readme,ovrfx,ovrht,supp);
init=function(
  ## doc parameters 
  doc='readme',                     # controls sim defaults, data, figure subdirs
  docx=match.arg(doc,doc.all),
  run.id=NULL,                      # to separate runs for tests
  ## simulation parameters. random-d simulation
  n.rand=switch(docx,                     # sample sizes
                readme=c(20,100),
                ovrfx=20,
                ovrht=seq(20,200,by=20)),
  m.rand=switch(docx,                     # number of populations
                readme=1e3,
                ovrfx=1e5,
                ovrht=1e4),
  d.gen=switch(docx,                      # function to generate population effect sizes
               readme=runif,
               ovrfx=runif,
               ovrht=rnorm),
  d.args=switch(docx,                     # arguments passed to d.gen
                readme=list(n=m.rand,min=-3,max=3),
                ovrfx=list(n=m.rand,min=-3,max=3),
                ovrht=list(n=m.rand,mean=0.3,sd=0.1)),
  d.rand=do.call(d.gen,d.args),            # population effect sizes. call generating function

  ## simulation parameters. fixed-d simulation
  n.fixd=switch(docx,                       # sample sizes
                readme=seq(20,100,by=20),
                ovrfx=seq(20,200,by=20),
                ovrht=seq(20,200,by=20)), 
  m.fixd=switch(docx,                       # number of studies per d
                readme=1e2,
                ovrfx=1e4,
                ovrht=NA),
  d.fixd=switch(docx,                       # population effect sizes
                readme=c(0.2,0.5,0.8),      # Cohen's small, medium, large
                ovrfx=c(0.3,0.5,0.7),
                ovrht=NA),

  ## mean effect size computation (domeand)
  n.meand=n.fixd,
  d.meand=d.fixd,

  ## analysis parameters
  sig.level=0.05,                   # for conventional significance

  ## program parameters, eg, for output files, error messages, etc.
  datadir=dirname('data',docx,run.id),       # directory for data files
  sim.rand.dir=dirname(datadir,'sim.rand'),  # directory for sim random-d files
  sim.fixd.dir=dirname(datadir,'sim.fixd'),  # directory for sim fixed-d files
  ## NG 18-10-18: figdir, tbldir moved to init_doc
  ## figdir=dirname('figure',docx,mdir), # directory for figures. default eg, figure/repwr/m=1e4
  ## tbldir=dirname('table',docx,mdir), # directory for tables. default eg, table/repwr/m=1e4
  verbose=F,                               # print progress messages

  ## program control
  ## TODO decide which are still needed
  must.exist=F,                  # must all sub-inits succeed?
  save=NA,                       # shorthand for other save params 
                                 #   NA means save unless file exists
                                 #   T, F mean always or never save
  save.sim=save,                 # save simulations (RData format)
  save.meand=save,               # save mean effect size results (RData & txt formats)
  save.txt=NA,                   # save results in txt format as well as RData
                                 #   NA means use default rule for type:
                                 #   F for all but top level data
  save.txt.sim=!is.na(save.txt)&save.txt,  # save txt simulations. default F
  save.txt.meand=is.na(save.txt)|save.txt, # save txt meand results. default T
                                 #    
  clean=F,                       # remove everything and start fresh
  clean.data=clean,              # remove datadir
  clean.sim=F,                   # clean simulations. default F
  clean.meand=F,                 # clean mean effect size data. default F
  end=NULL                       # placeholder for last parameter
  ) {
  doc=docx;                      # to avoid confusion later

  ## extend d.rand to cover m
  d.rand=rep(d.rand,len=m.rand);
  ## round d.fixd, d.meand to avoid imprecise decimals 
  d.fixd=round(d.fixd,digits=5);
  d.meand=round(d.meand,digits=5);
  ## assign parameters to param environment
  ## do it before calling any functions that rely on params
  init_param();
  ## clean and create output directories as needed
  if (clean.data) unlink(datadir,recursive=T);
  outdir=c(datadir,sim.rand.dir,sim.fixd.dir);
  ## create data subdirectories. nop if already exist
  sapply(outdir,function(dir) dir.create(dir,recursive=TRUE,showWarnings=FALSE));
  ## clean specific types if desired. cleans directories and cache
  if (clean.sim) {
    unlink(sim.rand.dir,recursive=T);
    unlink(sim.fixd.dir,recursive=T);
  }
  if (clean.meand) cleanq(meand);
  invisible();
}
## initialize doc parameters
## NG 19-01-11: abandon subdoc concept for 'supp' - not useful for effit
##              retain for xperiment just in case...
init_doc=function(
  subdoc=NULL,
  ## output directories. filename function ignores subdoc if NULL
  figdir=dirname('figure',param(doc),subdoc,param(run.id)), # directory for figures
  tbldir=dirname('table',param(doc),subdoc,param(run.id)),  # directory for tables
  ## output modifiers
  outpfx=switch(param(doc),supp='S',NULL),          # prefix before figure or table number
  outsfx=letters,                                   # suffix in figure and table blocks
  sectpfx=switch(param(doc),readme=T,F),            # add section number to prefix eg, S1
  sectnum=1,                                        # section number. incremented in docs
  ## figures
  figpfx=outpfx,
  figsfx=outsfx,
  fignum=1,
  figblk=NULL,                  # index into figsfx if in figure block
  ## tables
  tblpfx=outpfx,
  tblsfx=outsfx,
  tblnum=1,
  tblblk=NULL,                  # index into tblsfx if in table block
  ## xtra figures - not included in document
  xfigpfx='X',
  xfigsfx=outsfx,
  xfignum=1,
  xfigblk=NULL,                 # index into xfigsfx if in figure block
  ## error cutoffs for plots
  ## TODO decide which are still needed
  fpr.cutoff=0.05,              # false positive rate cutoff for plots
  fnr.cutoff=0.20,              # false negative rate cutoff for plots
  ## for pval colors
  steps.pvcol=100,                # number of colors in color ramp
  min.pvcol=1e-4,                  # min pval in ramp - smaller pvals mapped to min
  ## clean, save
  save.out=T,
  save.fig=save.out,            # save figures (when called via dofig)
  save.tbl=save.out,            # save tables (when called via dotbl)
  save.txt.tbl=T,               # save txt tables. default T
  clean.out=F,
  clean.fig=clean.out,          # remove figdir
  clean.tbl=clean.out,          # remove tbldir
  ## plot control
  figscreen=if(param(doc)=='readme') T else !save.fig,
                                 # plot figures on screen
  fignew=figscreen,              # plot each figure in new window
  figextra=T,                    # plot extra figures, too
  ## doc generation function
  docfun=get(paste(collapse='',c('doc_',param(doc),subdoc))),
  docsect=NULL,                  # all document sections. set by docfun
  end=NULL                       # placeholder for last parameter
  ) {
  ## assign parameters to param environment
  ## do it before calling any functions that rely on params
  assign_param();
  ## initialize pval colors
  init_pvcol();
  ## clean and create output directories
  outdir=c(figdir,tbldir);
  if (clean.fig) unlink(figdir,recursive=T);
  if (clean.tbl) unlink(tbldir,recursive=T);
  sapply(outdir,function(dir) dir.create(dir,recursive=TRUE,showWarnings=FALSE));
  invisible();
}

## clean specific data type. deletes directory, and any top level files
cleanq=function(what,cleandir=T) {
  what=as.character(pryr::subs(what));
  ## delete top level files if exist
  unlink(filename(datadir,list.files(datadir,pattern=paste(sep='','^',what,'\\.'))));
  if (cleandir) {
    whatdir=paste(sep='',what,'dir');
    ## delete directory if exists
    if (exists(whatdir,envir=param.env)) unlink(get(whatdir,envir=param.env),recursive=T);
  }
}
## setup pval colors. adapted from repwr/R/plot.R/heat_setup
init_pvcol=function() {
  param(steps.pvcol,min.pvcol,sig.level);
  ## reds=colorRampPalette(RColorBrewer::brewer.pal(4,'Reds'))(steps.pvcol);
  ## blues=colorRampPalette(RColorBrewer::brewer.pal(5,'Blues'))(steps.pvcol);
  reds=colorRampPalette(RColorBrewer::brewer.pal(4,'Reds')[2:4])(steps.pvcol);
  blues=colorRampPalette(RColorBrewer::brewer.pal(4,'Blues')[2:4])(steps.pvcol);
  col.pval=c(blues,rev(reds));
  hi.brk=seq(0,-log10(sig.level),length.out=(steps.pvcol)+1);
  lo.brk=seq(-log10(sig.level),-log10(min.pvcol),length.out=(steps.pvcol)+1);
  brk.pval=unique(c(hi.brk,lo.brk))
  ## param.env$col.pval=col.pval;
  ## param.env$brk.pval=brk.pval;
  param(col.pval=col.pval,brk.pval=brk.pval);
}
  
