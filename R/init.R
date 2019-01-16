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
doc.all=cq(readme,siglo,sighi,supp);
init=function(
  ## doc parameters 
  doc='readme',                     # controls sim defaults, data, figure subdirs
  docx=match.arg(doc,doc.all), 
  ## simulation parameters
  n=switch(docx,                            # sample sizes
           readme=20*2^(0:4),               # 20,40,80,160,320 (5 values)
           siglo=20,                        #
           sighi=c(20,seq(50,500,by=50)),   # 20,50,100,...,1000 (11 values)
           xperiment=NA),                   # xperiment must supply values
  m=switch(docx,                            # number of populations
           xperiment=NA,                    # xperiment must supply values
           1e4),                            # others
  d.gen=runif,                              # function to generate population effect sizes
  d.args=list(n=m,min=-3,max=3),            # arguments passed to d.gen
  d=switch(docx,                            # population effect sizes
           xperiment=NA,                    # xperiment must supply values
           do.call(d.gen,d.args)),          # others call generating function
  ## analysis parameters
  ## TODO decide which are still needed
  sig.level=0.05,                   # for conventional significance
  conf.level=0.95,                  # for confidence intervals
  pred.level=0.95,                  # for prediction intervals
                                    # grid for various precacluated data
  ## program parameters, eg, for output files, error messages, etc.
  ## TODO decide which are still needed
  scriptname='effit',                      #
  mdir=paste_nv(m,m_pretty(m)),            # m subdirectory
  datadir=file.path('data',docx,mdir),     # directory for data files. default eg, data/repwr/m=1e4
  simdir=file.path(datadir,'sim'),         # directory for sim files
  ## NG 18-10-18: figdir, tbldir moved to init_doc
  ## figdir=file.path('figure',docx,mdir), # directory for figures. default eg, figure/repwr/m=1e4
  ## tbldir=file.path('table',docx,mdir), # directory for tables. default eg, table/repwr/m=1e4
  verbose=F,                               # print progress messages
  ## program control
  ## TODO decide which are still needed
  must.exist=F,                  # must all sub-inits succeed?
  load=NA,                       # shorthand for other load params
                                 #   NA means load if file exists
                                 #   T, F mean always or never load
  load.sim=load,                 # load saved simulations
  save=NA,                       # shorthand for other save params 
                                 #   NA means save unless file exists
                                 #   T, F mean always or never save
  save.sim=save,                 # save simulations (RData format)
  save.data=save,                # save top level results (RData & txt formats)
  save.txt=NA,                   # save results in txt format as well as RData
                                 #   NA means use default rule for type:
                                 #   F for all but top level data
  save.txt.sim=!is.na(save.txt)&save.txt,  # save txt simulations. default F
  save.txt.data=is.na(save.txt)|save.txt,  # save txt top level results. default T
  keep=NA,                       # shorthand for other keep params 
                                 #   NA means use default keep rule for type:
                                 #   T for all but detl
                                 #   T, F mean always or never keep
  keep.sim=!is.na(keep)&keep,    # keep simulations. default F
  keep.data=is.na(keep)|keep,    # keep top-level data. default T
                                 #    
  clean=F,                       # remove everything and start fresh
  clean.data=clean,              # remove datadir & in-memory cache
  clean.cache=T,                 # clean in-memory cache - always safe
  clean.sim=F,                   # clean simulations. default F
  end=NULL                       # placeholder for last parameter
  ) {
  doc=docx;                      # to avoid confusion later
  if (doc=='xperiment'&any(is.na(c(n,m,d))))
    stop('doc=xperiment but no value provided for n, m, or d');
  ## round d to avoid imprecise decimals
  d=round(d,digits=5); 
  ## assign parameters to param environment
  ## do it before calling any functions that rely on params
  init_param();
  ## clean and create output directories as needed
  if (clean.data) unlink(datadir,recursive=T);
  outdir=c(datadir,simdir);
  ## create data subdirectories. nop if already exist
  sapply(outdir,function(dir) dir.create(dir,recursive=TRUE,showWarnings=FALSE));
  ## setup in-memory cache to hold simulations, etc. do carefully in case already setup
  init_cache()
  ## clean specific types if desired. cleans directories and cache
  if (clean.sim) cleanq(sim);
  invisible();
}
## initialize doc parameters
## NG 19-01-11: abandon subdoc concept for 'supp' - not useful for effit
##              retain for xperiment just in case...
init_doc=function(
  subdoc=NULL,
  ## output directories. filename function ignores subdoc if NULL
  figdir=filename('figure',param(doc),subdoc,param(mdir)), # directory for figures
  tbldir=filename('table',param(doc),subdoc,param(mdir)),  # directory for tables
  ## output modifiers
  outpfx=switch(param(doc),supp='S',NULL),          # prefix before figure or table number
  outsfx=switch(param(doc),xperiment=NULL,letters), # suffix in figure and table blocks
  sectnum=switch(param(doc),supp=T,NULL),           # add section number to prefix eg, S1
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

## setup in-memory cache to hold simulations, etc. do carefully in case already setup
init_cache=function(memlist=cq(sim)) {
  param(clean.data,clean.cache);
  if (clean.data||clean.cache||!exists('cache.env',envir=.GlobalEnv)) cache.env=new.env();
  sapply(memlist,function(what)
    if (!exists(what,envir=cache.env,inherit=F)) assign(what,list(),envir=cache.env));
  assign('cache.env',cache.env,envir=.GlobalEnv);
}
## clean specific data type. deletes directory, any top level files and in-memory list
cleanq=function(what,cleandir=T) {
  what=as.character(pryr::subs(what));
  ## delete top level files if exist
  unlink(filename(datadir,list.files(datadir,pattern=paste(sep='','^',what,'\\.'))));
  ## delete in-memory list
  if (exists(what,envir=cache.env)) rm(list=what,envir=cache.env);
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
  
