#################################################################################
##
## Author:  Nat Goodman
## Created: 19-01-01
##          from repwr/R/datman.R created 18-05-03
##
## Copyright (C) 201 Nat Goodman.
## 
## Data management -- files and cached objects -- for repwr.R: 
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- Save and Load ----
## call with file or software attempts to construct file name
##### for data keyed by n
## save data in RData and optionally txt formats
save_n=function(data,n,file=NULL,what,save,save.txt=F,keep) {
  if (is.null(file)) base=basename_n(n,what)
  else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    save(data,file=file);
    if (save.txt)
      write.table(data,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
  }
  if (keep) keep_n(data,n,what=what);
}
## load data from file
load_n=function(file=NULL,n=NULL,what) {
  if (is.null(file)) file=filename_n(n,what);
  what=load(file=file);               # what is name of saved data
  get(what);                          # return it
}
## get data already in memory (eg, in sim.list) or read from file
##   fail if data does not exist unless must.exist is FALSE
get_n=function(n=NULL,what,load,keep,must.exist=T) {
  case=casename_n(n,short=T);
  if (is.na(load)|load) {
    what.list=get(what,envir=cache.env);
    data=what.list[[case]];
    if (is.null(data)) {
      file=filename_n(n,what);
      if (file.exists(file)) data=load_n(file=file)
      else {
        if (must.exist)
          stop(paste(sep=' ',what,'for case',case,
                     'not in cache and file',file,'does not exist'));
        data=NULL;
      }
    if (keep) keep_n(data,n,what=what);
    }}
    else {
      if (must.exist)
        stop(paste(sep=' ','must.exist is TRUE but load is FALSE for',what,'case',case));
      data=NULL;
    }
  invisible(data);
}
keep_n=function(data,n=NULL,what) {
  case=casename_n(n,short=T);
  ## CAUTION: have to do get/assign in single statement to update env. sigh...
  get(what,envir=cache.env)[[case]]=data;
}

##### sim
save_sim=function(sim,n,file=NULL) {
  param(save.sim,save.txt.sim,keep.sim);
  save_n(sim,n,file,what='sim',save=save.sim,save.txt=save.txt.sim,keep=keep.sim);
  }
load_sim=function(file=NULL,n) load_n(file,n,what='sim');
get_sim=function(n,load=load.sim,keep=keep.sim,must.exist=T) {
  param(load.sim,keep.sim);
  get_n(n,what='sim',load=load.sim,keep=keep.sim,must.exist);
}

##### table - saved in tbldir
save_tbl=function(tbl,file) {
  param(save.tbl,save.txt.tbl);
  if (is.null(tbl)) stop('Trying to save NULL table. Is table name set correctly?');
  base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save.tbl)&!file.exists(file))|(!is.na(save.tbl)&save.tbl)) {
    save(tbl,file=file);
    if (save.txt.tbl) {
      file=filename(base=base,suffix='txt');
      if (length(dim(tbl))==2) write.table(tbl,file=file,sep='\t',quote=F,row.names=F)
      else if (is.vector(tbl)) writeLines(as.character(tbl),file)
      else stop('Trying to save object with more than 2 dimensions as text. Is table name set correctly?');
    }}
  invisible(tbl);
}

## ---- File Functions ----
##### names for data keyed by n
filename_n=function(n,what,suffix='RData')
  filename(basename_n(n,what),suffix=suffix);
basename_n=function(n,what)
  filename(param(list=paste(sep='',what,'dir')),base=what,tail=casename_n(n))
casename_n=function(n=NULL,short=F) if (short) as.character(n) else paste_nv(n);

##### sim
filename_sim=function(n=NULL,suffix='RData') filename(basename_sim(n),suffix=suffix);
basename_sim=function(n=NULL) {
  param(simdir);
  filename(simdir,base='sim',tail=casename_n(n))
}
casename_sim=casename_n;

## figure and table functions
filename_fig=function(figlabel,sect,figname,suffix='png')
  filename(param(figdir),paste(collapse='_',c('figure',figlabel,sect,figname)),suffix=suffix);
filename_tbl=function(tbllabel,sect,tblname,suffix='RData')
  filename(param(tbldir),paste(collapse='_',c('table',tbllabel,sect,tblname)),suffix=suffix);

## construct file or directory pathname from components
## wrapper for file.path with base, tail and suffix pasted on
##  base appended with '.'
##  tail components combined with '.'
##  suffix added unless already there
filename=function(...,base=NULL,tail=NULL,suffix=NULL) {
  if (!is.null(base)||!is.null(tail)) base=paste(collapse='.',c(base,tail));
  ## NG 18-10-15: remove NULL from ... before calling file.path
  ## do.call(f,as.list(unlist(list(...))))) from https://stackoverflow.com/questions/47360937/call-an-r-function-with-run-time-generated-ellipsis-arguments-dot-dot-dot-thr
  ##  if (is.null(base)) file=file.path(...) else file=file.path(...,base);
  file=do.call(file.path,as.list(unlist(list(...))));
  if (!is.null(base)) file=file.path(...,base);
  if (!is.null(suffix)) {
    ## remove leading '.' if present
    suffix=sub('^\\.','',suffix,perl=T);
    suffix.pattern=paste(collapse='|',paste(sep='','\\.',suffix,'$'));
    file=ifelse(grepl(suffix.pattern,file),file,paste(sep='.',file,suffix[1]));
  }
  file;
}
## remove suffix from filename
desuffix=function(file,suffix=c('RData','txt')) {
  if (!is.null(suffix)) {
    ## remove leading '.' if present
    suffix=sub('^\\.','',suffix,perl=T);
    suffix.pattern=paste(collapse='|',paste(sep='','\\.',suffix,'$'));
    file=sub(suffix.pattern,'',file);
  }
  file;
}
## filebasename same as filename but w/o suffix
filebasename=function(...) filename(...,suffix=NULL)
## construct directory pathname. synonym for filebasename
dirname=filebasename;

