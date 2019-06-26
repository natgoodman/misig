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
## save data in RData and optionally txt formats
save_=function(data,file,save,save.txt=F) {
  if ((is.na(save)&!file.exists(file))|(!is.na(save)&save)) {
    base=desuffix(file);
    save(data,file=filename(base=base,suffix='RData'));
    if (save.txt)
      write.table(data,file=filename(base=base,suffix='txt'),sep='\t',quote=F,row.names=F);
  }
}
## load data from RData file
load_=function(file,what) {
  base=desuffix(file);
  what=load(file=filename(base=base,suffix='RData')); # what is name of saved data
  get(what);                          # return it
}
##### data - top-level data saved in datadir
## base does not include path
save_data=function(what,file=NULL,data=NULL,base=NULL,
                   save=param(save.data),save.txt=param(save.txt)) {
  param(datadir,save.data,save.txt.data);
  what=as.character(pryr::subs(what));
  if (missing(data) && exists(what,envir=parent.frame(n=1)))
    data=get(what,envir=parent.frame(n=1));
  if (is.null(data)) stop('Trying to save NULL object. Is "what" set correctly?');
  if (is.null(file)) {
    if (is.null(base)) base=basename_data(what) else base=basename_data(base);
  } else base=desuffix(file);
  file=filename(base=base,suffix='RData');
  save_(data,file=file,save=save.data,save.txt=save.txt.data);
}
## load data from file
load_data=function(...,file=NULL,base=NULL,list=character()) {
  dots=match.call(expand.dots=FALSE)$...;
  if (length(dots) &&
      !all(vapply(dots,function(x) is.atomic(x)||is.symbol(x)||is.character(x),
                  NA,USE.NAMES=FALSE))) 
    stop("... must contain atomic data like names or character strings");
  parent.env=parent.frame(n=1);
  names=vapply(dots,as.character,"");
  if (length(names)==0L) names=character();
  names=c(list,names);
  if (length(names)==0) {
    ## easy case - load and return one file
    if (is.null(file)) stop('Cannot load data unless file or what is set');
    return(load_(file=file));
  }
  val=lapply(1:length(names),function(i) {
    what=names[i];
    if (length(file)>=i) file=file[i] else file=filename_data(what);
    val=load_(file=file);
    assign(what,val,envir=parent.env);
    val;
  });
  if (length(names)==1) return(val[[1]]) else names;
}
get_data=load_data;
##### sim_rand
save_sim_rand=function(sim,n,file=NULL) {
  param(save.sim,save.txt.sim);
  if (is.null(file)) file=filename(param(sim.rand.dir),base='sim',tail=nvq(n));
  save_(sim,file,save=save.sim,save.txt=save.txt.sim);
}
load_sim_rand=function(file=NULL,n) {
  if (is.null(file)) file=filename(param(sim.rand.dir),base='sim',tail=nvq(n));
  load_(file,'sim');
}
get_sim_rand=load_sim_rand;
##### sim_fixd
save_sim_fixd=function(sim,n,d,file=NULL) {
  param(save.sim,save.txt.sim);
  if (is.null(file))
    file=filename(param(sim.fixd.dir),base='sim',
                  tail=paste(sep=',',paste_nv(n),paste_nv(d,d_pretty(d))));
  save_(sim,file,save=save.sim,save.txt=save.txt.sim);
}
load_sim_fixd=function(file=NULL,n,d) {
  if (is.null(file))
    file=filename(param(sim.fixd.dir),base='sim',
                  tail=paste(sep=',',paste_nv(n),paste_nv(d,d_pretty(d))));
  load_(file,'sim');
}
get_sim_fixd=load_sim_fixd;

##### sim_hetd
save_sim_hetd=function(sim,n,d,sd,file=NULL) {
  param(save.sim,save.txt.sim);
  if (is.null(file))
    file=filename(param(sim.hetd.dir),base='sim',
                  tail=paste(sep=',',paste_nv(n),
                             paste_nv(d,d_pretty(d)),paste_nv(sd,sd_pretty(sd))));
  save_(sim,file,save=save.sim,save.txt=save.txt.sim);
}
load_sim_hetd=function(file=NULL,n,d,sd) {
  if (is.null(file))
    file=filename(param(sim.hetd.dir),base='sim',
                  tail=paste(sep=',',paste_nv(n),
                             paste_nv(d,d_pretty(d)),paste_nv(sd,sd_pretty(sd))));
  load_(file,'sim');
}
get_sim_hetd=load_sim_hetd;

##### interp_hetd. interp is approx data, NOT function, 'cuz it's dangerour to save approxfun
save_interp_hetd=function(interp,n,d,sd,file=NULL) {
  param(save.interp,save.txt.interp);
  if (is.null(file))
    file=filename(param(interp.hetd.dir),base='interp',
                  tail=paste(sep=',',paste_nv(n),paste_nv(d,d_pretty(d)),sd_pretty(sd)));
  save_(interp,file,save=save.interp,save.txt=save.txt.interp);
}
load_interp_hetd=function(file=NULL,n,d,sd) {
  if (is.null(file))
    file=filename(param(interp.hetd.dir),base='interp',
                  tail=paste(sep=',',paste_nv(n),paste_nv(d,d_pretty(d)),sd_pretty(sd)));
  interp=load_(file,'interp');
  approxfun(interp)
}
get_interp_hetd=load_interp_hetd;

## TODO 19-05-09: delete commented out code when I'm sure I'm happy with the new way
## ##### meand_empi
## save_meand_empi=function(meand,file=NULL) {
##   param(save.meand,save.txt.meand);
##   if (is.null(file)) file=filename(param(datadir),base='meand.empi');
##   save_(meand,file,save=save.meand,save.txt=save.txt.meand);
## }
## load_meand_empi=function(file=NULL) {
##   if (is.null(file)) file=filename(param(datadir),base='meand.empi');
##   load_(file,'meand');
## }
## get_meand_empi=load_meand_empi;
## ##### meand_theo
## save_meand_theo=function(meand,file=NULL) {
##   param(save.meand,save.txt.meand);
##   if (is.null(file)) file=filename(param(datadir),base='meand.theo');
##   save_(meand,file,save=save.meand,save.txt=save.txt.meand);
## }
## load_meand_theo=function(file=NULL) {
##   if (is.null(file)) file=filename(param(datadir),base='meand.theo');
##   load_(file,'meand');
## }
## get_meand_theo=load_meand_theo;

##### table - saved in tbldir
save_tbl=function(tbl,file,obj.ok=F) {
  param(save.tbl,save.txt.tbl);
  if (is.null(tbl)) stop('Trying to save NULL table. Is table name set correctly?');
  base=desuffix(file);
  file=filename(base=base,suffix='RData');
  if ((is.na(save.tbl)&!file.exists(file))|(!is.na(save.tbl)&save.tbl)) {
    save(tbl,file=file);
    if (save.txt.tbl) {
      file=filename(base=base,suffix='txt');
      if (length(dim(tbl))==2) write.table(tbl,file=file,sep='\t',quote=F,row.names=F)
      else if (is.list(tbl)) {
        sink(file);
        print(tbl);
        sink();
      }
      else if (is.vector(tbl)) {
        names=names(tbl);
        if (!is.null(names)) {
          tbl=data.frame(name=names,value=as.character(tbl));
          write.table(tbl,file=file,sep='\t',quote=F,row.names=F);
        } else writeLines(as.character(tbl),file);
      }
      else if (!obj.ok) stop('Trying to save generic object but obj.ok=F.');
    }}
  invisible(tbl);
}

##### data - arbitrary objects saved in datadir
filename_data=function(what,suffix='RData')
  filename(basename_data(what),suffix=suffix);
basename_data=function(what) filename(param(datadir),base=what)

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

