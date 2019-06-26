#################################################################################
##
## Author:  Nat Goodman
## Created: 19-01-01
##          from repwr/R/util.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Utility functions for effit
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################

## ---- Utility Functions ----
## generate name=value
paste_nv=function(name,value,sep='=') {
  name=as.character(pryr::subs(name));
  if (missing(value))
    if (exists(name,envir=parent.frame(n=1))) value=get(name,envir=parent.frame(n=1))
    else stop(paste('no value for',name,'in function call or parent environment'));
  paste(sep=sep,name,value); 
}
## generate list of name=value using values from parent environment. code adapted from base::rm
## ignore tells whether to ignore NULL and non-existant names
nvq=function(...,sep=' ',ignore=F) {
  dots=match.call(expand.dots=FALSE)$...
   if (length(dots) &&
     !all(vapply(dots,function(x) is.symbol(x) || is.character(x),NA,USE.NAMES=FALSE))) 
     stop("... must contain names or character strings");
  ## CAUTION: for some reason, doesn't work to use 'parent.frame(n=1)' inside sapply
  env=parent.frame(n=1);
  names=vapply(dots,as.character,"");
  values=sapply(names,function(name) {
    if (exists(name,envir=env)) get(name,envir=env)
    else if (!ignore) stop(paste('no value for',name,'in parent environment'));
  });
  ## values=sapply(names,function(name)
  ##   if (exists(name,envir=parent.frame(n=2))) get(name,envir=parent.frame(n=2))
  ##   else stop(paste('no value for',name,'in parent environment')));
  paste(collapse=sep,unlist(mapply(function(name,value)
    if (!is.null(value)|!ignore) paste(sep='=',name,value) else NULL,names,values)));
}
## tack id onto filebase if not NULL or NA
paste_id=function(base,id=NULL,sep='.') {
  ## test id this way to avoid running is.na when id=NULL 
  if (is.null(id)) return(base);
  if (is.na(id)) return(base);
  paste(sep=sep,base,id);
}  
## pretty print typical values of n, d, sd & m
n_pretty=function(n) as.character(n);
d_pretty=function(d) sprintf('%3.2f',d);
sd_pretty=function(sd) sprintf('%3.2f',sd);
m_pretty=function(m) {
  if (round(log10(m))==log10(m)) sub('e\\+0{0,1}','e',sprintf("%0.0e",m),perl=TRUE)
  else as.character(m);
}
## NG 19-01-31: CAUTION. doesn’t really work
##   supposed to search dynamic environment tree but does static instead
##   seemed to work “back in the day” because params were global and static predecessor
##   of most functions is the global environment 
## get value of variable from parent or set to default
## call with quoted or unquoted variable name
## if default missing, throws error if variable not found
parent=function(what,default) {
  what=as.character(pryr::subs(what));
  if (exists(what,envir=parent.frame(n=2))) return(get(what,envir=parent.frame(n=2)));
  if (!missing(default)) return(default);
  stop(paste(sep='',"object '",what,"' not found"));
}
## get value of variable from param environment and assign to same named variable in parent
## call with quoted or unquoted variable names
## adapted from base::rm
## set params using par-like notation, eg, param(m=1e4)
param=function(...,list=character()) {
  dots=match.call(expand.dots=FALSE)$...
  parent.env=parent.frame(n=1);
  ## handle params with new values
  names=names(dots);
  if (!is.null(names)) {
    dots=sapply(seq_along(dots),function(i) {
      if (nchar(names[i])==0) return(dots[i]);
      ## set new value in param.env
      what=names[i];
      val=eval(dots[[i]],envir=parent.env);
      assign(what,val,envir=param.env);
      ## replace element in dots with name so it'll get returned
      what;
    })}
  if (length(dots) &&
      !all(vapply(dots,function(x) is.atomic(x)||is.symbol(x)||is.character(x),
                  NA,USE.NAMES=FALSE))) 
    stop("... must contain atomic data like names or character strings");
  names=vapply(dots,as.character,"");
  if (length(names)==0L) names=character();
  names=c(list,names);
  ## make sure all params valid
  bad=which(sapply(names,function(name) !exists(name,envir=param.env)));
  if (any(bad)) stop(paste(sep=' ','Invalid param(s):',paste(collapse=', ',names(bad))))
  retval=lapply(names,function(name) assign(name,get(name,envir=param.env),envir=parent.env));
  ## fix up return value
  if (length(retval)==1) unlist(retval)
  else {
    names(retval)=names;
    retval;
  }
}
## copy local variables to global - to simplify init
## NG 19-01-11: used in repwr, not in effit
assign_global=function() {
  env=parent.frame(n=1);
  sapply(ls(envir=env),function(what) assign(what,get(what,envir=env),envir=.GlobalEnv));
}
## copy local variables to new or existing param environment - to simplify init
init_param=function() {
  param.env=new.env(parent=emptyenv());
  parent.env=parent.frame(n=1);
  sapply(ls(envir=parent.env),
        function(what) assign(what,get(what,envir=parent.env),envir=param.env));
  assign('param.env',param.env,envir=.GlobalEnv);
}
assign_param=function() {
  parent.env=parent.frame(n=1);
  sapply(ls(envir=parent.env),
        function(what) assign(what,get(what,envir=parent.env),envir=param.env));
}
## copy variable to parent
## NG 19-01-11: not used in effit. used once upon a time in dofig to update fignum
assign_parent=function(what,value) {
  what=as.character(pryr::subs(what));
  if (missing(value)) value=get(what,envir=parent.frame(n=1));
  assign(what,value,envir=parent.frame(n=2));
}
## NG 18-10-24: wrap function - propogate locals and ... then call function
##   funfun are additional functions called by fun with ... args
## TODO: handle partial matching of ... params
## adapted from stackoverflow.com/questions/4124900
wrap_fun=function(fun,funfun=NULL,...) {
  env=parent.frame(n=1);
  x=ls(envir=env);
  fx=do.call(c,lapply(c(fun,funfun),function(fun) names(formals(fun))));
  args=sapply(x[x%in%fx],function(x) get(x,envir=env),simplify=F);
  dots=list(...);
  args=c(args,dots[names(dots)%in%fx]);
  do.call(fun,args);
}

## like match.arg but uses general matching and, if several.ok, returns 'em all
pmatch_choice=
  function(arg,choices,several.ok=T,none.ok=F,start=T,ignore.case=T,perl=F,fixed=F,invert=F) {
    ## m=startsWith(choices,arg);
    pat=if(start) paste0('^',arg) else arg;
    m=grep(pat,choices,ignore.case=ignore.case,perl=perl,value=T,fixed=fixed,invert=invert);
    if (length(m)==0&&!none.ok)
      stop(paste(sep=' ',"'arg' matched none of",paste(collapse=', ',choices),
           "but 'none.ok' is FALSE"));
    if (length(m)>1&&!several.ok)
      stop(paste(sep=' ',"'arg' matched several of",paste(collapse=', ',choices),
                 "but 'several.ok' is FALSE"));
    if (length(m)==0) NULL else m;
  }
## quote names in paramter list. code adapted from base::rm
cq=function(...) {
 dots=match.call(expand.dots=FALSE)$...
 if (length(dots) &&
     !all(vapply(dots,function(x) is.atomic(x)||is.symbol(x)||is.character(x),
                 NA,USE.NAMES=FALSE))) 
   stop("... must contain atomic data like names or character strings");
 return(vapply(dots,as.character,""));
}
## upper case first letter of word. like Perl's ucfirst
## from https://stackoverflow.com/questions/18509527/first-letter-to-upper-case/18509816
ucfirst=function(word) paste0(toupper(substr(word,1,1)),substr(word,2,nchar(word)));

## extend akima::aspline for matrix
asplinem=function(x,y,xout,...) {
  if (is.vector(y)) y=data.frame(y=y);
  if (length(dim(y))!=2) stop('y must be vector or 2-dimensional matrix-like object');
  ## yout=apply(y,2,function(y) akima::aspline(x,y,xout,...)$y);
  ## extend y to correct number of rows if necessary. not really necessary for aspline
  ## CAUTION: perhaps this should be error...
  if (nrow(y)<length(x)) y=repr(y,length=length(x));
  yout=apply(y,2,function(y) {
    if (all(is.na(y))) rep(NA,length(xout))
    else if (length(which(!is.na(y)))==1) rep(y[which(!is.na(y))],length(xout))
    else akima::aspline(x,y,xout,...)$y;})
  ## if yout has single row (ie, xout has one element), R turns it into a vector...
  if (length(xout)==1) yout=t(yout);
  yout;
}
## extend loess.smooth for matrix - probably only useful for plotting
loessm=function(x,y,xout,...) {
  if (is.vector(y)) y=data.frame(y=y);
  if (length(dim(y))!=2) stop('y must be vector or 2-dimensional matrix-like object');
  ## extend y to correct number of rows if necessary
  ## CAUTION: perhaps this should be error...
  if (nrow(y)<length(x)) y=repr(y,length=length(x));
  data=data.frame(x=x,y);
  yout=do.call(data.frame,lapply(colnames(data)[2:ncol(data)],function(name) {
    fmla=as.formula(paste(name,'~ x'));
    yout=predict(loess(fmla,data=data),xout)
  }));
  ## if yout has single row (ie, xout has one element), R turns it into a vector...
  ## if (length(xout)==1) yout=t(yout);
  colnames(yout)=colnames(y);
  yout;
}
## extend smooth.spline for matrix - probably only useful for plotting
## NG 18-11-07: remove NAs (same as akima::aspline) else smooth.spline barfs
splinem=function(x,y,xout,spar=NULL,...) {
  if (is.null(spar)) spar=0.5;
  if (is.vector(y)) y=data.frame(y=y);
  if (length(dim(y))!=2) stop('y must be vector or 2-dimensional matrix-like object');
  ## extend y to correct number of rows if necessary
  ## CAUTION: perhaps this should be error...
  if (nrow(y)<length(x)) y=repr(y,length=length(x));
  yout=apply(y,2,function(y) {
    ## remove NAs. code adapted from akima::aspline
    ## CAUTION: must use '<-' not '=' or place assignment in extra parens ((na=is.na(y)))
    ##   see stackoverflow.com/questions/1741820 for explanation. gotta love R...
    if (any(na<-is.na(y))) {
      x=x[!na]; y=y[!na];
    }
    yout=predict(smooth.spline(x,y,spar=spar),xout)$y    
  });
  ## if yout has single row (ie, xout has one element), R turns it into a vector...
  ## if (length(xout)==1) yout=t(yout);
  if (length(xout)==1) yout=t(yout);
  colnames(yout)=colnames(y);
  yout;
}
## with case
## like 'with' but works on vectors. I use it inside apply(cases,1,function(case)...)
## note that plain 'with' works fine when applied to cases as a whole
## withcase=function(case,...) with(data.frame(t(case)),...)
withcase=function(case,...) {
  case=data.frame(t(case),stringsAsFactors=FALSE);
  assign('case',case,envir=parent.frame(n=1)); # so case will be data frame in called code
  with(case,...);
}  

## round up or down to nearest multiple of u. from https://grokbase.com/t/r/r-help/125c2v4e14/
round_up=function(x,u) ceiling(x/u)*u;
round_dn=function(x,u) floor(x/u)*u;
## x can be range or single number (lower bound)
round_rng=function(x,y,u) 
  if (missing(y)) c(round_dn(x[1],u),round_up(x[2],u)) else c(round_dn(x,u),round_up(y,u))

## pick n items from x approx evenly spaced
pick=function(x,n.want,n.min=1,rep.ok=FALSE,exclude=NULL) {
  x=x[x%notin%exclude];
  if (length(x)<n.min) stop('too few elements in x');
  if (length(x)<n.want&!rep.ok) x
  else {
    step=1/(n.want+1);
    probs=seq(step,by=step,len=n.want)
    unname(quantile(x,probs=probs,type=1))
  };
}

## repeat rows or columns of 2-dimensional matrix-like object. like rep
## like rep, ... can be times, length.out, or each
## based on StackOverflow https://stackoverflow.com/questions/11121385/repeat-rows-of-a-data-frame
repr=function(x,...) {
  i=rep(seq_len(nrow(x)),...);
  x=x[i,,drop=F];
  rownames(x)=NULL;
  x;
}
repc=function(x,...) {
  j=rep(seq_len(ncol(x)),...);
  x=x[,j,drop=F];
  colnames(x)=NULL;
  x;
}
## not in - based on example in RefMan - more intutive than !%in%
"%notin%"=function(x,table) match(x,table,nomatch=0)==0
## between, near - to subset sim results. closed on bottom, open on top
between=function(x,lo,hi) x>=lo&x<hi
near=function(x,target,tol=.01) between(x,target-tol,target+tol)

## debugging functions
## TODO: BREAKPOINT is sooo feeble :(
BREAKPOINT=browser;
## traceback with args I like
tback=function(max.lines=2) traceback(max.lines=max.lines)
devs.close=function() for (dev in dev.list()) dev.off(dev)
## display color palette
pal=function(col,border="light gray",...) {
 n=length(col)
 plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="",...)
 rect(0:(n-1)/n,0,1:n/n,1,col=col,border=border)
}
