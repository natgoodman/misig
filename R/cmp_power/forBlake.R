## Do it!
doit=
  function(n=c(21,51,310),mu=c(0.8,0.5,0.2),tau2=seq(0,0.125,by=0.025),m=1e5,
           gen.cases='default') {
    power=cmp_power(n=n,mu=mu,tau2=tau2,m=m,gen.cases=gen.cases);
    plot_power(power);
    return(power);
  }
## Power function from Blake
power.singleeffect.t.calc <- function(mu, Sigma, Sigma.null, df=Inf, alpha=0.05, alternative=NULL){				
	# alternative = "two.sided" or "one.sided"
	if(is.null(alternative)){alternative <- "one.sided"}
	tcut <- ifelse(alternative=="two.sided", qt(1-alpha/2,df=df), qt(1-alpha,df=df))
	s <- sqrt(diag(Sigma))
	s.null <- sqrt(diag(Sigma.null))
	if(alternative=="one.sided"){
		tmp <- (tcut*s.null - abs(mu)) / s
		power <- 1 - pt( tmp, df=df)
	}
		
	if(alternative=="two.sided"){
		tmp1 <- (tcut*s.null - abs(mu)) / s
		tmp2 <- (-tcut*s.null - abs(mu)) / s
		power <- 1 - pt( tmp1, df=df ) +
				     pt( tmp2, df=df )
	}
	
	list(power=power)
}
####################
## functions used to compare power and plot results
power_blake=function(n,mu,tau2,alpha=0.05,alternative=c('one.sided','two.sided')) {
  alternative=match.arg(alternative);
  X0 <- cbind(1, c(0,1))
  C <- nrow(X0)
  sigma2 <- 1
  Sigma.a <- matrix(0, C, C)
  Sigma.b <- tau2 * solve(t(X0)%*%X0)
  mu0 <- c(0, mu);
  Sigma.ebar <- sigma2 / n * solve(t(X0)%*%X0)
  S.sampling <- Sigma.a + Sigma.b + Sigma.ebar
  S.null <- Sigma.ebar
  effect=2;
  power.singleeffect.t.calc(mu0, S.sampling, S.null, df=(n-1)*C, alpha=alpha, alternative=alternative)$power[effect]  
}
vrnorm=Vectorize(rnorm,"mean");
pooled_sd=function(sd0,sd1) sqrt((sd0^2+sd1^2)/2);
power_sim=
  function(n,mu,sd.het,tau2,alpha=0.05,alternative=c('one.sided','two.sided'), strict=F,
           method=cq(blake,sem,nat),m=1e3) {
    ## Blake suggests this transformation instead
    ## if (missing(tau2)) tau2=sqrt(2)*sd.het^2 else sd.het=sqrt(tau2/sqrt(2));
    method=match.arg(method);
    if (missing(tau2)) tau2=sd2tau2(sd.het=sd.het,n=n,method=method)
    else sd.het=tau2sd(tau2=tau2,n=n,method=method);
    alternative=match.arg(alternative);
    if (alternative=='one.sided') alpha=alpha*2;
    ## print(paste(sep='','>>> power_sim: ',paste(sep=', ',n,mu,round(sd.het,digits=3),tau2)));
    d.pop=rnorm(m,mean=mu,sd=sd.het);
    group0=replicate(m,rnorm(n,mean=0));
    group1=vrnorm(n,mean=d.pop);
    mean0=colMeans(group0);
    mean1=colMeans(group1);
    d.raw=mean1-mean0;
    sd0=apply(group0,2,sd);
    sd1=apply(group1,2,sd);
    sd=pooled_sd(sd0,sd1);
    d.sdz=d.raw/sd;
    d.crit=d_crit(n,sig.level=alpha);
    d.htcrit=d_htcrit(n,sd.het,sig.level=alpha);
    if (strict) d.sdz=abs(d.sdz);
    power.sim=length(which(d.sdz>=d.htcrit))/m;
    power.simconv=length(which(d.sdz>=d.crit))/m;
    data.frame(n,mu,tau2,power.sim,power.simconv);
}
cmp_power=
  function(n,mu,sd.het,tau2,alpha=0.05,m=1e3,alternative=c('one.sided','two.sided'),strict=F,
           gen.cases=c('default','expand.grid','data.frame'),method=cq(blake,sem,nat)) {
    ## if (missing(tau2)) tau2=sd.het^2 else sd.het=sqrt(tau2);
    ## Blake suggests this transformation instead
    ## if (missing(tau2)) tau2=sqrt(2)*sd.het^2 else sd.het=sqrt(tau2/sqrt(2));
    method=match.arg(method);
    if (missing(tau2)) tau2=sd2tau2(sd.het=sd.het,n=n,method=method)
    else sd.het=tau2sd(tau2=tau2,n=n,method=method);
    alternative=match.arg(alternative);
    if (alternative=='one.sided') alpha=alpha*2;
    gen.cases=match.arg(gen.cases);
    cases=gen_cases(n,mu,tau2,sd.het,gen.cases);
    do.call(rbind,apply(cases,1,function(case) {
      n=case['n']; mu=case['mu']; tau2=case['tau2']; sd.het=case['sd.het'];
      power.blake=power_blake(n,mu,tau2);
      power.nat=power_d2ht(n=n,d=mu,sd.het=sd.het,sig.level=alpha,strict=strict);
      power.natconv=power_d2ht(n=n,d=mu,sd.het=sd.het,d.crit=d_crit(n,sig.level=alpha),
                               strict=strict);
      power.sim=power_sim(n,mu,tau2=tau2,sd.het=sd.het,alternative=alternative,strict=strict,
                          method=method,m=m);
      data.frame(n,mu,tau2,power.blake,power.nat,power.sim=power.sim$power.sim,
                 power.natconv,power.simconv=power.sim$power.simconv)
    }));
  }
gen_cases=function(n,mu,tau2,sd.het,gen.cases) {
  sd.het=setNames(sd.het,as.character(tau2));
  cases=switch(gen.cases,
               default={
                 ntau2=expand.grid(n=n,tau2=tau2);
                 mutau2=expand.grid(mu=mu,tau2=tau2);
                 data.frame(n=ntau2$n,mu=mutau2$mu,tau2=ntau2$tau2)},
               expand.grid=expand.grid(n=n,mu=mu,tau2=tau2),
               data.frame=data.frame(n=n,mu=mu,tau2=tau2));
  cases$sd.het=with(cases,sd.het[as.character(tau2)]);
  cases[with(cases,order(n,mu,tau2)),];
}
plot_power=
  function(power,title='power: nat vs blake',
           pch=c(20,1,20,1),cex=c(0.75,1,0.75,1),col=c('red','blue','salmon','cyan'),
           xlab='blake',ylab='nat',xlim=c(0,1),ylim=c(0,1),
           legend.where='topleft',legend.title='nat power',
           legend=c('analytic w/ correct pvals','simulated w/ correct pvals',
                    'analytic w/ conventional pvals','simulated w/ conventional pvals'),...) {
    dev.new();
    blake=power$power.blake;
    nat=subset(power,select=c(power.nat,power.sim,power.natconv,power.simconv));
    matplot(blake,nat,pch=pch,cex=cex,col=col,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=title,
            ...);
    abline(a=0,b=1,col='grey',lty='dotted');
    grid();
    legend(legend.where,title=legend.title,legend=legend,col=col,pch=pch,cex=0.8,bty='n');
  }
#####
## functions to convert between tau and sd.het
## 1) my original: tau=sd.het
## 2) Blake's suggestion: tau2=sqrt(2)*sd.het^2; sd.het=sqrt(tau2/sqrt(2));
## 3) tau=standard error of mean=sd.het/sqrt(n); sd.het=tau*sqrt(n)
tau2sd=function(tau=NA,tau2=NA,n,method=cq(blake,sem,nat)) {
  method=match.arg(method);
  if (missing(tau)&missing(tau2)) stop('Either tau or tau2 must be set');
  if (missing(tau)) tau=sqrt(tau2) else tau2=tau^2;
  switch(method,sem=tau*sqrt(n),blake=sqrt(tau2/sqrt(2)),nat=tau);
}
sd2tau=function(sd.het,n,method=cq(blake,sem,nat)) {
  method=match.arg(method);
  switch(method,sem=sd.het/sqrt(n),blake=sqrt(sqrt(2)*sd.het^2),nat=sd.het);
}
sd2tau2=function(sd.het,n,method=cq(blake,sem,nat)) {
  method=match.arg(method);
  switch(method,sem=sd.het^2/n,blake=sqrt(2)*sd.het^2,nat=sd.het^2);
}

####################
## functions extracted from R/stats.R.
## wrapper for R's t distribution specialized for 
##   two group difference-of-mean. Cohen's d.
##   same sample size each group. normal populations w/ unit variance
## handles conventional and het p-values
##
## functions for conventional p-values
## significance boundary for Cohen's d
d_crit=d_sig=function(n,sig.level=0.05) pval2d(n,pval=sig.level)

## pval to Cohen's d
pval2d=function(n,pval) q_d2t(n,p=pval/2,lower.tail=F)

## noncentrality parameter
ncp=function(n,d) sqrt(n/2)*d

## t statistic to/from Cohen's d & pval
t2d=function(n,t) t*sqrt(2/n);
t2pval=function(n,t) 2*pt(-abs(t),df=2*(n-1))
d2t=function(n,d) d*sqrt(n/2);
d2pval=function(n,d) t2pval(n,d2t(n,d))
pval2t=function(n,pval) qt(pval/2,df=2*(n-1),lower.tail=F)
pval2d=function(n,pval) q_d2t(n,p=pval/2,lower.tail=F)

## quantile function
q_d2t=function(n,p,d0=NULL,lower.tail=TRUE) {
  df=2*(n-1);
  if (!is.null(d0)) t=suppressWarnings(qt(p,df=df,ncp=ncp(n,d0),lower.tail=lower.tail))
    else t=qt(p,df=df,lower.tail=lower.tail)
  t2d(n,t);
}
## cumulative distribution function
d_d2t=function(n,d,d0=NULL) {
  df=2*(n-1);
  t=d2t(n,d);
  sqrt(n/2)*
    if (!is.null(d0)) suppressWarnings(dt(t,df=df,ncp=ncp(n,d0)))
    else dt(t,df=df)
}
## functions for het p-values
## Kenny & Judd scale factor
q_=function(n,sd.het) sqrt((2/n)/((2/n)+sd.het^2));

## significance boundary for Cohen's d
d_htcrit=function(n,sd.het,sig.level=0.05) htpval2d(n,sd.het,pval=sig.level)

## pval to Cohen's d
htpval2d=function(n,sd.het,pval) q_d2ht(n,sd.het=sd.het,p=pval/2,lower.tail=F)

## noncentrality parameter
ncp_ht=function(n,d.het,sd.het) sqrt(n/(2+(n*sd.het^2)))*d.het;

## het t statistic to/from Cohen's d, conventional t, pval
t2ht=function(n,sd.het,t) t*q_(n,sd.het);
ht2t=function(n,sd.het,ht) ht/q_(n,sd.het);
d2ht=function(n,sd.het,d) t2ht(n,sd.het,d2t(n,d));
ht2d=function(n,sd.het,ht) t2d(n,ht2t(n,sd.het,ht));
d2htpval=function(n,sd.het,d) t2pval(n,d2ht(n,sd.het,d))
htpval2d=function(n,sd.het,pval) q_d2ht(n,sd.het=sd.het,p=pval/2,lower.tail=F)

## quantile function
q_d2ht=function(n,d.het=NULL,sd.het,p,lower.tail=TRUE) {
  df=2*(n-1);
  if (!is.null(d.het))
    ht=suppressWarnings(qt(p,df=df,ncp=ncp_ht(n,d.het,sd.het),lower.tail=lower.tail))
  else ht=qt(p,df=df,lower.tail=lower.tail)
  ht2d(n,sd.het,ht);
}
## cumulative distribution function
p_d2ht=function(n,d.het=NULL,sd.het,d,lower.tail=TRUE) {
  df=2*(n-1);
  ht=d2ht(n,sd.het,d);
  if (!is.null(d.het))
    suppressWarnings(pt(ht,df=df,ncp=ncp_ht(n,d.het,sd.het),lower.tail=lower.tail))
  else pt(ht,df=df,lower.tail=lower.tail)
}
## power. strict=F only considers positive tail. default for power.t.test
power_d2ht=function(n,d,sd.het,d.crit,sig.level=0.05,strict=F) {
  if (missing(d.crit)) d.crit=d_htcrit(n=n,sd.het=sd.het,sig.level=sig.level);
  p_d2ht(n=n,d.het=d,sd.het=sd.het,d=d.crit,lower.tail=F)+
    ifelse(strict,p_d2ht(n=n,d.het=d,sd.het=sd.het,d=-d.crit,lower.tail=T),0)
}
