# Power function
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


# Gpower results: sample sizes are per cell
n.gpower <- c(21, 51, 310)		# One-sided
power.gpower <- c(0.8167878, 0.8058986, 0.8002178)

# Model parameters:
X0 <- cbind(1, c(0,1))
C <- nrow(X0)
mu.vec <- c(0.8, 0.5, 0.2)
sigma2 <- 1
tau2.vec <- seq(0, 3*0.20^2, length=201) #seq(0,2,length=21) * sigma2/n.orig.per.cell
alpha <- 0.05

# Compute sample size required for 80% power
# Total sample size required is optimized value * C
tmpf <- function(n.per.cell, target=0.80, effect=1){
	Sigma.ebar <- sigma2 / n.per.cell * solve(t(X0)%*%X0)
	S.sampling <- Sigma.a + Sigma.b + Sigma.ebar
	S.null <- Sigma.ebar
	tmp.power <- power.singleeffect.t.calc(mu, S.sampling, S.null, df=(n.per.cell-1)*C, alpha=alpha, alternative="one.sided")$power[effect]
	(tmp.power - target)^2
}

pp <- array(NA, dim=c(length(tau2.vec), length(mu.vec)))
nn <- array(NA, dim=c(length(tau2.vec), length(mu.vec)))
for(j in 1:length(tau2.vec)){
	Sigma.a <- matrix(0, C, C)
	Sigma.b <- tau2.vec[j] * solve(t(X0)%*%X0)
	for(k in 1:length(mu.vec)){
		mu <- c(0, mu.vec[k])
		pp[j,k] <- sqrt(tmpf(n.gpower[k], target=0, effect=2))
		
		tmpnM <- 10
		tmpn <- optimize(tmpf, interval=c(0,tmpnM), effect=2)$minimum
		while(ceiling(tmpn)==tmpnM & tmpnM<1000){
			tmpnM <- 10*tmpnM
			tmpn <- optimize(tmpf, interval=c(0,tmpnM), effect=2)$minimum
		}
		nn[j,k] <- ifelse(ceiling(tmpn)==tmpnM, NA, tmpn)
	}
}

# Make plots
library(ggplot2)
library(reshape2)
tmpgg <- melt(nn)
tmpgg$nceiling <- ceiling(tmpgg$value)
tmpgg$powerg <- melt(pp)$value
tmpgg$tau2 <- tau2.vec[tmpgg$Var1]
tmpgg$Effect <- c("Large (d=0.8)", "Medium (d=0.5)", "Small (d=0.2)")[tmpgg$Var2]
tmpgg$n.gpower <- n.gpower[tmpgg$Var2]
tmpgg$power.gpower <- power.gpower[tmpgg$Var2]
tmpgg$mu <- paste0("mu==",mu.vec[tmpgg$Var2])
tmpgg$ratio <- tmpgg$nceiling/tmpgg$n.gpower
tmpgg$rawmu <- mu.vec[tmpgg$Var2]


pdf("graph_twocell_n.pdf", height=9.1, width=3.9)
ggplot(tmpgg, aes(x=tau2,y=nceiling)) + geom_line() +
	facet_grid(mu~., labeller=label_parsed, scales="free") +
	geom_hline(aes(yintercept=n.gpower), linetype=2, size=0.25) +
	xlab(expression(tau^2)) + ylab("n") 
dev.off()

pdf("graph_twocell_power_v2.pdf", height=3.5, width=4.5)
ggplot(tmpgg, aes(x=tau2,y=powerg,linetype=factor(rawmu))) + geom_line() +
	scale_linetype_discrete(expression(mu)) +
	xlab(expression(tau^2)) + ylab("Power") 
dev.off()

