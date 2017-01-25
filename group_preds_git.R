R
library(rethinking)
library(lubridate)
library(RColorBrewer)

dd <- read.csv("~/group_innovation_data_git.csv" , header=TRUE)

m <- map2stan(
    alist(
    inno ~ dpois( lambda ),
    log(lambda) <- a + ag  + logsampeff,
    #c(ad,ba_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    ag[group_index] ~ dnorm( 0 , sigma_group ),
	a ~ dnorm(0,3),
    sigma_group ~ dcauchy(0,2) 
),
data=list(
	inno = dd$GYsum,
	#domain_index = i$domain_index,
	group_index = dd$group_index,
	logsampeff=log(dd$sampeff/365)
	) , 
warmup=2000 , iter=4000 , chains=3 , cores=3 
)

col.pal <- brewer.pal(4, "Dark2")

dom_zeros <- matrix(0,1000,length(unique(d$domain_index)))
group_zeros <- matrix(0,1000,length(unique(d$group_index)))
gd_zeros <- matrix(0,1000,length(unique(d$gd_index)))
grouplist <- sort(unique(d$MainGroup))

pdf("Fig1.groupannualinnoest.pdf", width=8.5 , height=6 )
par(mfrow = c(2,5))
par(cex = .6)
par(mar = c(0, 0, 0, 0), oma = c(4.5, 4.5, 2,5, 1) )

for(i in 1){
d.pred <- list(
	domain_index=1,
	group_index=i,
	logsampeff=mean(log(dd$sampeff/365))
	)

link2 <- link(m, n=1000 , data=d.pred)
pred.med <- apply(link2  , 2 , median)
dens(link2, col="grey", main="" ,  show.HPDI=0.99999 , ylim=c(0,1) , xlim=c(0,12) , xlab="" ,ylab="" ,xaxt="n")
abline(v=pred.med)
 title(grouplist[i], line = -1 , cex=.5)
}

 axis(3, seq(from=0 , to=12 , by = 2), labels=seq(from=0 , to=12 , by = 2) , tck=-0.02)

for(i in 2:5){
d.pred <- list(
	domain_index=1,
	group_index=i,
	logsampeff=mean(log(dd$sampeff/365))
	)

link2 <- link(m, n=1000 , data=d.pred)
pred.med <- apply(link2  , 2 , median)
dens(link2, col="grey", main="" ,  show.HPDI=0.99999 , ylim=c(0,1) , xlim=c(0,12) , xlab="" ,ylab="" ,xaxt="n",yaxt="n")
abline(v=pred.med)
 title(grouplist[i], line = -1 , cex=.5)
  axis(3, seq(from=0 , to=12 , by = 2), labels=seq(from=0 , to=12 , by = 2) , tck=-0.02)

}
 axis(3, seq(from=0 , to=12 , by = 2), labels=seq(from=0 , to=12 , by = 2) , tck=-0.02)

#par(mar = c(1, 0, 0, 0), oma = c(3, 3, 3, 1) )
mfrow=c(1,11)
for(i in 6){
d.pred <- list(
	domain_index=1,
	group_index=i,
	logsampeff=mean(log(dd$sampeff/365))
	)

link2 <- link(m, n=1000 , data=d.pred)
pred.med <- apply(link2  , 2 , median)
dens(link2, col="grey", main="" ,  show.HPDI=0.99999 , ylim=c(0,1) , xlim=c(0,12) , xlab="" ,ylab="")
abline(v=pred.med)
 title(grouplist[i], line = -1 , cex=.5)
}

mfrow=c(1,11)
for(i in 7:10){
d.pred <- list(
	domain_index=1,
	group_index=i,
	logsampeff=mean(log(dd$sampeff/365))
	)

link2 <- link(m, n=1000 , data=d.pred)
pred.med <- apply(link2  , 2 , median)
dens(link2, col="grey", main="" ,  show.HPDI=0.99999 , ylim=c(0,1) , xlim=c(0,12) , xlab="" ,ylab="" , yaxt="n")
abline(v=pred.med)
 title(grouplist[i], line = -1 , cex=.5)
}

mtext("estimated annual innovation rate per group", outer =TRUE, cex = 1.4, side=1, line=3)
mtext("posterior density", outer =TRUE, cex = 1.4, side=2, line=2.5)


dev.off()








