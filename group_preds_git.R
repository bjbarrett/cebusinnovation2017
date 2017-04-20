R
library(rethinking)
library(lubridate)
library(RColorBrewer)

c <- read.csv("PNASfinaldomainindividual_git.csv" , header=TRUE)
length(unique(c$id[c$innovator==1]))
length(unique(c$id[c$inno_count==1]))
length(unique(c$id[c$inno_count==2]))
length(unique(c$id[c$inno_count==3]))
length(unique(c$id[c$inno_count==4]))
length(unique(c$id[c$inno_count==5]))
length(unique(c$id[c$inno_count==0]))

max(c$inno_count)

df <- subset(c, select=c("year" , "sampyearsum" , "domain" , "domain_index" , "MainGroup"  , "group_index" , "year_index" , "inno_count", "inno_count_ND"))

df$GDYsum <- df$GYsum <- df$sampeff <- -1

for (i in 1:max(df$year_index)){
	for (j in 1:max(df$group_index)){
		for (k in 1:max(df$domain_index)){
			df$GDYsum <- ifelse(df$year_index==i & df$group_index==j & df$domain_index==k, sum(df$inno_count[df$year_index==i & df$group_index==j & df$domain_index==k]) , df$GDYsum)
		}
		df$GYsum <-ifelse( df$year_index==i & df$group_index==j , sum(df$inno_count[df$year_index==i & df$group_index==j ]) , df$GYsum )
		df$sampeff <- ifelse( df$year_index==i & df$group_index==j , sum(df$sampyearsum[df$year_index==i & df$group_index==j ]) , df$sampeff )
	}
}


d <- subset(df , select= -c(inno_count,inno_count_ND,sampyearsum) )

d <- d[!duplicated(d), ]

dd <- subset(d , select=c(GYsum,sampeff,year_index,group_index,MainGroup,year ) )
dd <- dd[!duplicated(dd), ]

##add 1 to SP 2010 for XB who was excluded from main analysis due to uncertainty about ID

dd$GYsum[dd$MainGroup=="SP" & dd$year==2010] <- dd$GYsum[dd$MainGroup=="SP" & dd$year==2010] +1

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



pdf("Fig1.groupannualinnoestnoyaxis.pdf", width= 11.4/2.54, height=5.25/2.54 )

par(mfrow = c(2,5))
par(cex = .6)
par(mar = c(0, 0, 0, 0), oma = c(2.75, 0.25, 0.25, 0.25) ,  mgp=c(3, .25, 0))

for(i in 1){
d.pred <- list(
	domain_index=1,
	group_index=i,
	logsampeff=mean(log(dd$sampeff/365))
	)

link2 <- link(m, n=1000 , data=d.pred)
pred.med <- apply(link2  , 2 , median)
dens(link2, col="white", main="" , ylim=c(0,0.9) , xlim=c(0,10) , xlab="" ,ylab="" , yaxt="n", xaxt="n" ,zero.line = FALSE , cex.axis=0.75)
shade( density(link2) , lim= as.vector(HPDI(link2, prob=0.99999)) , col = col.alpha("orange", 0.75))

abline(v=pred.med)
 title(grouplist[i], line = -1 , cex=.5)
}

 axis(1, seq(from=0 , to=10 , by = 2), labels=FALSE , tck=0.02)
 axis(1, seq(from=0 , to=10 , by = 1), labels=FALSE , tck=0.01)

for(i in 2:5){
d.pred <- list(
	domain_index=1,
	group_index=i,
	logsampeff=mean(log(dd$sampeff/365))
	)

link2 <- link(m, n=1000 , data=d.pred)
pred.med <- apply(link2  , 2 , median)
dens(link2, col="white", main="" , ylim=c(0,0.9) , xlim=c(0,10) , xlab="" ,ylab="" , yaxt="n", xaxt="n" ,zero.line = FALSE)
shade( density(link2) , lim= as.vector(HPDI(link2, prob=0.99999)) , col = col.alpha("orange", 0.75))
abline(v=pred.med)
 title(grouplist[i], line = -1 , cex=.5)
 axis(1, seq(from=0 , to=10 , by = 2), labels=FALSE , tck=0.02)
 axis(1, seq(from=0 , to=10 , by = 1), labels=FALSE , tck=0.01)
}
 axis(1, seq(from=0 , to=10 , by = 2), labels=FALSE , tck=0.02)
 axis(1, seq(from=0 , to=10 , by = 1), labels=FALSE , tck=0.01)
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
#dens(link2, col="grey", main="" ,  show.HPDI=0.99999 , ylim=c(0,1) , xlim=c(0,12) , xlab="" ,ylab="")
dens(link2, col="white", main="" , ylim=c(0,0.9) , xlim=c(0,10) , xlab="" ,ylab="", xaxt="n" , yaxt="n", zero.line = FALSE , cex.axis=0.75)
shade( density(link2) , lim= as.vector(HPDI(link2, prob=0.99999)) , col = col.alpha("orange", 0.75))
axis(1, seq(from=0 , to=10 , by = 1), labels=FALSE , tck=-0.03)
axis(1, seq(from=0 , to=10 , by = 2), labels=c("",2,4,6,8,10) , tck=-0.06 ,  cex.axis=1)

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
dens(link2, col="white", main="" , ylim=c(0,0.9) , xlim=c(0,10) , xlab="" ,ylab="" , yaxt="n",xaxt="n", zero.line = FALSE)
shade( density(link2) , lim= as.vector(HPDI(link2, prob=0.99999)) , col = col.alpha("orange", 0.75))

abline(v=pred.med)
 title(grouplist[i], line = -1 , cex=.5 )
 axis(1, seq(from=0 , to=10 , by = 1), labels=FALSE , tck=-0.03, cex.axis=1 , cex.lab=1.5)
 axis(1, seq(from=0 , to=10 , by = 2), labels=c("",2,4,6,8,10) , tck=-0.06 , cex.axis=1, cex.lab=1.5)

}

mtext("estimated annual innovation rate per group", outer =TRUE, cex = 1, side=1, line=1.5)
#mtext("posterior density", outer =TRUE, cex = .8, side=2, line=2)


dev.off()


