R
library(rethinking)
library(lubridate)
library(RColorBrewer)


i <- read.csv("~/PNASinnovationlistused_git.csv" , header=TRUE)


 m1 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + bg*groupsize ,

    a_dom[domain_index] ~ dnorm( 0 , sigma_dom ),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
	c(a,bg) ~ dnorm(0,3),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	groupsize=i$groupsize.c

 

	) , 
warmup=1000 , iter=2000 , chains=3 , cores=3 ,types=list(adapt.delta=0.95)

)

m2 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + (ba + ba_dom)*age + bg*groupsize,

    c(a_dom,ba_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,ba,bg) ~ dnorm(0,3),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	#behavior_index = i$behavior_index,
	age = i$age.c,
	groupsize=i$groupsize.c


	) , 
warmup=1000 , iter=2000 , chains=3 , cores=3 ,types=list(adapt.delta=0.95)

)

m3 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group  + (bs + bs_dom)*soc + bg*groupsize,

    c(a_dom,bs_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,bs,bg) ~ dnorm(0,3),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	#behavior_index = i$behavior_index,
	soc= i$sociality.s,
	groupsize=i$groupsize.c


	) , 
warmup=1000 , iter=2000 , chains=3 , cores=3 ,types=list(adapt.delta=0.95)

)

m4 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group  + (bm + bm_dom)*male + bg*groupsize,

    c(a_dom,bm_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),
	c(a,bm,bg) ~ dnorm(0,3),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	#behavior_index = i$behavior_index,
	male= i$male,
	groupsize=i$groupsize.c

	) , 
warmup=1000 , iter=2000 , chains=3 , cores=3 ,types=list(adapt.delta=0.95)

)

m5 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + (brhi +brhi_dom)*rankhi + (brlo + brlo_dom)*ranklo + bg*groupsize,

    c(a_dom,brhi_dom,brlo_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,brhi,brlo,bg) ~ dnorm(0,3),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	#behavior_index = i$behavior_index,
	rankhi= i$rankhi,
	ranklo= i$ranklo,
	groupsize=i$groupsize.c

	) , 
warmup=1000 , iter=2000 , chains=3 , cores=3 ,types=list(adapt.delta=0.95)

)

m6 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + (ba + ba_dom)*age +  (bs + bs_dom)*soc +
    (bm + bm_dom)*male + (brhi +brhi_dom)*rankhi + (brlo + brlo_dom)*ranklo + bg*groupsize ,
    c(a_dom,ba_dom,bs_dom,bm_dom,brhi_dom,brlo_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,ba,bs,bm,brhi,brlo,bg) ~ dnorm(0,3),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	male= i$male,
	age=i$age.c,
	#behavior_index = i$behavior_index,
	soc= i$sociality.s,
	rankhi= i$rankhi,
	ranklo= i$ranklo,
	groupsize=i$groupsize.c

	) , 
warmup=1000 , iter=2000 , chains=3 , cores=3 ,types=list(adapt.delta=0.95)

)

m7 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + (ba + ba_dom)*age +  (bs + bs_dom)*soc + bg*groupsize ,
    c(a_dom,ba_dom,bs_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,ba,bs,bg) ~ dnorm(0,3),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	male= i$male,
	age=i$age.c,
	#behavior_index = i$behavior_index,
	soc= i$sociality.s,
	rankhi= i$rankhi,
	ranklo= i$ranklo,
	groupsize=i$groupsize.c

	) , 
warmup=1000 , iter=2000 , chains=3 , cores=3 ,types=list(adapt.delta=0.95)

)



write.csv(compare(m1,m2,m3,m4,m5,m6,m7)@output , file="WAICsocialinno.csv"
save(i,m1,m2,m3,m4,m5,m6,m7, file="predspredinno.RData")
save(i,m1,m2,m6, file="predspredinno.RData")



####plots
pdf("probobservedinothersage.pdf", width=8.5 , height=8.5 )

col.pal <- brewer.pal(4, "Dark2")

dom_zeros <- matrix(0,1000,length(unique(i$domain_index)))
group_zeros <- matrix(0,1000,length(unique(i$group_index)))
beh_zeros <- matrix(0,1000,length(unique(i$behavior_index)))

age.seq=seq(from=min(i$age.c), to=max(i$age.c) , length=30)

age.seq2=seq(from=min(i$age.c)-1, to=max(i$age.c)+3 , by=5)

plot(others~age.c , data=i , col="white" , ylab="probability of observing innovation in other group members" , 
	xlab="age of initial innovator" , ylim=c(0,1.1) , xaxt="n", cex.lab=1.2 , xlim=c(min(i$age.c)-2 , max(i$age.c)+3) )
for (j in 1:4){
d.pred <- list(
	domain_index=rep(j,30),
	group_index=rep(1,30),
	male=rep(mean(i$male),30),
	age=age.seq,
	soc=rep(0,30),
	rankhi= rep(mean(i$rankhi),30),
	ranklo= rep(mean(i$ranklo),30),
	groupsize=rep(0,30)
	#behavior_index=rep(1,30)
	)
replacelist <- list(a_group=group_zeros , a_dom=dom_zeros , ba_dom=dom_zeros , bs_dom=dom_zeros,
	bm_dom=dom_zeros , brhi_dom=dom_zeros , brlo_dom=dom_zeros)

ens2 <- ensemble(m1,m2,m3,m4,m5,m6,m7, n=1000 , data=d.pred, replace=list(a_group=group_zeros , bs_dom=dom_zeros,
	bm_dom=dom_zeros , brhi_dom=dom_zeros , brlo_dom=dom_zeros))

    pred.p.med <- apply(ens2$link  , 2 , median)
    pred.p.PI <- apply(ens2$link  , 2 , PI , prob=0.5)
    shade(pred.p.PI , age.seq, col=col.alpha(col.pal[j], alpha=0.1))
	lines(pred.p.med ~ age.seq , lw=3, col=col.pal[j] , lty=1)

}
domain_list <- c("foraging","investigative", "self-directed","social")
 axis(1, at = seq(from=min(i$age.c)-2 , to=max(i$age.c)+3, by = 5) , labels=seq(from=0 , to=40 , by = 5), tck=-0.02)

legend("top", inset=.01, title="Behavioral Domains", cex=1.2, bty="n",
    domain_list, fill=col.pal, horiz=TRUE)
dev.off()


###plots
pdf("probobservedinothersoc.pdf", width=8.5 , height=8.5 )

col.pal <- brewer.pal(4, "Dark2")

dom_zeros <- matrix(0,1000,length(unique(i$domain_index)))
group_zeros <- matrix(0,1000,length(unique(i$group_index)))
age.seq=seq(from=min(i$age.c), to=max(i$age.c) , length=30)
soc.seq=seq(from=min(i$sociality.s), to=max(i$sociality.s) , length=30)

plot(others~sociality.s , data=i , col="white" , ylab="probability of observing innovation in other group members" , xlab="sociality of initial innovator (standardized)", cex.lab=1.2 , ylim=c(0,1.1))
for (j in 1:4){
d.pred <- list(
	domain_index=rep(j,30),
	group_index=rep(1,30),
	male=rep(mean(i$male),30),
	age=rep(0,30),
	soc=soc.seq,
	rankhi= rep(mean(i$rankhi),30),
	ranklo= rep(mean(i$ranklo),30),
	groupsize=rep(0,30)
	)
replacelist <- list(a_group=group_zeros , a_dom=dom_zeros , ba_dom=dom_zeros , bs_dom=dom_zeros,
	bm_dom=dom_zeros , brhi_dom=dom_zeros , brlo_dom=dom_zeros)

ens2 <- ensemble(m1,m2,m3,m4,m5,m6,m7, n=1000 , data=d.pred, replace=list(a_group=group_zeros , a_dom=dom_zeros , ba_dom=dom_zeros ,
	bm_dom=dom_zeros , brhi_dom=dom_zeros , brlo_dom=dom_zeros))

    pred.p.med <- apply(ens2$link  , 2 , median)
    pred.p.PI <- apply(ens2$link  , 2 , PI , prob=0.5)
    shade(pred.p.PI , soc.seq, col=col.alpha(col.pal[j], alpha=0.1))
	lines(pred.p.med ~ soc.seq , lw=3, col=col.pal[j] , lty=1)

}
domain_list <- c("foraging","investigative", "self-directed","social")

legend("top", inset=.01, title="Behavioral Domains", cex=1.2, bty="n",
    domain_list, fill=col.pal, horiz=TRUE)
dev.off()


###plots
pdf("probobservedinothergroup.pdf", width=8.5 , height=8.5 )

col.pal <- brewer.pal(4, "Dark2")

dom_zeros <- matrix(0,1000,length(unique(i$domain_index)))
group_zeros <- matrix(0,1000,length(unique(i$group_index)))
age.seq=seq(from=min(i$age.c), to=max(i$age.c) , length=30)
soc.seq=seq(from=min(i$sociality.s), to=max(i$sociality.s) , length=30)
size.seq=seq(from=min(i$groupsize.c), to=max(i$groupsize.c) , length=30)
plot(others~groupsize.c , data=i , col="white" , ylab="probability of observing innovation in other group members" , 
	xlab="group size", cex.lab=1.2 , ylim=c(0,1.1) , xaxt="n")

d.pred <- list(
	domain_index=rep(1,30),
	group_index=rep(1,30),
	male=rep(mean(i$male),30),
	age=rep(0,30),
	soc=rep(0,30),
	rankhi= rep(mean(i$rankhi),30),
	ranklo= rep(mean(i$ranklo),30),
	groupsize=size.seq
	)
replacelist <- list(a_group=group_zeros , a_dom=dom_zeros , ba_dom=dom_zeros , bs_dom=dom_zeros,
	bm_dom=dom_zeros , brhi_dom=dom_zeros , brlo_dom=dom_zeros)

ens2 <- ensemble(m1,m2,m3,m4,m5,m6,m7, n=1000 , data=d.pred, replace=list(a_group=group_zeros , a_dom=dom_zeros , ba_dom=dom_zeros , bs_dom=dom_zeros,
	bm_dom=dom_zeros , brhi_dom=dom_zeros , brlo_dom=dom_zeros))

    pred.p.med <- apply(ens2$link  , 2 , median)
    pred.p.PI <- apply(ens2$link  , 2 , PI , prob=0.5)
    shade(pred.p.PI , size.seq, col=col.alpha("black", alpha=0.1))
	lines(pred.p.med ~ size.seq , lw=3, col="black" , lty=1)
 axis(1, at = seq(from=min(i$groupsize.c)-2 , to=max(i$groupsize.c)+3, by = 5) , labels=seq(from=5 , to=45 , by = 5), tck=-0.02)


domain_list <- c("foraging","investigative", "self-directed","social")


dev.off()



