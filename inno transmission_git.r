R
library(rethinking)
library(lubridate)
library(RColorBrewer)


d <- read.csv("PNASfinaldomainindividual_git.csv" , header=TRUE) #direct to local directory
inno <- read.csv("PNASinnovationlistused_git.csv" , header=TRUE)


###innoused

 dsoc=subset(d, select=c("id","year","male","sampyearsum" , "sociality" , "age" , 
    "rank" , "max_group" , "rank_class", "sociality.s", "logage" , "logage.c", "adult", "rankhi",
    "ranklo","groupsize.c" , "inno_count"))

 dsoc <- dsoc[dsoc$inno_count>0,]

# dsoc=subset(dsoc, select=c("id","year","male","sampyearsum" , "sociality" , "age" , "age.c" , 
#    "rank" , "max_group" , "rank_class", "sociality.s", "logage" , "logage.c", "adult", "rankhi",
#    "ranklo","groupsize.c"  ))
 
dsoc <- droplevels(dsoc)
innosub <- subset(inno, select=c("group","date","id","behav_code","domain","date2","year",
    "domain_index","behavior_index" ,"others"))

 innosub$YEARmatch <- innosub$year
 df <- merge( dsoc,innosub , by=c("year","id") )

sort(unique(inno$id)) %in% sort(unique(d$id[d$inno_count>0]))
innosub$group_index <- as.integer(innosub$group)
innosub$mono_index <- as.integer(innosub$id)
innosub$sociality.s <- innosub$male <-innosub$rankhi <- innosub$ranklo <- innosub$groupsize.c <- innosub$logage.c <- innosub$logage <- innosub$age<-0
for(i in unique(innosub$id)){

	for(j in 2007:2011){
	innosub$sociality.s <- ifelse(innosub$id==i & innosub$year==j , dsoc$sociality.s[dsoc$id==i & dsoc$year==j] , innosub$sociality.s)
	innosub$logage.c  <- ifelse(innosub$id==i & innosub$year==j , dsoc$logage.c[dsoc$id==i & dsoc$year==j] , innosub$logage.c )
	innosub$logage  <- ifelse(innosub$id==i & innosub$year==j , dsoc$logage[dsoc$id==i & dsoc$year==j] , innosub$logage)
	innosub$rankhi <- ifelse(innosub$id==i & innosub$year==j , dsoc$rankhi[dsoc$id==i & dsoc$year==j] , innosub$rankhi)
	innosub$ranklo <- ifelse(innosub$id==i & innosub$year==j , dsoc$ranklo[dsoc$id==i & dsoc$year==j] , innosub$ranklo)
	innosub$groupsize.c <- ifelse(innosub$id==i & innosub$year==j , dsoc$groupsize.c[dsoc$id==i & dsoc$year==j] , innosub$groupsize.c)
	innosub$age <- ifelse(innosub$id==i & innosub$year==j , dsoc$age[dsoc$id==i & dsoc$year==j] , innosub$age)
	innosub$groupsize <- ifelse(innosub$id==i & innosub$year==j , dsoc$groupsize[dsoc$id==i & dsoc$year==j] , innosub$groupsize)

	}
		innosub$male <- ifelse(innosub$id==i , dsoc$male[dsoc$id==i ] , innosub$male)

}

i <- innosub[complete.cases(innosub),]
i$behavior_index <- as.integer(as.factor(i$behavior))
i$group_index <- as.integer(i$group)
i$domain_index <- as.integer(i$domain)
i$age.c=i$age - mean(i$age)
#write.csv(i,"PNASinnovationlistused_git.csv")

 m1 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + bg*groupsize ,

    a_dom[domain_index] ~ dnorm( 0 , sigma_dom ),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
	c(a,bg) ~ dnorm(0,1),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	groupsize=i$groupsize.c

 

	) , 
warmup=2000 , iter=4000 , chains=3 , cores=3 ,types=list(adapt.delta=0.99)

)

m2 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + (ba + ba_dom)*age + bg*groupsize,

    c(a_dom,ba_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,ba,bg) ~ dnorm(0,1),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	#behavior_index = i$behavior_index,
	age = i$logage.c,
	groupsize=i$groupsize.c


	) , 
warmup=2000 , iter=4000 , chains=3 , cores=3 ,types=list(adapt.delta=0.99)

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
warmup=2000 , iter=4000 , chains=2 , cores=2 ,types=list(adapt.delta=0.99)

)

m4 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group  + (bm + bm_dom)*male + bg*groupsize,

    c(a_dom,bm_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),
	c(a,bm,bg) ~ dnorm(0,1),
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
warmup=2000 , iter=4000 , chains=3 , cores=3 ,types=list(adapt.delta=0.99)

)

m5 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + (brhi +brhi_dom)*rankhi + (brlo + brlo_dom)*ranklo + bg*groupsize,

    c(a_dom,brhi_dom,brlo_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,brhi,brlo,bg) ~ dnorm(0,1),
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
warmup=2000 , iter=4000 , chains=3 , cores=3 ,types=list(adapt.delta=0.99)

)

m6 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + (ba + ba_dom)*age +  (bs + bs_dom)*soc +
    (bm + bm_dom)*male + (brhi +brhi_dom)*rankhi + (brlo + brlo_dom)*ranklo + bg*groupsize ,
    c(a_dom,ba_dom,bs_dom,bm_dom,brhi_dom,brlo_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,ba,bs,bm,brhi,brlo,bg) ~ dnorm(0,1),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	male= i$male,
	age=i$logage.c,
	#behavior_index = i$behavior_index,
	soc= i$sociality.s,
	rankhi= i$rankhi,
	ranklo= i$ranklo,
	groupsize=i$groupsize.c

	) , 
warmup=2000 , iter=4000 , chains=3 , cores=3 ,types=list(adapt.delta=0.99)

)

m7 <- map2stan(
    alist(
    others ~ dbinom( 1,p ),
    logit(p) <- a + a_dom + a_group + (ba + ba_dom)*age +  (bs + bs_dom)*soc + bg*groupsize ,
    c(a_dom,ba_dom,bs_dom)[domain_index] ~  dmvnorm2(0, sigma_dom, Rho_dom),
    a_group[group_index] ~ dnorm( 0 , sigma_group ),
    #a_beh[behavior_index] ~ dnorm( 0 , sigma_beh ),

	c(a,ba,bs,bg) ~ dnorm(0,1),
    c(sigma_dom,sigma_group) ~ dcauchy(0,2),
    Rho_dom ~ dlkjcorr(3)
),

data=list(
	others = i$others,
	domain_index = i$domain_index,
	group_index = i$group_index,
	male= i$male,
	age=i$logage.c,
	#behavior_index = i$behavior_index,
	soc= i$sociality.s,
	rankhi= i$rankhi,
	ranklo= i$ranklo,
	groupsize=i$groupsize.c

	) , 
warmup=2000 , iter=4000 , chains=3 , cores=3 ,types=list(adapt.delta=0.99)

)



write.csv(compare(m1,m2,m3,m4,m5,m6,m7)@output , file="WAICsocialinno2.csv")
save(i,m1,m2,m3,m4,m5,m6,m7, file="predspredinno2.RData")



####plots
pdf("probobservedinothersage2.pdf", width=8.5 , height=8.5 )

col.pal <- brewer.pal(4, "Dark2")

dom_zeros <- matrix(0,1000,length(unique(i$domain_index)))
group_zeros <- matrix(0,1000,length(unique(i$group_index)))
beh_zeros <- matrix(0,1000,length(unique(i$behavior_index)))

age.seq=seq(from=min(i$logage.c), to=max(i$logage.c) , length=30)

age.seq2=seq(from=min(i$logage.c)-1, to=max(i$logage.c)+2 , by=5)

plot(others~logage.c , data=i , col="white" , ylab="probability of observing innovation in other group members" , 
	xlab="log(age) of initial innovator" , ylim=c(0,1.1) , cex.lab=1.2 , xaxt='n' )
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
axis(1, at = seq(from=0-mean(i$logage) , to=4-mean(i$logage) , by = 1) , labels=seq(from=0 , to=4 , by = 1), tck=-0.02)

legend("top", inset=.01, title="Behavioral Domains", cex=1.2, bty="n",
    domain_list, fill=col.pal, horiz=TRUE)
dev.off()


########plot all ############
pdf("probobservedinothers.pdf", width=8.5 , height=3.25 )
par(mfrow=c(1,3) , mar=c(2.5,.1,.1,.1) , oma=c(2,3.5,.5,.5)  )
col.pal <- brewer.pal(4, "Dark2")

col.pal <- brewer.pal(4, "Dark2")

dom_zeros <- matrix(0,1000,length(unique(i$domain_index)))
group_zeros <- matrix(0,1000,length(unique(i$group_index)))
beh_zeros <- matrix(0,1000,length(unique(i$behavior_index)))

age.seq=seq(from=min(i$logage.c), to=max(i$logage.c) , length=30)

age.seq2=seq(from=min(i$logage.c)-1, to=max(i$logage.c)+2 , by=5)

plot(others~logage.c , data=i , col="white" , ylab="probability of observing innovation in other group members" , 
	xlab="log(age) of initial innovator" , ylim=c(0,1.1) , cex.lab=1 , xaxt='n' )
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
	lines(pred.p.med ~ age.seq , lw=2, col=col.pal[j] , lty=1)

}
domain_list <- c("foraging","investigative", "self-directed","social")
axis(1, at = seq(from=0-mean(i$logage) , to=4-mean(i$logage) , by = 1) , labels=seq(from=0 , to=4 , by = 1), tck=-0.02)
axis(1, at = seq(from=0-mean(i$logage) , to=4-mean(i$logage) , by = 0.5) , labels=FALSE, tck=-0.01)

legend("top", inset=.05, title="Behavioral Domains", cex=.7, bty="n",
domain_list, fill=col.pal, horiz=TRUE)
mtext( text="log(age) of initial innovator" , side=1, line=2.5, outer=FALSE , cex=.9)
legend("topleft", "a.", bty="n") 


###plots

col.pal <- brewer.pal(4, "Dark2")

dom_zeros <- matrix(0,1000,length(unique(i$domain_index)))
group_zeros <- matrix(0,1000,length(unique(i$group_index)))
age.seq=seq(from=min(i$age.c), to=max(i$age.c) , length=30)
soc.seq=seq(from=min(i$sociality.s), to=max(i$sociality.s) , length=30)

plot(others~sociality.s , data=i , col="white" , yaxt='n' , xlab="sociality of initial innovator (standardized)", cex.lab=1 , ylim=c(0,1.1))
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
	lines(pred.p.med ~ soc.seq , lw=2, col=col.pal[j] , lty=1)

}
domain_list <- c("foraging","investigative", "self-directed","social")

legend("top", title="Behavioral Domains", cex=.7, bty="n",
    domain_list, fill=col.pal, horiz=TRUE, inset=0.05)
mtext( text="sociality of initial innovator (standardized)" , side=1, line=2.5, outer=FALSE , cex=.9)
legend("topleft", "b.", bty="n") 

###plot group size

col.pal <- brewer.pal(4, "Dark2")

dom_zeros <- matrix(0,1000,length(unique(i$domain_index)))
group_zeros <- matrix(0,1000,length(unique(i$group_index)))
age.seq=seq(from=min(i$age.c), to=max(i$age.c) , length=30)
soc.seq=seq(from=min(i$sociality.s), to=max(i$sociality.s) , length=30)
size.seq=seq(from=min(i$groupsize.c), to=max(i$groupsize.c) , length=30)
plot(others~groupsize.c , data=i , col="white" , yaxt='n' , 
	xlab="group size", cex.lab=1 , ylim=c(0,1.1) , xaxt="n")

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
	lines(pred.p.med ~ size.seq , lw=2, col="black" , lty=1)
 axis(1, at = seq(from=min(i$groupsize.c)-2 , to=max(i$groupsize.c)+3, by = 5) , labels=seq(from=5 , to=45 , by = 5), tck=-0.02)
legend("topleft", "c.", bty="n") 

mtext( text="group size" , side=1, line=2.5, outer=FALSE , cex=.9)


mtext( text="probability of observing innovation in others" , side=2, line=2, outer=TRUE , cex=.9)

dev.off()
