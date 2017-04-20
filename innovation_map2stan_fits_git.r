R # load R

require(rethinking)
require(lubridate)
require(RColorBrewer)
require(Cairo)

d <- read.csv("PNASfinaldomainindividual_git.csv" , header=TRUE) #direct to local directory


###"null" model with intercepts only
 mD <- map2stan(
     alist(
         inno_count ~ dzipois( p , lambda ),
         logit(p) <- ap + logsampeff + ap_id + ap_group + ap_dom ,
         log(lambda) <- al + logsampeff + al_id + al_group + al_dom ,
         
         c(al_id,ap_id)[mono_index] ~ dmvnorm2(0, sigma_id, Rho_id),
         c(al_group,ap_group)[group_index] ~ dmvnorm2(0, sigma_group, Rho_group),
         c(al_dom,ap_dom)[domain_index] ~ dmvnorm2(0, sigma_dom, Rho_dom),
         c(ap,al)~ dstudent(2,0,1),
         c(sigma_id,sigma_group,sigma_dom) ~ dexp(1),
         c(Rho_id,Rho_group,Rho_dom) ~ dlkjcorr(3)
  
         ),
 data=list(
     inno_count=d$inno_count,
     age=d$logage.c,
     logsampeff=log(d$sampyearsum),
     mono_index=d$mono_index,
     year_index=d$year_index,
     group_index=d$group_index,
     domain_index=d$domain_index
 
 ),
     cores=3 , chains=3 , warmup=1500, iter=3500, WAIC=TRUE, types=list(adapt.delta=0.99)
 )

###age model
mD.age <- map2stan(
    alist(
        inno_count ~ dzipois( p , lambda ),
        logit(p) <- ap + logsampeff + ap_id + ap_group + ap_dom + (bagep + bagep_id + bagep_group + bagep_dom)*age,
        log(lambda) <- al + logsampeff + al_id + al_group + al_dom + (bagel + bagel_id + bagel_group + bagel_dom)*age,


        c(al_id,bagel_id,ap_id,bagep_id)[mono_index] ~ dmvnorm2(0, sigma_id, Rho_id),
        c(al_group,bagel_group,ap_group,bagep_group)[group_index] ~ dmvnorm2(0, sigma_group, Rho_group),
        c(al_dom,bagel_dom,ap_dom,bagep_dom)[domain_index] ~ dmvnorm2(0, sigma_dom, Rho_dom),
        c(ap,al)~ dstudent(2,0,1),
        c(bagep,bagel)~ dnorm(0,1),
        c(sigma_id,sigma_group,sigma_dom) ~ dexp(1),
        c(Rho_id,Rho_group,Rho_dom) ~ dlkjcorr(3)

        ),
data=list(
    inno_count=d$inno_count,
    age=d$logage.c,
    logsampeff=log(d$sampyearsum/365),
    mono_index=d$mono_index,
    year_index=d$year_index,
    group_index=d$group_index,
    domain_index=d$domain_index

),
    cores=3 , chains=3 , warmup=1500, iter=3500, WAIC=TRUE, types=list(adapt.delta=0.99)
)

##male model
mD.male <- map2stan(
    alist(
        inno_count ~ dzipois( p , lambda ),
        logit(p) <- ap + logsampeff + ap_group + ap_dom + 
        (bmalep  + bmalep_group + bmalep_dom)*male,

        log(lambda) <- al + logsampeff + al_group  + al_dom + 
        (bmalel + bmalel_group + bmalel_dom)*male,


        c(al_id,ap_id)[mono_index] ~ dmvnorm2(0, sigma_id, Rho_id),
        c(al_group,bmalel_group,ap_group,bmalep_group)[group_index] ~ dmvnorm2(0, sigma_group, Rho_group),
        c(al_dom,bsocl_dom,bagel_dom,bmalel_dom,brankhil_dom,branklol_dom,ap_dom,bsocp_dom,bagep_dom,bmalep_dom,brankhip_dom,branklop_dom)[domain_index] ~ dmvnorm2(0, sigma_dom, Rho_dom),
        c(ap,al)~ dstudent(2,0,1),
        c(bmalep,bmalel)~ dnorm(0,1),
        c(sigma_id,sigma_group,sigma_dom) ~ dexp(1),
        c(Rho_id,Rho_group,Rho_dom) ~ dlkjcorr(3)

        ),
data=list(
    inno_count=d$inno_count,
    male=d$male,
    logsampeff=log(d$sampyearsum/365),
    mono_index=d$mono_index,
    group_index=d$group_index,
    domain_index=d$domain_index

),
    cores=3 , chains=3 , warmup=1500, iter=3500, WAIC=TRUE, types=list(adapt.delta=0.99)
)

##rank model
mD.rank <- map2stan(
    alist(
        inno_count ~ dzipois( p , lambda ),
        logit(p) <- ap + logsampeff + ap_group + ap_dom + 
        (brankhip  + brankhip_group + brankhip_dom)*rankhi +
        (branklop  + branklop_group + branklop_dom)*ranklo + bsizep*groupsize,

        log(lambda) <- al + logsampeff + al_group  + al_dom + 
        (brankhil + brankhil_group + brankhil_dom)*rankhi +
        (branklol + branklol_group + branklol_dom)*ranklo + bsizel*groupsize,

        c(al_id,ap_id)[mono_index] ~ dmvnorm2(0, sigma_id, Rho_id),
        c(al_group,brankhil_group,branklol_group,ap_group,brankhip_group,branklop_group)[group_index] ~ dmvnorm2(0, sigma_group, Rho_group),
        c(al_dom,brankhil_dom,branklol_dom,ap_dom,brankhip_dom,branklop_dom)[domain_index] ~ dmvnorm2(0, sigma_dom, Rho_dom),
        c(ap,al)~ dstudent(2,0,1),
        c(brankhip,brankhil,branklop,branklol,bsizel,bsizep)~ dnorm(0,1),
        c(sigma_id,sigma_group,sigma_dom) ~ dexp(1),
        c(Rho_id,Rho_group,Rho_dom) ~ dlkjcorr(3)

        ),
data=list(
    inno_count=d$inno_count,
    rankhi=d$rankhi,
    ranklo=d$ranklo,
    logsampeff=log(d$sampyearsum/365),
    group_index=d$group_index,
    groupsize=d$groupsize.c,
    mono_index=d$mono_index,
    group_index=d$group_index,
    domain_index=d$domain_index

),
    cores=3 , chains=3 , warmup=1500, iter=3500, WAIC=TRUE, types=list(adapt.delta=0.99)
)

##sociaility model
mD.soc <- map2stan(
    alist(
        inno_count ~ dzipois( p , lambda ),
        logit(p) <- ap + logsampeff + ap_id + ap_group + ap_dom + 
        (bsocp + bsocp_id + bsocp_group + bsocp_dom)*sociality,
        log(lambda) <- al + logsampeff + al_id + al_group + al_dom + 
        (bsocl + bsocl_id + bsocl_group + bsocl_dom)*sociality,

        c(al_id,bsocl_id,ap_id,bsocp_id)[mono_index] ~ dmvnorm2(0, sigma_id, Rho_id),
        c(al_group,bsocl_group,ap_group,bsocp_group)[group_index] ~ dmvnorm2(0, sigma_group, Rho_group),
        c(al_dom,bsocl_dom,ap_dom,bsocp_dom)[domain_index] ~ dmvnorm2(0, sigma_dom, Rho_dom),
        c(ap,al)~ dstudent(2,0,1),
        c(bsocp,bsocl)~ dnorm(0,1),
        c(sigma_id,sigma_group,sigma_dom) ~ dexp(1),
        c(Rho_id,Rho_group,Rho_dom) ~ dlkjcorr(3)

        ),
data=list(
    inno_count=d$inno_count,
    sociality=d$sociality.s,
    logsampeff=log(d$sampyearsum/365),
    mono_index=d$mono_index,
    group_index=d$group_index,
    domain_index=d$domain_index

),
    cores=3 , chains=3 , warmup=1500, iter=3500, WAIC=TRUE, types=list(adapt.delta=0.99)
)

##global model
mD.ageglobal <- map2stan(
    alist(
        inno_count ~ dzipois( p , lambda ),
        logit(p) <- ap + logsampeff + ap_id + ap_group  + ap_dom + 
        (bsocp + bsocp_id  + bsocp_group + bsocp_dom)*sociality
        + (bagep + bagep_id + bagep_group + bagep_dom)*age
        + (bmalep + bmalep_group + bmalep_dom)*male +
        (brankhip  + brankhip_group + brankhip_dom)*rankhi +
        (branklop  + branklop_group + branklop_dom)*ranklo + bsizep*groupsize
,
        log(lambda) <- al + logsampeff + al_id + al_group + al_dom
        + (bagel + bagel_id  + bagel_group + bagel_dom)*age + 
        (bsocl + bsocl_id + bsocl_group + bsocl_dom)*sociality
        + (bmalel + bmalel_group + bmalel_dom)*male +
        (brankhil + brankhil_group + brankhil_dom)*rankhi +
        (branklol + branklol_group + branklol_dom)*ranklo + bsizel*groupsize
,

        c(al_id,bsocl_id,bagel_id,ap_id,bsocp_id,bagep_id)[mono_index] ~ dmvnorm2(0, sigma_id, Rho_id),
        c(al_group,bsocl_group,bagel_group,bmalel_group,brankhil_group,branklol_group,ap_group,bsocp_group,bagep_group,bmalep_group,brankhip_group,branklop_group)[group_index] ~ dmvnorm2(0, sigma_group, Rho_group),
        c(al_dom,bsocl_dom,bagel_dom,bmalel_dom,brankhil_dom,branklol_dom,ap_dom,bsocp_dom,bagep_dom,bmalep_dom,brankhip_dom,branklop_dom)[domain_index] ~ dmvnorm2(0, sigma_dom, Rho_dom),
        c(ap,al)~ dstudent(2,0,1),
        c(bsocp,bsocl,bagep,bagel,bmalep,bmalel,brankhip,brankhil,branklop,branklol,bsizel,bsizep)~ dnorm(0,1),
        c(sigma_id,sigma_group,sigma_dom) ~ dexp(1),
        c(Rho_id,Rho_group,Rho_dom) ~ dlkjcorr(3)

        ),
data=list(
    inno_count=d$inno_count,
    sociality=d$sociality.s,
    male=d$male,
    logsampeff=log(d$sampyearsum/365),
    mono_index=d$mono_index,
    age=d$logage.c,
    rankhi=d$rankhi,
    ranklo=d$ranklo,
    groupsize=d$groupsize.c,
    group_index=d$group_index,
    domain_index=d$domain_index

),
    cores=3 , chains=3 , warmup=1500, iter=3500, WAIC=TRUE, types=list(adapt.delta=0.99)
)

compare(mD,mD.age,mD.soc,mD.rank,mD.male,mD.ageglobal)
write.csv(compare(mD,mD.age,mD.soc,mD.rank,mD.male,mD.ageglobal)@output, "WAICtableinnovation.csv") #save WAIC outputs
write.csv(precis(mD.ageglobal, depth=2, digits=2  )@output, "22marglobalmodeloutputall.csv") #save model summary output for parameter estimates

#ALL GRAPHING CODE FOR MAIN MODELS BELOW#

library(rethinking)
library(lubridate)
library(RColorBrewer)

setwd("~/Dropbox/InnovationProject/Analysis")
id_zeros <- matrix(0,1000,length(unique(d$mono_index))) #zeros for varef
#year_zeros <- matrix(0,1000,length(unique(d$year))) #zeros for varef
group_zeros <- matrix(0,1000,length(unique(d$MainGroup))) #zeros for varef
domain_zeros <- matrix(0,1000,length(unique(d$domain_index))) #zeros for varef

age.seq <- seq(from=min(d$logage.c),to=max(d$logage.c),length=30) #sequencee variable of interest
soc.seq <- seq(from=min(d$sociality.s),to=max(d$sociality.s),length=30) #sequencee variable of interest
group.seq <- seq(from=min(d$groupsize.c),to=max(d$groupsize.c),length=30) #sequencee variable of interest

mono_list <- sort(unique(d$id))   

d$innovator <- 0
 for (i in 1:max(d$mono_index)){
    d$innovator <- ifelse(d$mono_index==i & sum(d$inno_count[d$mono_index==i]) > 0 , 1 , d$innovator)
 }

#display.brewer.all(colorblindFriendly=TRUE)
domain_listp <- c("a. foraging","b. investigative", "c. self-directed","d. social")
domain_listl <- c("e. foraging","f. investigative", "g. self-directed","h. social")
domain_list <- c("foraging","investigative", "self-directed","social")


####group size
d.pred <- list(
    domain_index=rep(i,30),
    mono_index=rep(1,30),
    year_index=rep(1,30),
    group_index=rep(1,30),
    sociality=rep(0,30),
    age=rep(0,30),
    groupsize=group.seq,
    male=rep(mean(d$male),30),
    rankhi=rep(mean(d$rankhi),30),
    ranklo=rep(mean(d$ranklo),30),
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros,
    ap_group=group_zeros, al_group=group_zeros,bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, 
    bagep_id=id_zeros,bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros,bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)

plot(inno_count~groupsize.c, data=d, xlab="group size centered", ylab="# of innovations per year",  
        col=col.alpha("black", alpha=0.05) , pch=19 , cex.lab=2 , ylim=c(0,5) , cex.axis=0.8 )

    pred.p.med <- apply(link2$lambda  , 2 , median)
    pred.p.PI <- apply(link2$lambda  , 2 , PI , prob=0.89)
    #shade(pred.p.PI , age.seq, col=col.alpha(col.pal[i], alpha=0.2))
    for (j in sample( c(1:1000) , 100) ){
            lines(link2$lambda[j,] ~ group.seq , lw=3, col=col.alpha(col.pal[i], alpha=0.05) , lty=1)
    }
    lines(pred.p.med ~ group.seq , lw=1, col=col.pal[i] , lty=1)

Cairo(file="Fig2.vareff_age_domain_joint.pdf" ,width=11.4 , height=11.4 , units="cm" ,type="pdf" )

col.pal <- brewer.pal(4, "Dark2")
par(mfrow = c(2,2))
par(cex = 0.4)
par(mar = c(0, 0, 0, 0), oma = c(5.1, 5.1, 0.5, 0.5))

for(i in 1:max(d$domain_index)){
        plot(inno_count~logage.c, data=d[d$domain_index==i,], xlab="log(age)", ylab="# of innovations per year",  
        col=col.alpha("white", alpha=0.05) , pch=19 , cex.lab=2.2 , ylim=c(0,1.75) ,xaxt="n", yaxt="n", cex.axis=0.8)
        title(domain_listp[i], line = -2 , cex.main=2)
####Xaxes###
if(i==1){        
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , tck=0.01, cex.axis=2, labels=FALSE )
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = .5) , labels=FALSE ,  tck=0.005, cex.axis=2)
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , tck=-0.01, cex.axis=2, labels=FALSE )
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = .5) , labels=FALSE ,  tck=-0.005, cex.axis=2)
}
if(i==2){        
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , tck=0.01, cex.axis=2, labels=FALSE )
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = .5) , labels=FALSE ,  tck=0.005, cex.axis=2)
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , tck=-0.01, cex.axis=2, labels=FALSE )
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = .5) , labels=FALSE ,  tck=-0.005, cex.axis=2)
}
if(i==3){        
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , labels=seq(from=0 , to=4 , by = 1), tck=-0.02, cex.axis=1.75)
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = .5) , labels=FALSE ,  tck=-0.01, cex.axis=1.75)
}
if(i==4){        
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , labels=seq(from=0 , to=4 , by = 1), tck=-0.02 , cex.axis=1.75)
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = .5) , labels=FALSE ,  tck=-0.01, cex.axis=1.75)
}
##yaxes
if(i==1){        
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=seq(from=0 , to=1.5 , by = .5), tck=-0.02 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=-0.01 , cex.axis=1.75)

}
if(i==3){        
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=seq(from=0 , to=1.5 , by = .5), tck=-0.02 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=-0.01 , cex.axis=1.75)
}
if(i==2){        
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=FALSE, tck=-0.01 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=FALSE, tck=0.01 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=-0.005 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=0.005 , cex.axis=1.75)
}

if(i==4){        
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=FALSE, tck=-0.01 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=FALSE, tck=0.01 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=-0.005 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=0.005 , cex.axis=1.75)
}

d.pred <- list(
    domain_index=rep(i,30),
    mono_index=rep(1,30),
    year_index=rep(1,30),
    group_index=rep(1,30),
    sociality=rep(0,30),
    age=age.seq,
    groupsize=rep(0,30),
    male=rep(mean(d$male),30),
    rankhi=rep(mean(d$rankhi),30),
    ranklo=rep(mean(d$ranklo),30),
    logsampeff=rep(mean(log(d$sampyearsum/365)),30)
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, 
    bagep_id=id_zeros,bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros,bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)
    pred.p.med <- apply(link2$lambda*(1-link2$p)  , 2 , median)
    pred.p.PI <- apply(link2$lambda*(1-link2$p)  , 2 , PI , prob=0.89)
    for (j in sample( c(1:1000) , 100) ){
            lines((1-link2$p[j,])*link2$lambda[j,] ~ age.seq , lw=3, col=col.alpha(col.pal[i], alpha=0.03) , lty=1)
    }
    lines(pred.p.med ~ age.seq , lw=2, col=col.pal[i] , lty=1)

}
    mtext("individual innovation rate", outer =TRUE, cex = 1.2, side=2, line=3)
    mtext("log(age)", outer =TRUE, cex = 1.2, side=1, line=3)
dev.off()



Cairo(file="Fig3.vareff_soc_domain_joint.pdf" ,width=11.4 , height=11.4 , units="cm" ,type="pdf" )

par(mfrow = c(2,2))
par(cex = 0.4)
par(mar = c(0, 0, 0, 0), oma = c(5.1, 5.1, 0.5, 0.5))

for(i in 1:max(d$domain_index)){
      plot(inno_count~sociality.s, data=d[d$domain_index==i,] , xlab='n', ylab='n', col=col.alpha("white" , alpha=0.00005) , pch=19 , ylim=c(0,1.75), cex.axis=0.8,xaxt="n", yaxt="n")
      title(domain_listp[i], line = -2 , cex.main=2)
####Xaxes###
if(i==1){        
        axis(1, at = seq(from=-3 , to=3, by = 1) , tck=-0.01, cex.axis=2, labels=FALSE )
        axis(1, at = seq(from=-3 , to=3, by = .5) , labels=FALSE ,  tck=-0.005, cex.axis=1.75)
        axis(1, at = seq(from=-3 , to=3, by = 1) , tck=0.01, cex.axis=2, labels=FALSE )
        axis(1, at = seq(from=-3 , to=3, by = .5) , labels=FALSE ,  tck=0.005, cex.axis=1.75)
}
if(i==2){        
         axis(1, at = seq(from=-3 , to=3, by = 1) , tck=-0.01, cex.axis=2, labels=FALSE )
        axis(1, at = seq(from=-3 , to=3, by = .5) , labels=FALSE ,  tck=-0.005, cex.axis=1.75)
        axis(1, at = seq(from=-3 , to=3, by = 1) , tck=0.01, cex.axis=2, labels=FALSE )
        axis(1, at = seq(from=-3 , to=3, by = .5) , labels=FALSE ,  tck=0.005, cex.axis=1.75)
}
if(i==3){        
        axis(1, at = seq(from=-3 , to=3, by = 1) , labels=seq(from=-3 , to=3, by = 1), tck=-0.02, cex.axis=1.75)
        axis(1, at = seq(from=-3 , to=3, by = .5) , labels=FALSE ,  tck=-0.01, cex.axis=1.75)
}
if(i==4){        
               axis(1, at = seq(from=-3 , to=3, by = 1) , labels=seq(from=-3 , to=3, by = 1), tck=-0.02, cex.axis=1.75)
        axis(1, at = seq(from=-3 , to=3, by = .5) , labels=FALSE ,  tck=-0.01, cex.axis=1.75)
}
##yaxes
if(i==1){        
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=seq(from=0 , to=1.5 , by = .5), tck=-0.02 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=-0.01 , cex.axis=1.75)
}
if(i==3){        
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=seq(from=0 , to=1.5 , by = .5), tck=-0.02 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=-0.01 , cex.axis=1.75)
}
if(i==2){        
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=FALSE, tck=-0.01 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=FALSE, tck=0.01 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=-0.005 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=0.005 , cex.axis=1.75)
}

if(i==4){        
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=FALSE, tck=-0.01 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.5, by = 0.5) , labels=FALSE, tck=0.01 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=-0.005 , cex.axis=1.75)
        axis(2, at = seq(from=0 , to=1.75, by = 0.25) , labels=FALSE, tck=0.005 , cex.axis=1.75)
}

d.pred <- list(
    domain_index=rep(i,30),
    mono_index=rep(1,30),
    year_index=rep(1,30),
    group_index=rep(1,30),
    sociality=soc.seq,
    age=rep(0,30),
    male=rep(mean(d$male),30),
    rankhi=rep(mean(d$rankhi),30),
    ranklo=rep(mean(d$ranklo),30),
    groupsize=rep(0,30),
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros  ), WAIC=TRUE)
#link2$lambda*link2$p 
    pred.p.med <- apply(link2$lambda*(1-link2$p)  , 2 , median)
    pred.p.PI <- apply(link2$lambda*(1-link2$p)  , 2 , PI , prob=0.89)
    lines(pred.p.med ~ soc.seq , lw=2, col=col.pal[i] , lty=1)
        for (j in sample( c(1:1000) , 100) ){
            lines((1-link2$p[j,])*link2$lambda[j,] ~ soc.seq , lw=3, col=col.alpha(col.pal[i], alpha=0.03) , lty=1)
        }
}

    mtext("individual innovation rate", outer =TRUE, cex = 1.2, side=2, line=3)
    mtext("sociality scores (standardized)", outer =TRUE, cex = 1.2, side=1, line=3)
dev.off()

Cairo(file="SuppFig1.vareff_age_domain_joint.pdf" ,width=17.8 , height=17.8/1.8, units="cm" ,type="pdf" )

par(mfrow = c(2,4))
par(cex = 0.5)
par(mar = c(2, 2, 2, 0), oma = c(3, 3, 0, .5))
for(i in 1:max(d$domain_index)){
        plot(innobinom~logage.c, data=d[d$domain_index==i ,] , xlab="log(age)", ylab="individual probability of innovating",  col=col.alpha("black", alpha=0.05) , pch=19 , cex.lab=2 , ylim=c(0,1.1) , xaxt='n', cex.axis=0.8)
        title(domain_listp[i], line = -1 , cex=.5)
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , labels=seq(from=0 , to=4 , by = 1), tck=-0.02, cex.axis=0.8)
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 0.5), labels=FALSE ,  tck=-0.01, cex.axis=0.8)
        #if(i=1){axis(2, at = 0 , to=5, by = 1 , labels=FALSE ,  tck=-0.01)}
        #axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , labels=seq(from=0 , to=4 , by = 1), tck=-0.02)
        #axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = .5) , labels=FALSE ,  tck=-0.01)
d.pred <- list(
    domain_index=rep(i,30),
    mono_index=rep(1,30),
    year_index=rep(1,30),
    group_index=rep(1,30),
    sociality=rep(0,30),
    age=age.seq,
    male=rep(mean(d$male),30),
    rankhi=rep(mean(d$rankhi),30),
    ranklo=rep(mean(d$ranklo),30),
    groupsize=rep(0,30),
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, 
    bagep_id=id_zeros,bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros,bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)

    pred.p.med <- apply((1-link2$p)  , 2 , median)
    pred.p.PI <- apply((1-link2$p)  , 2 , PI , prob=0.89)
    #shade(pred.p.PI , age.seq, col=col.alpha(col.pal[i], alpha=0.2))
    for (j in sample( c(1:1000) , 100) ){
            lines((1-link2$p[j,]) ~ age.seq , lw=3, col=col.alpha(col.pal[i], alpha=0.05) , lty=1)
    }
lines(pred.p.med ~ age.seq , lw=2, col=col.pal[i] , lty=1)

    #if(i==1){ mtext("probability of innovating", outer =FALSE, cex = .8, side=2, line=2)}
}
par(mar = c(2, 2, 1, 0))
for(i in 1:max(d$domain_index)){
        plot(inno_count~logage.c, data=d[d$domain_index==i & d$innovator==1,], xlab="log(age)", ylab="probability of innovating",  col=col.alpha("black", alpha=0.1) , pch=19 , cex.lab=2 , ylim=c(0,5) , xaxt="n", cex.axis=0.8)
        title(domain_listl[i], line = -1 , cex=1.6)
        #axis(1, at = seq(from=(min(d$logage.c)) , to=(exp(4)-mean(d$logage)), by = 1) , labels=seq(from=(min(d$logage)) , to=exp(4) , by = 1), tck=-0.02, cex.axis=0.8)
        #axis(1, at = seq(from=(min(d$logage.c)) , to=(exp(4)-mean(d$logage)), by = .5) , labels=FALSE ,  tck=-0.01, cex.axis=0.8)
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 1) , labels=seq(from=0 , to=4 , by = 1), tck=-0.02, cex.axis=0.8)
        axis(1, at = seq(from=(0-mean(d$logage)) , to=(4-mean(d$logage)), by = 0.5), labels=FALSE ,  tck=-0.01, cex.axis=0.8)
d.pred <- list(
    domain_index=rep(i,30),
    mono_index=rep(1,30),
    year_index=rep(1,30),
    group_index=rep(1,30),
    sociality=rep(0,30),
    age=age.seq,
    male=rep(mean(d$male),30),
    rankhi=rep(mean(d$rankhi),30),
    ranklo=rep(mean(d$ranklo),30),
    groupsize=rep(0,30),
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, 
    bagep_id=id_zeros,bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros,bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)
    pred.p.med <- apply(link2$lambda  , 2 , median)
    pred.p.PI <- apply(link2$lambda  , 2 , PI , prob=0.89)
    lines(pred.p.med ~ age.seq , lw=2, col=col.pal[i] , lty=1)
       for (j in sample( c(1:1000) , 100) ){
            lines((link2$lambda[j,]) ~ age.seq , lw=3, col=col.alpha(col.pal[i], alpha=0.05) , lty=1)
        }

}
    #mtext("# of innovations per year conditional on being an innovator", outer =TRUE, cex = 1.2, side=2, line=1)
    mtext("log(age) ", outer =TRUE, cex = 1.2, side=1, line=1)
    #mtext("# of innovations per year conditional", outer =TRUE, cex = .8, side=2, line=2 ,  at=c(0,.25))
    #mtext("on being an innovator", outer =TRUE, cex = .8, side=2, line=1 , at=c(0,.25) )
    mtext("# of innovations/year | innovator       probability of innovating       ", outer=TRUE, cex = 0.7, side=2, line=1.5 )

dev.off()


##########
Cairo(file="SuppFig2.vareff_soc_domain_joint.pdf" ,width=17.8 , height=17.8/1.8, units="cm" ,type="pdf" )

par(mfrow = c(2,4))
par(cex = 0.5)
par(mar = c(2, 2, 2, 0), oma = c(3, 3, 0, .5))
for(i in 1:max(d$domain_index)){

d.pred <- list(
    domain_index=rep(i,30),
    mono_index=rep(1,30),
    year_index=rep(1,30),
    group_index=rep(1,30),
    sociality=soc.seq,
    age=rep(0,30),
    male=rep(mean(d$male),30),
    rankhi=rep(mean(d$rankhi),30),
    ranklo=rep(mean(d$ranklo),30),
    groupsize=rep(0,30),
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros  ), WAIC=TRUE)
#link2$lambda*link2$p 
    pred.p.med <- apply((1-link2$p)  , 2 , median)
    pred.p.PI <- apply((1-link2$p)  , 2 , PI , prob=0.89)
    plot(innobinom~sociality.s, data=d[d$domain_index==i,], xlab="sociality scores (standardized)", ylab="# of innovations per year",  col=col.alpha("black", alpha=0.05) , pch=19 , cex.lab=2 , ylim=c(0,1.1), cex.axis=0.8)
    title(domain_listp[i], line = -1 , cex=.5)
    lines(pred.p.med ~ soc.seq , lw=2, col=col.pal[i] , lty=1)
    #shade(pred.p.PI , soc.seq, col=col.alpha(col.pal[i], alpha=0.2))
        for (j in sample( c(1:1000) , 100) ){
            lines((1-link2$p[j,]) ~ soc.seq , lw=3, col=col.alpha(col.pal[i], alpha=0.05) , lty=1)
        }
}
  
#dev.off()
#pdf("vareff_soc_domain_lambda.pdf", width=8.5 , height=8.5 )

for(i in 1:max(d$domain_index)){

d.pred <- list(
    domain_index=rep(i,30),
    mono_index=rep(1,30),
    year_index=rep(1,30),
    group_index=rep(1,30),
    sociality=soc.seq,
    age=rep(0,30),
    male=rep(mean(d$male),30),
    rankhi=rep(mean(d$rankhi),30),
    ranklo=rep(mean(d$ranklo),30),
    groupsize=rep(0,30),
    logsampeff=mean(log(d$sampyearsum/365))
    )


link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros  ), WAIC=TRUE)

    pred.p.med <- apply(link2$lambda  , 2 , median)
    pred.p.PI <- apply(link2$lambda  , 2 , PI , prob=0.89)
    plot(inno_count~sociality.s, data=d[d$domain_index==i & d$innovator==1,], ylim=c(0,5), xlab="sociality scores (standardized)", ylab="# of innovations per year",  col=col.alpha("black", alpha=0.05) , pch=19 , cex.lab=2, cex.axis=0.8)
    title(domain_listl[i], line = -1 , cex=.5)
    lines(pred.p.med ~ soc.seq , lw=2, col=col.pal[i] , lty=1)
    #shade(pred.p.PI , soc.seq, col=col.alpha(col.pal[i], alpha=0.2))
            for (j in sample( c(1:1000) , 100) ){
            lines(link2$lambda[j,] ~ soc.seq , lw=3, col=col.alpha(col.pal[i], alpha=0.05) , lty=1)
        }
}
    mtext("sociality scores (standardized)", outer =TRUE, cex = 1.2, side=1, line=1)
    #mtext("# of innovations per year conditional", outer =TRUE, cex = .8, side=2, line=2 ,  at=c(0,.25))
    #mtext("on being an innovator", outer =TRUE, cex = .8, side=2, line=1 , at=c(0,.25) )
    mtext("# of innovations/year | innovator       probability of innovating       ", outer=TRUE, cex = 0.7, side=2, line=1.5 )

dev.off()


##joint for male and female
Cairo(file="SuppFig3.vareff_male_domain.pdf" ,width=17.8 , height=17.8/2.9 , units="cm" ,type="pdf" )

##joint for male and female
par(mfrow=c(1,3),mar = c(2, 2, .5, 2) , oma=c(.5,3,.5,.5) , cex=0.5)

plot(1, type="n", xlab="", ylab="", xlim=c(.5, 4.5), ylim=c(0, 1) , xaxt='n', cex.axis=0.6)
points(0.5, 1 , pch="a",cex=1.2)

for(i in 1:max(d$domain_index)){

d.pred <- list(
    domain_index=c(i,i),
    mono_index=c(1,1),
    year_index=c(1,1),
    group_index=c(1,1),
    sociality=c(1,1),
    age=c(0,0),
    male=c(0,1),
    rankhi=rep(mean(d$rankhi),2),
    ranklo=rep(mean(d$ranklo),2),
    groupsize=rep(0,2),
    logsampeff=rep(mean(log(d$sampyearsum/365)),2)
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)
    pointz <- c(18,17)
    pred.p.med <- apply(link2$lambda*(1-link2$p)  , 2 , median)
    pred.p.PI <- apply(link2$lambda*(1-link2$p)  , 2 , PI , prob=0.89)

    points(i-.1,pred.p.med[1] , col=col.pal[i] , pch=pointz[1], cex=1.5)
    points(i+.1,pred.p.med[2] , col=col.pal[i] , pch=pointz[2], cex=1.5)

    lines (c(i-.1,i-.1) , pred.p.PI[,1], col=col.pal[i] )
    lines (c(i+.1,i+.1) , pred.p.PI[,2], col=col.pal[i] )

    fem <- unique(d$mono_index[d$male==0])
    mal <- unique(d$mono_index[d$male==1])


    d.pred <- list(
    domain_index=rep(i,length(fem)),
    mono_index=fem,
    year_index=rep(1,length(fem)),
    group_index=rep(1,length(fem)),
    sociality=rep(0,length(fem)),
    age=rep(0,length(fem)),
    male=rep(0,length(fem)),
    rankhi=rep(mean(d$rankhi),length(fem)),
    ranklo=rep(mean(d$ranklo),length(fem)),
    groupsize=rep(0,length(fem)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(fem))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)

    pred.p.med <- apply((link2$lambda)*(1-link2$p)   , 2 , median)
    points(rep(i-.1,length(fem)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))


d.pred <- list(
    domain_index=rep(i,length(mal)),
    mono_index=mal,
    year_index=rep(1,length(mal)),
    group_index=rep(1,length(mal)),
    sociality=rep(0,length(mal)),
    age=rep(0,length(mal)),
    male=rep(1,length(mal)),
    rankhi=rep(mean(d$rankhi),length(mal)),
    ranklo=rep(mean(d$ranklo),length(mal)),
    groupsize=rep(0,length(mal)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(mal))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)
    pred.p.med <- apply((link2$lambda)*(1-link2$p) , 2 , median)

    points(rep(i+.1,length(mal)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

    #kf <- d$inno_count[d$male==0 & d$domain_index==i]
    #points(rep(i-.1,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[1], cex=1.2)

    #km <- d$inno_count[d$male==1 & d$domain_index==i]
    #points(rep(i+.1,length(km)), km , col=col.alpha("black", 0.1) , pch=pointz[2], cex=1.2)

}

    mtext("individual innovation rate", outer =FALSE, cex = 0.6 , side=2, line=2)
domlab <- c("foraging" , "investigative" , "self-directed" , "social")
axis(1, at = c(1:4) , labels=domlab, tick=FALSE,cex.axis=0.7)

########prob of inno male female


plot(1, type="n", xlab="", ylab="", xlim=c(.5, 4.5), ylim=c(0, 1) , xaxt='n', cex.axis=0.6)
points(0.5, 1 , pch="b",cex=1.2)

for(i in 1:max(d$domain_index)){

d.pred <- list(
    domain_index=c(i,i),
    mono_index=c(1,1),
    year_index=c(1,1),
    group_index=c(1,1),
    sociality=c(1,1),
    age=c(0,0),
    male=c(0,1),
    rankhi=rep(mean(d$rankhi),2),
    ranklo=rep(mean(d$ranklo),2),
    groupsize=rep(0,2),
    logsampeff=rep(mean(log(d$sampyearsum/365)),2)
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)
    pointz <- c(18,17)
    pred.p.med <- apply((1-link2$p)  , 2 , median)
    pred.p.PI <- apply((1-link2$p)  , 2 , PI , prob=0.89)

    points(i-.1,pred.p.med[1] , col=col.pal[i] , pch=pointz[1], cex=1.5)
    points(i+.1,pred.p.med[2] , col=col.pal[i] , pch=pointz[2], cex=1.5)

    lines (c(i-.1,i-.1) , pred.p.PI[,1], col=col.pal[i] )
    lines (c(i+.1,i+.1) , pred.p.PI[,2], col=col.pal[i] )

    fem <- unique(d$mono_index[d$male==0])
    mal <- unique(d$mono_index[d$male==1])

d.pred <- list(
    domain_index=rep(i,length(fem)),
    mono_index=fem,
    year_index=rep(1,length(fem)),
    group_index=rep(1,length(fem)),
    sociality=rep(0,length(fem)),
    age=rep(0,length(fem)),
    male=rep(0,length(fem)),
    rankhi=rep(mean(d$rankhi),length(fem)),
    ranklo=rep(mean(d$ranklo),length(fem)),
    groupsize=rep(0,length(fem)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(fem))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)

    pred.p.med <- apply((1-link2$p)  , 2 , median)
    points(rep(i-.1,length(fem)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))


d.pred <- list(
    domain_index=rep(i,length(mal)),
    mono_index=mal,
    year_index=rep(1,length(mal)),
    group_index=rep(1,length(mal)),
    sociality=rep(0,length(mal)),
    age=rep(0,length(mal)),
    male=rep(1,length(mal)),
    rankhi=rep(mean(d$rankhi),length(mal)),
    ranklo=rep(mean(d$ranklo),length(mal)),
    groupsize=rep(0,length(mal)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(mal))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)
    pred.p.med <- apply((1-link2$p) , 2 , median)

    points(rep(i+.1,length(mal)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

}
    mtext("probability of innovating", outer =FALSE, cex = 0.6, side=2, line=2)
    domlab <- c("foraging" , "investigative" , "self-directed" , "social")
legend("bottom", inset=.01, title="Sex", cex=0.8,
    c("female" , "male"), pch=pointz , horiz=TRUE , bty="n")
axis(1, at = c(1:4) , labels=domlab, tick=FALSE, cex.axis=0.7)


##lambda for male and female

plot(1, type="n", xlab="", ylab="", xlim=c(.5, 4.5), ylim=c(0, 1) , xaxt='n', cex.axis=0.6)
points(0.5, 1 , pch="c",cex=1.2)

for(i in 1:max(d$domain_index)){

d.pred <- list(
    domain_index=c(i,i),
    mono_index=c(1,1),
    year_index=c(1,1),
    group_index=c(1,1),
    sociality=c(1,1),
    age=c(0,0),
    male=c(0,1),
    rankhi=rep(mean(d$rankhi),2),
    ranklo=rep(mean(d$ranklo),2),
    groupsize=rep(0,2),
    logsampeff=rep(mean(log(d$sampyearsum/365)),2)
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)
    pointz <- c(18,17)
    pred.p.med <- apply(link2$lambda  , 2 , median)
    pred.p.PI <- apply(link2$lambda  , 2 , PI , prob=0.89)

    points(i-.1,pred.p.med[1] , col=col.pal[i] , pch=pointz[1], cex=1.5)
    points(i+.1,pred.p.med[2] , col=col.pal[i] , pch=pointz[2], cex=1.5)

    lines (c(i-.1,i-.1) , pred.p.PI[,1], col=col.pal[i] )
    lines (c(i+.1,i+.1) , pred.p.PI[,2], col=col.pal[i] )

    #kf <- d$inno_count[d$male==0 & d$domain_index==i]
    #points(rep(i-.1,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[1], cex=1.2)

    #km <- d$inno_count[d$male==1 & d$domain_index==i]
    #points(rep(i+.1,length(km)), km , col=col.alpha("black", 0.1) , pch=pointz[2], cex=1.2)varef

    d.pred <- list(
    domain_index=rep(i,length(fem)),
    mono_index=fem,
    year_index=rep(1,length(fem)),
    group_index=rep(1,length(fem)),
    sociality=rep(0,length(fem)),
    age=rep(0,length(fem)),
    male=rep(0,length(fem)),
    rankhi=rep(mean(d$rankhi),length(fem)),
    ranklo=rep(mean(d$ranklo),length(fem)),
    groupsize=rep(0,length(fem)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(fem))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)

    pred.p.med <- apply((link2$lambda)  , 2 , median)
    points(rep(i-.1,length(fem)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))


d.pred <- list(
    domain_index=rep(i,length(mal)),
    mono_index=mal,
    year_index=rep(1,length(mal)),
    group_index=rep(1,length(mal)),
    sociality=rep(0,length(mal)),
    age=rep(0,length(mal)),
    male=rep(1,length(mal)),
    rankhi=rep(mean(d$rankhi),length(mal)),
    ranklo=rep(mean(d$ranklo),length(mal)),
    groupsize=rep(0,length(mal)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(mal))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros ), WAIC=TRUE)
    pred.p.med <- apply((link2$lambda) , 2 , median)

    points(rep(i+.1,length(mal)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

}
mtext("# of innovations/year | innovator", outer =FALSE, cex = 0.6, side=2, line=2)
mtext("Behavioral Domain", outer =TRUE, cex=0.8 , side=1, line=1)
domlab <- c("foraging" , "investigative" , "self-directed" , "social")
axis(1, at = c(1:4) , labels=domlab, tick=FALSE , cex.axis=0.7)

dev.off()



########RANKKK
#THINGS TO CHANGE
#LEGEND- ONLY MIDDLE, TOP, CEX=0.8
#CEX.AXIS FOR LABEL=0.7
MTEXYAXIS=CEX=0.6
XAXIS CEX=0.8
POINT SHADE ALPHA=0.03
CORBETT




########RANK

Cairo(file="SuppFig4.vareff_rank_domain.pdf" ,width=17.8 , height=17.8/2.9 , units="cm" ,type="pdf" )

##joint for male and female
par(mfrow=c(1,3),mar = c(2, 2, .5, 2) , oma=c(.5,3,.5,.5) , cex=0.5)


plot(1, type="n", xlab="", ylab="", xlim=c(.5, 4.5), ylim=c(0, 1) , xaxt='n')
points(0.5, 1 , pch="a",cex=1.2)

for(i in 1:max(d$domain_index)) {
d.predmid <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=0,
    ranklo=0,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2mid <- link(mD.ageglobal  , data=d.predmid,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

d.predlo <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=0,
    ranklo=1,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )


link2lo <- link(mD.ageglobal  , data=d.predlo,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,brankhip_dom=domain_zeros,
    brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)


d.predhi <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=1,
    ranklo=0,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2hi <- link(mD.ageglobal  , data=d.predhi,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    branklol_dom=domain_zeros,branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

 pointz <- c(17,18,19)
 #mid
    pred.p.med.mid <- apply(link2mid$lambda*(1-link2mid$p)  , 2 , median)
    pred.p.PI.mid <- apply(link2mid$lambda*(1-link2mid$p)   , 2 , PI , prob=0.89)
    points(i,pred.p.med.mid , col=col.pal[i] , pch=pointz[2], cex=1.5)
    lines (c(i,i) , pred.p.PI.mid, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==0 & d$ranklo==0 & d$domain_index==i]
    #points(rep(i,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[2], cex=1)
##lo
    pred.p.med.lo <- apply(link2lo$lambda*(1-link2lo$p)   , 2 , median)
    pred.p.PI.lo <- apply(link2lo$lambda*(1-link2mid$p)   , 2 , PI , prob=0.89)
    points(i-.2,pred.p.med.lo , col=col.pal[i] , pch=pointz[1], cex=1.5)
    lines (c(i-.2,i-.2) , pred.p.PI.lo, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==0 & d$ranklo==1 & d$domain_index==i]
    #points(rep(i-.2,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[1], cex=1)

##hi
    pred.p.med.hi <- apply(link2hi$lambda*(1-link2hi$p)   , 2 , median)
    pred.p.PI.hi <- apply(link2hi$lambda*(1-link2hi$p)   , 2 , PI , prob=0.89)
    points(i+.2,pred.p.med.hi , col=col.pal[i] , pch=pointz[3], cex=1.5)
    lines (c(i+.2,i+.2) , pred.p.PI.hi, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==1 & d$ranklo==0 & d$domain_index==i]
    #points(rep(i+.2,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[3], cex=1)


    hi <- unique(d$mono_index[d$rankhi==1])
    
    d.pred <- list(
    domain_index=rep(i,length(hi)),
    mono_index=hi,
    year_index=rep(1,length(hi)),
    group_index=rep(1,length(hi)),
    sociality=rep(0,length(hi)),
    age=rep(0,length(hi)),
    male=rep(mean(d$male),length(hi)),
    rankhi=rep(1,length(hi)),
    ranklo=rep(0,length(hi)),
    groupsize=rep(0,length(hi)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(hi))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    branklol_dom=domain_zeros,branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

    pred.p.med <- apply((link2$lambda)*(1-link2$p)  , 2 , median)
    points(rep(i+.2,length(hi)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

    lo <- unique(d$mono_index[d$ranklo==1])
    
    d.pred <- list(
    domain_index=rep(i,length(lo)),
    mono_index=lo,
    year_index=rep(1,length(lo)),
    group_index=rep(1,length(lo)),
    sociality=rep(0,length(lo)),
    age=rep(0,length(lo)),
    male=rep(mean(d$male),length(lo)),
    rankhi=rep(1,length(lo)),
    ranklo=rep(0,length(lo)),
    groupsize=rep(0,length(lo)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(lo))
    )
link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,brankhip_dom=domain_zeros,
    brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)


    pred.p.med <- apply((link2$lambda)*(1-link2$p)  , 2 , median)
    points(rep(i-.2,length(lo)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

    mid <- unique(d$mono_index[d$rankmid==1])
    
    d.pred <- list(
    domain_index=rep(i,length(mid)),
    mono_index=mid,
    year_index=rep(1,length(mid)),
    group_index=rep(1,length(mid)),
    sociality=rep(0,length(mid)),
    age=rep(0,length(mid)),
    male=rep(mean(d$male),length(mid)),
    rankhi=rep(0,length(mid)),
    ranklo=rep(0,length(mid)),
    groupsize=rep(0,length(mid)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(mid))
    )


link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)


    pred.p.med <- apply((link2$lambda)*(1-link2$p)  , 2 , median)
    points(rep(i,length(mid)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))


}

    mtext("individual innovation rate", outer =FALSE, cex = 0.6 , side=2, line=2)
domlab <- c("foraging" , "investigative" , "self-directed" , "social")
axis(1, at = c(1:4) , labels=domlab, tick=FALSE,cex.axis=0.7)


plot(1, type="n", xlab="", ylab="", xlim=c(.5, 4.5), ylim=c(0, 1) , xaxt='n')
points(0.5, 1 , pch="b",cex=1.2)

for(i in 1:max(d$domain_index)) {
d.predmid <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=0,
    ranklo=0,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2mid <- link(mD.ageglobal  , data=d.predmid,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

d.predlo <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=0,
    ranklo=1,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )


link2lo <- link(mD.ageglobal  , data=d.predlo,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,brankhip_dom=domain_zeros,
    brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)


d.predhi <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=1,
    ranklo=0,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2hi <- link(mD.ageglobal  , data=d.predhi,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    branklol_dom=domain_zeros,branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

 pointz <- c(17,18,19)
 #mid
    pred.p.med.mid <- apply((1-link2mid$p)  , 2 , median)
    pred.p.PI.mid <- apply((1-link2mid$p)   , 2 , PI , prob=0.89)
    points(i,pred.p.med.mid , col=col.pal[i] , pch=pointz[2], cex=1.5)
    lines (c(i,i) , pred.p.PI.mid, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==0 & d$ranklo==0 & d$domain_index==i]
  #  points(rep(i,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[2], cex=1)
##lo
    pred.p.med.lo <- apply((1-link2lo$p)   , 2 , median)
    pred.p.PI.lo <- apply((1-link2mid$p)   , 2 , PI , prob=0.89)
    points(i-.2,pred.p.med.lo , col=col.pal[i] , pch=pointz[1], cex=1.5)
    lines (c(i-.2,i-.2) , pred.p.PI.lo, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==0 & d$ranklo==1 & d$domain_index==i]
   # points(rep(i-.2,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[1], cex=1)

##hi
    pred.p.med.hi <- apply((1-link2hi$p)   , 2 , median)
    pred.p.PI.hi <- apply((1-link2hi$p)   , 2 , PI , prob=0.89)
    points(i+.2,pred.p.med.hi , col=col.pal[i] , pch=pointz[3], cex=1.5)
    lines (c(i+.2,i+.2) , pred.p.PI.hi, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==1 & d$ranklo==0 & d$domain_index==i]
   # points(rep(i+.2,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[3], cex=1)

    hi <- unique(d$mono_index[d$rankhi==1])
    
    d.pred <- list(
    domain_index=rep(i,length(hi)),
    mono_index=hi,
    year_index=rep(1,length(hi)),
    group_index=rep(1,length(hi)),
    sociality=rep(0,length(hi)),
    age=rep(0,length(hi)),
    male=rep(mean(d$male),length(hi)),
    rankhi=rep(1,length(hi)),
    ranklo=rep(0,length(hi)),
    groupsize=rep(0,length(hi)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(hi))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    branklol_dom=domain_zeros,branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

    pred.p.med <- apply((1-link2$p)  , 2 , median)
    points(rep(i+.2,length(hi)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

    lo <- unique(d$mono_index[d$ranklo==1])
    
    d.pred <- list(
    domain_index=rep(i,length(lo)),
    mono_index=lo,
    year_index=rep(1,length(lo)),
    group_index=rep(1,length(lo)),
    sociality=rep(0,length(lo)),
    age=rep(0,length(lo)),
    male=rep(mean(d$male),length(lo)),
    rankhi=rep(0,length(lo)),
    ranklo=rep(mean(d$ranklo),length(lo)),
    groupsize=rep(0,length(lo)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(lo))
    )
link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,brankhip_dom=domain_zeros,
    brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)


    pred.p.med <- apply((1-link2$p)  , 2 , median)
    points(rep(i-.2,length(lo)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

    mid <- unique(d$mono_index[d$rankmid==1])
    
    d.pred <- list(
    domain_index=rep(i,length(mid)),
    mono_index=mid,
    year_index=rep(1,length(mid)),
    group_index=rep(1,length(mid)),
    sociality=rep(0,length(mid)),
    age=rep(0,length(mid)),
    male=rep(mean(d$male),length(mid)),
    rankhi=rep(0,length(mid)),
    ranklo=rep(0,length(mid)),
    groupsize=rep(0,length(mid)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(mid))
    )


link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)


    pred.p.med <- apply((1-link2$p)  , 2 , median)
    points(rep(i,length(mid)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

}

    mtext("probability of innovating", outer =FALSE, cex = 0.6, side=2, line=2)
    domlab <- c("foraging" , "investigative" , "self-directed" , "social")
axis(1, at = c(1:4) , labels=domlab, tick=FALSE, cex.axis=0.7)

legend("bottom", inset=.01, title="Rank", cex=1,
    c("low" , "mid" , "hi"), pch=pointz , horiz=TRUE , bty="n")



plot(1, type="n", xlab="", ylab="", xlim=c(.5, 4.5), ylim=c(0, 1) , xaxt='n')
points(0.5, 1 , pch="c",cex=1.2)

for(i in 1:max(d$domain_index)) {
d.predmid <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=0,
    ranklo=0,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2mid <- link(mD.ageglobal  , data=d.predmid,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

d.predlo <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=0,
    ranklo=1,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )


link2lo <- link(mD.ageglobal  , data=d.predlo,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,brankhip_dom=domain_zeros,
    brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)


d.predhi <- list(
    domain_index=i,
    mono_index=1,
    year_index=1,
    group_index=1,
    sociality=0,
    age=0,
    male=mean(d$male),
    rankhi=1,
    ranklo=0,
    groupsize=0,
    logsampeff=mean(log(d$sampyearsum/365))
    )

link2hi <- link(mD.ageglobal  , data=d.predhi,replace=list(ap_id=id_zeros, al_id=id_zeros,
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    branklol_dom=domain_zeros,branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

 pointz <- c(17,18,19)
 #mid
    pred.p.med.mid <- apply(link2mid$lambda  , 2 , median)
    pred.p.PI.mid <- apply(link2mid$lambda   , 2 , PI , prob=0.89)
    points(i,pred.p.med.mid , col=col.pal[i] , pch=pointz[2], cex=1.5)
    lines (c(i,i) , pred.p.PI.mid, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==0 & d$ranklo==0 & d$domain_index==i]
    #points(rep(i,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[2], cex=1)
##lo
    pred.p.med.lo <- apply(link2lo$lambda  , 2 , median)
    pred.p.PI.lo <- apply(link2lo$lambda   , 2 , PI , prob=0.89)
    points(i-.2,pred.p.med.lo , col=col.pal[i] , pch=pointz[1], cex=1.5)
    lines (c(i-.2,i-.2) , pred.p.PI.lo, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==0 & d$ranklo==1 & d$domain_index==i]
    #points(rep(i-.2,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[1], cex=1)

##hi
    pred.p.med.hi <- apply(link2hi$lambda  , 2 , median)
    pred.p.PI.hi <- apply(link2hi$lambda   , 2 , PI , prob=0.89)
    points(i+.2,pred.p.med.hi , col=col.pal[i] , pch=pointz[3], cex=1.5)
    lines (c(i+.2,i+.2) , pred.p.PI.hi, col=col.pal[i] )
    kf <- d$inno_count[d$rankhi==1 & d$ranklo==0 & d$domain_index==i]
    #points(rep(i+.2,length(kf)), kf , col=col.alpha("black", 0.1) , pch=pointz[3], cex=1)


    hi <- unique(d$mono_index[d$rankhi==1])
    
    d.pred <- list(
    domain_index=rep(i,length(hi)),
    mono_index=hi,
    year_index=rep(1,length(hi)),
    group_index=rep(1,length(hi)),
    sociality=rep(0,length(hi)),
    age=rep(0,length(hi)),
    male=rep(mean(d$male),length(hi)),
    rankhi=rep(1,length(hi)),
    ranklo=rep(0,length(hi)),
    groupsize=rep(0,length(hi)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(hi))
    )

link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    branklol_dom=domain_zeros,branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

    pred.p.med <- apply((link2$lambda)  , 2 , median)
    points(rep(i+.2,length(hi)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

    lo <- unique(d$mono_index[d$ranklo==1])
    
    d.pred <- list(
    domain_index=rep(i,length(lo)),
    mono_index=lo,
    year_index=rep(1,length(lo)),
    group_index=rep(1,length(lo)),
    sociality=rep(0,length(lo)),
    age=rep(0,length(lo)),
    male=rep(mean(d$male),length(lo)),
    rankhi=rep(1,length(lo)),
    ranklo=rep(0,length(lo)),
    groupsize=rep(0,length(lo)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(lo))
    )
link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,brankhip_dom=domain_zeros,
    brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)


    pred.p.med <- apply((link2$lambda)  , 2 , median)
    points(rep(i-.2,length(lo)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))

    mid <- unique(d$mono_index[d$rankmid==1])
    
    d.pred <- list(
    domain_index=rep(i,length(mid)),
    mono_index=mid,
    year_index=rep(1,length(mid)),
    group_index=rep(1,length(mid)),
    sociality=rep(0,length(mid)),
    age=rep(0,length(mid)),
    male=rep(mean(d$male),length(mid)),
    rankhi=rep(0,length(mid)),
    ranklo=rep(0,length(mid)),
    groupsize=rep(0,length(mid)),
    logsampeff=rep(mean(log(d$sampyearsum/365)) , length(mid))
    )


link2 <- link(mD.ageglobal  , data=d.pred,replace=list(
    ap_group=group_zeros, al_group=group_zeros ,
    bagel_id=id_zeros, bsocl_id=id_zeros, bsocp_id=id_zeros, bagep_id=id_zeros, 
    bagel_group=group_zeros,bagep_group=group_zeros, bsocp_group=group_zeros,bmalep_group=group_zeros,
    bsocl_group=group_zeros,bmalel_group=group_zeros, bagep_dom=domain_zeros,bagel_dom=domain_zeros, 
    bsocp_dom=domain_zeros, bsocl_dom=domain_zeros, bmalep_dom=domain_zeros, bmalel_dom=domain_zeros, 
    brankhil_dom=domain_zeros,branklol_dom=domain_zeros,brankhip_dom=domain_zeros,
    branklop_dom=domain_zeros, brankhil_group=group_zeros,branklol_group=group_zeros,
    brankhip_group=group_zeros,branklop_group=group_zeros), WAIC=TRUE)

    pred.p.med <- apply((link2$lambda)  , 2 , median)
    points(rep(i,length(mid)),pred.p.med , pch=19 , col=col.alpha(col.pal[i], 0.03))


}

mtext("# of innovations/year | innovator", outer =FALSE, cex = 0.6, side=2, line=2)
mtext("Behavioral Domain", outer =TRUE, cex=0.8 , side=1, line=1)
domlab <- c("foraging" , "investigative" , "self-directed" , "social")
axis(1, at = c(1:4) , labels=domlab, tick=FALSE , cex.axis=0.7)

dev.off()
