#
# Fitting NCC data using R's msm package with bootstrap
# on correlated data
#
# refer to the toturial for details
#
# 2017-10-15 Hongbin Zhang
# 2018-02-12 including location, age, gender interaction (without 95% CI, i.e., no bootstrap yet)
# 2018-02-20 bootstrap interaction
# 2018-02-21 fix a bug on age interaction
# 2018-02-23 add LR likelihood ratio test on interaction term
# 2018-02-25 add descriptive for treatment under the final 221 paths
# 2018-02-25 follow-up summary
# 2018-02-27 age range, average path length
# 2018-03-02 path summary; modified the plot labels, etc.
# 2018-03-05 merge a file with number of cyst variable, try to adjust that

library(msm)
library(haven)

raw <- read.csv("../data/NCClong.csv", header=T)[,c(1:7,44)]
raw <- raw[order(raw[,1], raw[,2], raw[,3]), ]

dim(raw)
names(raw)
raw.sub <- subset(raw, Tissue != "Parenchymal")
str(raw)

location <- paste(as.character(raw.sub$Region), as.character(raw.sub$Area), as.character(raw.sub$Tissue), sep="-")

table(location)

raw1 <- subset(raw, drug.x == 1)
raw0 <- subset(raw, drug.x == 0 )

myraw <- raw1
myraw <- raw0

table(subset(myraw, MONTH==0)$STATUS)
sum(table(subset(myraw, MONTH==0)$STATUS))
prop.table(table(subset(myraw, MONTH==0)$STATUS))
length(unique(subset(myraw, MONTH==0)$ID))

table(subset(myraw, MONTH==1)$STATUS)
sum(table(subset(myraw, MONTH==1)$STATUS))
prop.table(table(subset(myraw, MONTH==1)$STATUS))
length(unique(subset(myraw, MONTH==1)$ID))
sum(!is.na(subset(myraw, MONTH==1)$STATUS))
length(unique(subset(myraw, MONTH==1 & !is.na(STATUS))$ID))


table(subset(myraw, MONTH==6)$STATUS)
sum(table(subset(myraw, MONTH==6)$STATUS))
prop.table(table(subset(myraw, MONTH==6)$STATUS))
length(unique(subset(myraw, MONTH==12 & !is.na(STATUS))$ID))


table(subset(myraw, MONTH==12)$STATUS)
sum(table(subset(myraw, MONTH==12)$STATUS))
prop.table(table(subset(myraw, MONTH==12)$STATUS))
length(unique(subset(myraw, MONTH==12 & !is.na(STATUS))$ID))

table(subset(myraw, MONTH==24)$STATUS)
sum(table(subset(myraw, MONTH==24)$STATUS))
prop.table(table(subset(myraw, MONTH==24)$STATUS))
length(unique(subset(myraw, MONTH==24 & !is.na(STATUS))$ID))



raw$newid <- paste(raw$ID, raw$CombCode, sep="-")
raw$drug.x <- ifelse(raw$drug.x==1, 1, 0)

length(unique(raw$ID))
length(unique(raw$CombCode))

ids <- unique(raw$newid)
len <- length(ids); len


ncc <- NULL

c.drop <- 0

for(i in 1:len){
   ttt <- subset(raw, newid %in% ids[i])
   OO  <- ttt$STATUS
   # remove the trailing '4'
   flag <- 0
   for(j in 1:dim(ttt)[1]){
      if(!is.na(OO[j])){
         if(OO[j] == 4){
            flag <- 1;
            break;
         }
      }
   }
   if(j==dim(ttt)[1]){
      tt <- ttt
   }else{
      tt <- ttt[1:j,]
   }

   O  <- ttt$STATUS
   iNA<- sum(is.na(O))
   nv <- dim(tt)[1]
   if(iNA == 0){
      ncc <- rbind(ncc, tt)
   }else{
      if( ! (is.na(O[nv])) ){
         ttt <- subset(tt, !is.na(STATUS) )
         ncc <- rbind(ncc, ttt)
         cat("intermittee missing, i=", i, "nv=", dim(tt)[1], "iNA=", iNA, "\n");
      }else{
         for(j in (nv-1):1){
            if( !is.na(O[j]) )
               break
         }
         if(j == nv - 1){
             miss <- nv
             cat("drop out 1, i=", i, "nv=", dim(tt)[1], "iNA=", iNA, "miss=", miss, "\n");
             tt$STATUS[nv] <- 99
             ttt <- subset(tt, !is.na(STATUS) )
             ncc <- rbind(ncc, ttt)
             c.drop <- c.drop + 1
         }else{
             miss <- seq(j+1, nv)
             cat("drop out 2, i=", i, "nv=", dim(tt)[1], "iNA=", iNA, "miss=", miss, "\n");
             tt$STATUS[j+1] <- 99
             ttt <- subset(tt, !is.na(STATUS) )
             ncc <- rbind(ncc, ttt)
             c.drop <- c.drop + 1
         }

      }
   }
}
cat("c.drop=", c.drop, "ratio=", c.drop/len, "\n");
#write.csv(ncc, "ncc.csv", row.names=F)

dim(raw); dim(ncc);

table(subset(ncc, MONTH==0)$STATUS)

table(subset(ncc, MONTH==1)$STATUS)

table(subset(ncc, MONTH==6)$STATUS)

table(subset(ncc, MONTH==12)$STATUS)

table(subset(ncc, MONTH==24)$STATUS)

head(ncc, 20)

length(unique(ncc$ID))
length(unique(ncc$newid))

length(unique(subset(ncc, drug.x==1)$ID))
# 62

#62/117  53%

length(unique(subset(ncc, drug.x==0)$ID))
# 55

# summary the paths

ids <- unique(ncc$newid)
len <- length(ids); len

ncc.tmp <- NULL
for(i in 1:len){
  tt <- subset(ncc, newid %in% ids[i])
  nv <- dim(tt)[1]
  STATE <- tt$STATUS
  c.STATE <- as.character(STATE[1])
  MON  <- tt$MONTH
  c.MON <- as.character(MON[1])
  for(j in 2:nv){
    c.STATE <- paste(c.STATE, as.character(STATE[j]), sep="-")
    c.MON   <- paste(c.MON, as.character(MON[j]), sep="-")
  }
  tt$c.STATE <- c.STATE
  tt$c.MON     <- c.MON
  first <- numeric(nv)
  first[1] <- 1
  tt$first <- first
  ncc.tmp <- rbind(ncc.tmp, tt)

}
tmp <- subset(ncc.tmp, first==1)
table(tmp$c.STATE)
cbind(tmp$c.MON, tmp$c.STATE)

tmp.drug1 <- subset(tmp, drug.x == 1)
table(tmp.drug1$c.STATE)


statetable.msm(STATUS,newid, data=subset(ncc, STATUS!=99))
sum(statetable.msm(STATUS,newid, data=subset(ncc, STATUS!=99)) )
(treated <- statetable.msm(STATUS,newid, data=subset(ncc, STATUS!= 99 & drug.x==1)))
(placebo <- statetable.msm(STATUS,newid, data=subset(ncc, STATUS!= 99 & drug.x==0)))

Q <- rbind ( c(0, 0.25, 0,    0.25),
             c(0, 0,    0.25, 0.25),
             c(0, 0,    0,    0.25),
             c(0, 0,    0,    0) )

Q

Q.ini <- crudeinits.msm(STATUS ~ MONTH, newid, data=ncc, qmatrix=Q, censor=99, censor.states=c(2,3,4))
Q.ini

###########################
# hazard ratio estimation #
###########################

fit0 <- msm(STATUS ~ MONTH, subject=newid, data = ncc, censor=99, censor.states=c(2,3,4),
              qmatrix = Q.ini,  covariates = ~ drug.x)

hazard.msm(fit0)



# note, above fitting assume independence among patients with multiple
# locations. To correct that, we use bootstrap

# get bootstrap based hazard 95% CI

dobootstrap <- 0

if(dobootstrp == 1){
   B <- 1000 # the bootstrap size (WARNING: this takes 10 minutes on my PC)

   ids <- unique(ncc$ID) # we get bootstrap from re-sampling on the patients
   N   <- length(ids)

   set.seed(123)
   hm.boot <- NULL
   for(s in 1:B){
      bt.ids <- sample(ids, N, replace=T) # important to specify replace=T for bootstrap
      ncc.boot <- NULL

      for(i in 1:N){
         tt <- subset(ncc, ID %in% bt.ids[i])
         tt$id2 <- paste(tt$newid, as.character(i), sep="-") #need to define a new id for msm
         ncc.boot <- rbind(ncc.boot, tt)
      }
      fit <- msm(STATUS ~ MONTH, subject=id2, data = ncc.boot,
                censor=99, censor.states=c(2,3,4),
              qmatrix = Q.ini,  covariates = ~ drug.x)

      hm <- hazard.msm(fit)
      hm.boot <- rbind(hm.boot, hm$drug.x[, 1])
   }

   my.quantile <- function(x){
      quantile(x, probs = c(0.025, 0.975))
   }

   hm.median <- apply(hm.boot, 2, median);
   hm.CI    <- apply(hm.boot, 2, function(x) my.quantile(x))
   round(t(rbind(hm.median, hm.CI)),2)
}

#                   hm.median 2.5% 97.5%
# State 1 - State 2      2.09 0.68  7.66
# State 1 - State 4      3.00 1.43  7.23
# State 2 - State 3      1.05 0.00  6.14
# State 2 - State 4      1.26 0.66  2.53
# State 3 - State 4      1.08 0.45  2.23

#############################
# plot the survival curve   #
#############################

plotsurvival <- 0

if(plotsurvival == 1){

 png('transition_prob.png', width = 8, height = 8, units = 'in', res = 300)

 par(mfrow=c(2,2), oma = c(4, 3, 0, 0), mar = c(1,3,2,2), xpd = NA)


 plot(fit0, covariate = list(0), legend.pos=c(10,1), ylab="Transition probability",
         xlab="", las=1)
 text(2, 0.1, "Placebo")

 plot(c(0,0))

 plot(fit0, covariate = list(1), legend.pos=c(10,1), ylab="Transition probability",
      xlab="Time after starting of the trial (in month)", las=1)
 text(2, 0.1, "ALB")

 dev.off()

 plot.prevalence.msm(fit0)
 plot.survfit.msm(fit0, mark.time = FALSE)

 tt <- seq(0,24, by=0.5)
 len <- length(tt); len
 p12 <- numeric(len); p12.L <- numeric(len); p12.U <- numeric(len)
 p14 <- numeric(len); p14.L <- numeric(len); p14.U <- numeric(len)
 p23 <- numeric(len); p23.L <- numeric(len); p23.U <- numeric(len)
 p24 <- numeric(len); p24.L <- numeric(len); p24.U <- numeric(len)
 p34 <- numeric(len); p34.L <- numeric(len); p34.U <- numeric(len)

 for(i in 1:len){
    pp <- pmatrix.msm(fit0, t=tt[i], covariate = list(0), ci="normal")
    p12[i] <- pp$estimate[1,2]; p12.L[i] <- pp$L[1,2]; p12.U[i] <- pp$U[1,2]
    p14[i] <- pp$estimate[1,4]; p14.L[i] <- pp$L[1,4]; p14.U[i] <- pp$U[1,4]
    p23[i] <- pp$estimate[2,3]; p23.L[i] <- pp$L[2,3]; p23.U[i] <- pp$U[2,3]
    p24[i] <- pp$estimate[2,4]; p24.L[i] <- pp$L[2,4]; p24.U[i] <- pp$U[2,4]
    p34[i] <- pp$estimate[3,4]; p34.L[i] <- pp$L[3,4]; p34.U[i] <- pp$U[3,4]
 }

 plot(tt, p12, ylim=c(0,1), type="l")
 lines(tt, p14)
 lines(tt, p23)
 lines(tt, p24)
 lines(tt, p34)

 plot.prevalence.msm(fit0)

 my.pp <- function(pp){
 round(100*rbind(
   c(pp$estimates[1,2], pp$L[1,2], pp$U[1,2]),
   c(pp$estimates[1,4], pp$L[1,4], pp$U[1,4]),
   c(pp$estimates[2,3], pp$L[2,3], pp$U[2,3]),
   c(pp$estimates[2,4], pp$L[2,4], pp$U[2,4]),
   c(pp$estimates[3,4], pp$L[3,4], pp$U[3,4])),1)
 }

 pp1c <- pmatrix.msm(fit0, t=1, covariate = list(0), ci="normal")

 my.pp(pp1c)

 pp6c <- pmatrix.msm(fit0, t=6, covariate = list(0), ci="normal")

 my.pp(pp6c)

 pp1 <- pmatrix.msm(fit0, t=12, covariate = list(0), ci="normal")

 my.pp(pp1)

 pp2 <- pmatrix.msm(fit0, t=24, covariate = list(0), ci="normal")

 my.pp(pp2)

 #> my.pp(pp1)
 #     [,1] [,2] [,3]
 #[1,] 10.0  4.7 17.9
 #[2,] 68.2 56.2 80.5
 #[3,]  6.6  2.6 14.2
 #[4,] 75.6 63.8 85.0
 #[5,] 68.5 50.1 85.9
 #> pp2 <- pmatrix.msm(fit0, t=24, covariate = list(0), ci="normal")
 #> my.pp(pp2)
 #     [,1] [,2] [,3]
 #[1,]  3.8  1.5  7.5
 #[2,] 90.6 83.4 95.6
 #[3,]  3.3  1.1  7.9
 #[4,] 93.6 87.0 97.4
 #[5,] 90.1 74.4 97.8

 pp3 <- pmatrix.msm(fit0, t=12, covariate = list(1), ci="normal")

 my.pp(pp3)

 pp4 <- pmatrix.msm(fit0, t=24, covariate = list(1), ci="normal")

 my.pp(pp4)

 #> my.pp(pp3)
 #    [,1] [,2] [,3]
 #[1,]  5.3  2.6  9.3
 #[2,] 91.8 86.7 95.4
 #[3,]  5.4  2.0 12.5
 #[4,] 82.4 71.0 90.3
 #[5,] 71.3 56.1 85.3
 #> pp4 <- pmatrix.msm(fit0, t=24, covariate = list(1), ci="normal")
 #> my.pp(pp4)
 #     [,1] [,2] [,3]
 #[1,]  0.7  0.2  2.1
 #[2,] 98.5 96.7 99.4
 #[3,]  2.2  0.6  5.5
 #[4,] 96.3 91.3 98.6
 #[5,] 91.8 79.4 97.8
}

 #################################
 #                               #
 #  subgroup analysis            #
 #                               #
 #  assuming a group variable    #
 #  is coded 0/1                 #
 #                               #
 #  refer to tutorial for details#
 #################################

 ncc$loc <- ifelse(ncc$Tissue != 'Parenchymal' | is.na(ncc$Tissue), 0, 1)

 fit1 <- msm(STATUS ~ MONTH, subject=newid, data = ncc, censor=99, censor.states=c(2,3,4),
             qmatrix = Q.ini,  covariates = ~ drug.x * loc)

 res <- hazard.msm(fit1)
 is.list(res)
 names(res)

 beta.drug.x <- log(res$drug.x[,1])

 beta.loc <- log(res$loc[,1])

 beta.drug.x.loc <- log(res$`drug.x:loc`[,1])

 # hazard for drug.x = 1, and loc = 1
 hazardA <- exp(beta.drug.x + beta.loc + beta.drug.x.loc) #note I did not include the baseline beta since it will
                                                          #be cancelled anyway

 # hazard for drug.x = 0, and loc = 1
 hazardB <- exp(0 + beta.loc + 0)

 # hazard ratio (treatment vs. placebo ) for loc = 1
 HR.loc1 <- hazardA/hazardB

 # hazard for drug.x = 1, and loc = 0
 hazardC <- exp(beta.drug.x + 0 + 0)

 # hazard for drug.x = 0, and loc = 0
 hazardD <- exp(0 + 0 + 0)

 # hazard ratio (treatment vs. placebo ) for group = 1
 HR.loc0 <- hazardC/hazardD

 # put results together
 (interaction.loc <- t(rbind(HR.loc0, HR.loc1)))

 # LR test on interaction term

 fit1a <- msm(STATUS ~ MONTH, subject=newid, data = ncc, censor=99, censor.states=c(2,3,4),
             qmatrix = Q.ini,  covariates = ~ drug.x  + loc)

 lrtest.msm(fit1a, fit1)

 #> lrtest.msm(fit1a, fit1)
 #     -2 log LR df        p
 #fit1  8.253796  5 0.142792

 # readin the patient info data

 patinfo <- read_spss("patinfo.sav")
 names(patinfo)

 age.raw <- as.numeric((patinfo$PIADATE - patinfo$PIA3)/365)
 summary(age.raw)
# >  summary(age.raw)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3.142  27.529  40.295  40.582  50.958  82.068


 hist(age.raw)
 age <- ifelse(age.raw > 40, 1, 0)

 male <- ifelse(patinfo$PIA2==1, 1, 0);

 sub.patinfo <- data.frame(cbind(patinfo$ID, age, male))
 names(sub.patinfo) <- c("ID", "age", "male")


 ncc2 <- merge(ncc, sub.patinfo, by="ID")

 # descriptive statistics
 loc2 <- NULL;
 drug2 <- NULL
 ids2 <- unique(ncc2$newid)
 len2 <- length(ids2)
 for(i in 1:len2){
   tt <- subset(ncc2, newid %in% ids2[i])
   loc2 <- c(loc2, tt$loc[1])
   drug2 <- c(drug2, tt$drug.x[1])
 }
 table(loc2)
 summary(loc2)

# > table(loc2)
# loc2
# 0   1
# 89 132
# > summary(loc2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.0000  1.0000  0.5973  1.0000  1.0000
 table(drug2)
 summary(drug2)

 #>  table(drug2)
 #drug2
 #0   1
 #97 124
 #> summary(drug2)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 #0.0000  0.0000  1.0000  0.5611  1.0000  1.0000


 fit2 <- msm(STATUS ~ MONTH, subject=newid, data = ncc2, censor=99, censor.states=c(2,3,4),
             qmatrix = Q.ini,  covariates = ~ drug.x * age)

 res <- hazard.msm(fit2)

 beta.drug.x <- log(res$drug.x[,1])

 beta.age <- log(res$age[,1])

 beta.drug.x.age <- log(res$`drug.x:age`[,1])

 # hazard for drug.x = 1, and group = 1
 hazardA <- exp(beta.drug.x + beta.age + beta.drug.x.age)

 # hazard for drug.x = 0, and group = 1
 hazardB <- exp(0 + beta.age + 0)

 # hazard ratio (treatment vs. placebo ) for group = 1
 HR.age1 <- hazardA/hazardB

 # hazard for drug.x = 1, and group = 0
 hazardC <- exp(beta.drug.x + 0 + 0)

 # hazard for drug.x = 0, and group = 0
 hazardD <- exp(0 + 0 + 0)

 # hazard ratio (treatment vs. placebo ) for group = 1
 HR.age0 <- hazardC/hazardD

 # put results together
 (interaction.age <- t(rbind(HR.age0, HR.age1)))

 fit2a <- msm(STATUS ~ MONTH, subject=newid, data = ncc2, censor=99, censor.states=c(2,3,4),
             qmatrix = Q.ini,  covariates = ~ drug.x + age)
 lrtest.msm(fit2a, fit2)

# > lrtest.msm(fit2a, fit2)
#      -2 log LR df          p
# fit2  10.75516  5 0.05645648


 fit3 <- msm(STATUS ~ MONTH, subject=newid, data = ncc2, censor=99, censor.states=c(2,3,4),
             qmatrix = Q.ini,  covariates = ~ drug.x * male)

 res <- hazard.msm(fit3)

 beta.drug.x <- log(res$drug.x[,1])

 beta.male <- log(res$male[,1])

 beta.drug.x.male <- log(res$`drug.x:male`[,1])

 # hazard for drug.x = 1, and group = 1
 hazardA <- exp(beta.drug.x + beta.male + beta.drug.x.male)

 # hazard for drug.x = 0, and group = 1
 hazardB <- exp(0 + beta.age + 0)

 # hazard ratio (treatment vs. placebo ) for group = 1
 HR.male1 <- hazardA/hazardB

 # hazard for drug.x = 1, and group = 0
 hazardC <- exp(beta.drug.x + 0 + 0)

 # hazard for drug.x = 0, and group = 0
 hazardD <- exp(0 + 0 + 0)

 # hazard ratio (treatment vs. placebo ) for group = 1
 HR.male0 <- hazardC/hazardD

 # put results together
 (interaction.male <- t(rbind(HR.male0, HR.male1)))

 fit3a <- msm(STATUS ~ MONTH, subject=newid, data = ncc2, censor=99, censor.states=c(2,3,4),
             qmatrix = Q.ini,  covariates = ~ drug.x + male)
 lrtest.msm(fit3a, fit3)

# > lrtest.msm(fit3a, fit3)
#      -2 log LR df           p
# fit3   17.1567  5 0.004211998


 # read-in a data with ncyst, pparen, ptran info
 nccSeed123 <- read.csv("NccSeed123.csv", header=T)

 ids <- unique(ncc2$newid)
 len <- length(ids)

 ncc3  <- NULL
 for(i in 1:len){
   tt <- subset(ncc2, newid %in% ids[i])
   nv <- dim(tt)[1]
   ttt <- subset(nccSeed123, ID == tt$ID[1] & visit == tt$MONTH[1])
   cat("i=", i, "dim(ttt)=", dim(ttt), "\n");
   if(dim(ttt)[1] > 0){
      tt$ncyst <- rep(ttt$ncyst, nv) -1
      tt$pparen <- rep(ttt$pparen, nv)
      tt$ptrans <- rep(ttt$ptrans, nv)
      ncc3 <- rbind(ncc3, tt)
   }
 }
dim(ncc2);
dim(ncc3)


fit0 <- msm(STATUS ~ MONTH, subject=newid, data = ncc3, censor=99, censor.states=c(2,3,4),
            qmatrix = Q.ini,  covariates = ~ 1 )

fit1 <- msm(STATUS ~ MONTH, subject=newid, data = ncc3, censor=99, censor.states=c(2,3,4),
            qmatrix = Q.ini,  covariates = ~ drug.x )

fit2 <- msm(STATUS ~ MONTH, subject=newid, data = ncc3, censor=99, censor.states=c(2,3,4),
            qmatrix = Q.ini,  covariates = ~ loc)
hazard.msm(fit2)

fit3 <- msm(STATUS ~ MONTH, subject=newid, data = ncc3, censor=99, censor.states=c(2,3,4),
            qmatrix = Q.ini,  covariates = ~ ncyst)

fit4 <- msm(STATUS ~ MONTH, subject=newid, data = ncc3, censor=99, censor.states=c(2,3,4),
            qmatrix = Q.ini,  covariates = ~ pparen)

fit5 <- msm(STATUS ~ MONTH, subject=newid, data = ncc3, censor=99, censor.states=c(2,3,4),
            qmatrix = Q.ini,  covariates = ~ ptrans)

lrtest.msm(fit0, fit1)

lrtest.msm(fit0, fit2)

lrtest.msm(fit0, fit3)

lrtest.msm(fit0, fit4)

lrtest.msm(fit0, fit5)

fit1a <- msm(STATUS ~ MONTH, subject=newid, data = ncc3, censor=99, censor.states=c(2,3,4),
            qmatrix = Q.ini,  covariates = ~ drug.x + pparen + ptrans)
hazard.msm(fit1a)

fit1b <- msm(STATUS ~ MONTH, subject=newid, data = ncc3, censor=99, censor.states=c(2,3,4),
             qmatrix = Q.ini,  covariates = ~ drug.x*loc + pparen + ptrans)
hazard.msm(fit1b)

#> lrtest.msm(fit0, fit1)
#-2 log LR df            p
#fit1  22.71512  5 0.0003826277
#> lrtest.msm(fit0, fit2)
#-2 log LR df            p
#fit2  40.72464  5 1.066462e-07

hazard.msm(fit0)
hazard.msm(fit1)

res <- hazard.msm(fit1)
is.list(res)
names(res)

beta.drug.x <- log(res$drug.x[,1])

beta.loc <- log(res$loc[,1])

beta.drug.x.loc <- log(res$`drug.x:loc`[,1])

# hazard for drug.x = 1, and loc = 1
hazardA <- exp(beta.drug.x + beta.loc + beta.drug.x.loc) #note I did not include the baseline beta since it will
#be cancelled anyway

# hazard for drug.x = 0, and loc = 1
hazardB <- exp(0 + beta.loc + 0)

# hazard ratio (treatment vs. placebo ) for loc = 1
HR.loc1 <- hazardA/hazardB

# hazard for drug.x = 1, and loc = 0
hazardC <- exp(beta.drug.x + 0 + 0)

# hazard for drug.x = 0, and loc = 0
hazardD <- exp(0 + 0 + 0)

# hazard ratio (treatment vs. placebo ) for group = 1
HR.loc0 <- hazardC/hazardD

# put results together
(interaction.loc <- t(rbind(HR.loc0, HR.loc1)))

# LR test on interaction term

fit1a <- msm(STATUS ~ MONTH, subject=newid, data = ncc, censor=99, censor.states=c(2,3,4),
             qmatrix = Q.ini,  covariates = ~ drug.x  + loc)

lrtest.msm(fit1a, fit1)



 # bootstrap for interaction analysis

 # get bootstrap based hazard 95% CI

 dobootstrap <- 0

 if(dobootstrpa == 1){
    B <- 1000 # the bootstrap size (WARNING: this takes 10 minutes on my PC)

    ids <- unique(ncc2$ID) # we get bootstrap from re-sampling on the patients
    N   <- length(ids)

    set.seed(123)
    hm.boot0 <- NULL; hm.boot1 <- NULL
    hm.boot.age0 <- NULL; hm.boot.age1 <- NULL
    hm.boot.male0 <- NULL; hm.boot.male1 <- NULL

   for(s in 1:B){
      bt.ids <- sample(ids, N, replace=T) # important to specify replace=T for bootstrap
      ncc.boot <- NULL

      for(i in 1:N){
         tt <- subset(ncc2, ID %in% bt.ids[i])
         tt$id2 <- paste(tt$newid, as.character(i), sep="-") #need to define a new id for msm
         ncc.boot <- rbind(ncc.boot, tt)
      }

      fit <- msm(STATUS ~ MONTH, subject=id2, data = ncc.boot,
              censor=99, censor.states=c(2,3,4),
              qmatrix = Q.ini,  covariates = ~ drug.x * loc)

      res <- hazard.msm(fit)

      beta.drug.x <- log(res$drug.x[,1]); beta.loc <- log(res$loc[,1]); beta.drug.x.loc <- log(res$`drug.x:loc`[,1])

      hazardA <- exp(beta.drug.x + beta.loc + beta.drug.x.loc) #note I did not include the baseline beta since it will
      hazardB <- exp(0 + beta.loc + 0)

      HR.loc1 <- hazardA/hazardB

      hazardC <- exp(beta.drug.x + 0 + 0)
      hazardD <- exp(0 + 0 + 0)

      HR.loc0 <- hazardC/hazardD

      hm.boot0 <- rbind(hm.boot0, HR.loc0); hm.boot1 <- rbind(hm.boot1, HR.loc1)

      # drug * age

      fit <- msm(STATUS ~ MONTH, subject=id2, data = ncc.boot,
              censor=99, censor.states=c(2,3,4),
              qmatrix = Q.ini,  covariates = ~ drug.x * age)

      res <- hazard.msm(fit)

      beta.drug.x <- log(res$drug.x[,1]); beta.age <- log(res$age[,1]); beta.drug.x.age <- log(res$`drug.x:age`[,1])

      hazardA <- exp(beta.drug.x + beta.age + beta.drug.x.age) #note I did not include the baseline beta since it will
      hazardB <- exp(0 + beta.age + 0)

      HR.age1 <- hazardA/hazardB

      hazardC <- exp(beta.drug.x + 0 + 0)
      hazardD <- exp(0 + 0 + 0)

      HR.age0 <- hazardC/hazardD

      hm.boot.age0 <- rbind(hm.boot.age0, HR.age0); hm.boot.age1 <- rbind(hm.boot.age1, HR.age1)

      # drug * male

      fit <- msm(STATUS ~ MONTH, subject=id2, data = ncc.boot,
              censor=99, censor.states=c(2,3,4),
              qmatrix = Q.ini,  covariates = ~ drug.x * male)

      res <- hazard.msm(fit)

      beta.drug.x <- log(res$drug.x[,1]); beta.male <- log(res$male[,1]); beta.drug.x.male <- log(res$`drug.x:male`[,1])

      hazardA <- exp(beta.drug.x + beta.male + beta.drug.x.male) #note I did not include the baseline beta since it will
      hazardB <- exp(0 + beta.male + 0)

      HR.male1 <- hazardA/hazardB

      hazardC <- exp(beta.drug.x + 0 + 0)
      hazardD <- exp(0 + 0 + 0)

      HR.male0 <- hazardC/hazardD

      hm.boot.male0 <- rbind(hm.boot.male0, HR.male0); hm.boot.male1 <- rbind(hm.boot.male1, HR.male1)

    }

    my.quantile <- function(x){
       quantile(x, probs = c(0.025, 0.975))
    }

    # drug * loc

    HR.loc0 <- apply(hm.boot0, 2, median);
    hm.CI0    <- apply(hm.boot0, 2, function(x) my.quantile(x))


    HR.loc1 <- apply(hm.boot1, 2, median);
    hm.CI1    <- apply(hm.boot1, 2, function(x) my.quantile(x))

    druglocint <- round(cbind(t(rbind(HR.loc0, hm.CI0)), t(rbind(HR.loc1, hm.CI1))),2)

    # drug * age

    HR.age0 <- apply(hm.boot.age0, 2, median);
    hm.CI.age0    <- apply(hm.boot.age0, 2, function(x) my.quantile(x))

    HR.age1 <- apply(hm.boot.age1, 2, median);
    hm.CI.age1    <- apply(hm.boot.age1, 2, function(x) my.quantile(x))

    drugageint <- round(cbind(t(rbind(HR.age0, hm.CI.age0)), t(rbind(HR.age1, hm.CI.age1))),2)

    # drug * male

    HR.male0 <- apply(hm.boot.male0, 2, median);
    hm.CI.male0    <- apply(hm.boot.male0, 2, function(x) my.quantile(x))


    HR.male1 <- apply(hm.boot.male1, 2, median);
    hm.CI.male1    <- apply(hm.boot.male1, 2, function(x) my.quantile(x))

    drugmaleint <- round(cbind(t(rbind(HR.male0, hm.CI.male0)), t(rbind(HR.male1, hm.CI.male1))),2)

    druglocint

    drugageint

    drugmaleint

 }


#  >  druglocint
#                    HR.loc0 2.5%      97.5% HR.loc1 2.5%      97.5%
#  State 1 - State 2    0.36 0.00      16.31    4.44 0.97      >1000
#  State 1 - State 4    2.33 0.60      11.42    2.70 1.07      11.40
#  State 2 - State 3    0.87 0.00      >1000    0.94 0.00      10.10
#  State 2 - State 4    1.96 0.49      >1000    1.16 0.50       2.71
#  State 3 - State 4    0.67 0.00       3.08    1.10 0.32       3.28



#  >  drugageint
#                    HR.age0 2.5%      97.5%   HR.age1 2.5%        97.5%
#  State 1 - State 2    3.56 0.54      >1000      1.09 0.00         6.01
#  State 1 - State 4    3.80 1.30      13.04      2.36 0.96         7.12
#  State 2 - State 3    0.22 0.00      >1000      2.26 0.00        >1000
#  State 2 - State 4    0.68 0.17       1.99      2.61 0.74        11.31
#  State 3 - State 4    0.62 0.12       2.06      1.61 0.51         7.20

#   >  drugmaleint
#                    HR.male0 2.5%      97.5% HR.male1 2.5%       97.5%
#  State 1 - State 2     1.37 0.06      >1000     3.79 0.98       20.64
#  State 1 - State 4     1.73 0.25       5.91     4.84 1.81       15.70
#  State 2 - State 3     0.27 0.00      >1000     3.64 0.00       >1000
#  State 2 - State 4     0.65 0.18       1.77     3.05 1.04       14.54
#  State 3 - State 4     1.25 0.20      >1000     1.34 0.57        3.52

 #medhis  <- read_spss("medhis.sav")

 # adjusting ncyst, pparen, ptrans

 B <- 1000 # the bootstrap size

 ids <- unique(ncc3$ID) # we get bootstrap from re-sampling on the patients
 N   <- length(ids)

 set.seed(123)
 hm.boot0 <- NULL; hm.boot1 <- NULL
 hm.boot.age0 <- NULL; hm.boot.age1 <- NULL
 hm.boot.male0 <- NULL; hm.boot.male1 <- NULL

 for(s in 1:B){
   bt.ids <- sample(ids, N, replace=T) # important to specify replace=T for bootstrap
   ncc.boot <- NULL

   for(i in 1:N){
     tt <- subset(ncc3, ID %in% bt.ids[i])
     tt$id2 <- paste(tt$newid, as.character(i), sep="-") #need to define a new id for msm
     ncc.boot <- rbind(ncc.boot, tt)
   }

   fit <- msm(STATUS ~ MONTH, subject=id2, data = ncc.boot,
              censor=99, censor.states=c(2,3,4),
              qmatrix = Q.ini,  covariates = ~ drug.x * loc + ncyst + pparen + ptrans)

   res <- hazard.msm(fit)

   beta.drug.x <- log(res$drug.x[,1]); beta.loc <- log(res$loc[,1]); beta.drug.x.loc <- log(res$`drug.x:loc`[,1])

   hazardA <- exp(beta.drug.x + beta.loc + beta.drug.x.loc) #note I did not include the baseline beta since it will
   hazardB <- exp(0 + beta.loc + 0)

   HR.loc1 <- hazardA/hazardB

   hazardC <- exp(beta.drug.x + 0 + 0)
   hazardD <- exp(0 + 0 + 0)

   HR.loc0 <- hazardC/hazardD

   hm.boot0 <- rbind(hm.boot0, HR.loc0); hm.boot1 <- rbind(hm.boot1, HR.loc1)

   # drug * age

   fit <- msm(STATUS ~ MONTH, subject=id2, data = ncc.boot,
              censor=99, censor.states=c(2,3,4),
              qmatrix = Q.ini,  covariates = ~ drug.x * age+ ncyst + pparen + ptrans)

   res <- hazard.msm(fit)

   beta.drug.x <- log(res$drug.x[,1]); beta.age <- log(res$age[,1]); beta.drug.x.age <- log(res$`drug.x:age`[,1])

   hazardA <- exp(beta.drug.x + beta.age + beta.drug.x.age) #note I did not include the baseline beta since it will
   hazardB <- exp(0 + beta.age + 0)

   HR.age1 <- hazardA/hazardB

   hazardC <- exp(beta.drug.x + 0 + 0)
   hazardD <- exp(0 + 0 + 0)

   HR.age0 <- hazardC/hazardD

   hm.boot.age0 <- rbind(hm.boot.age0, HR.age0); hm.boot.age1 <- rbind(hm.boot.age1, HR.age1)

   # drug * male

   fit <- msm(STATUS ~ MONTH, subject=id2, data = ncc.boot,
              censor=99, censor.states=c(2,3,4),
              qmatrix = Q.ini,  covariates = ~ drug.x * male+ ncyst + pparen + ptrans)

   res <- hazard.msm(fit)

   beta.drug.x <- log(res$drug.x[,1]); beta.male <- log(res$male[,1]); beta.drug.x.male <- log(res$`drug.x:male`[,1])

   hazardA <- exp(beta.drug.x + beta.male + beta.drug.x.male) #note I did not include the baseline beta since it will
   hazardB <- exp(0 + beta.male + 0)

   HR.male1 <- hazardA/hazardB

   hazardC <- exp(beta.drug.x + 0 + 0)
   hazardD <- exp(0 + 0 + 0)

   HR.male0 <- hazardC/hazardD

   hm.boot.male0 <- rbind(hm.boot.male0, HR.male0); hm.boot.male1 <- rbind(hm.boot.male1, HR.male1)

 }

 my.quantile <- function(x){
   quantile(x, probs = c(0.025, 0.975))
 }

 # drug * loc

 HR.loc0 <- apply(hm.boot0, 2, median);
 hm.CI0    <- apply(hm.boot0, 2, function(x) my.quantile(x))


 HR.loc1 <- apply(hm.boot1, 2, median);
 hm.CI1    <- apply(hm.boot1, 2, function(x) my.quantile(x))

 druglocint <- round(cbind(t(rbind(HR.loc0, hm.CI0)), t(rbind(HR.loc1, hm.CI1))),2)

 # drug * age

 HR.age0 <- apply(hm.boot.age0, 2, median);
 hm.CI.age0    <- apply(hm.boot.age0, 2, function(x) my.quantile(x))


 HR.age1 <- apply(hm.boot.age1, 2, median);
 hm.CI.age1    <- apply(hm.boot.age1, 2, function(x) my.quantile(x))

 drugageint <- round(cbind(t(rbind(HR.age0, hm.CI.age0)), t(rbind(HR.age1, hm.CI.age1))),2)

 # drug * male

 HR.male0 <- apply(hm.boot.male0, 2, median);
 hm.CI.male0    <- apply(hm.boot.male0, 2, function(x) my.quantile(x))


 HR.male1 <- apply(hm.boot.male1, 2, median);
 hm.CI.male1    <- apply(hm.boot.male1, 2, function(x) my.quantile(x))

 drugmaleint <- round(cbind(t(rbind(HR.male0, hm.CI.male0)), t(rbind(HR.male1, hm.CI.male1))),2)

 druglocint

 drugageint

 drugmaleint

# >  druglocint
#                   HR.loc0 2.5%        97.5% HR.loc1 2.5%        97.5%
# State 1 - State 2    0.39 0.00         5.38    8.80 1.10 1.219805e+12
# State 1 - State 4    2.51 0.08        39.59    2.29 0.50 9.146158e+06
# State 2 - State 3    0.57 0.00 996680271.38    0.44 0.05 4.160000e+00
# State 2 - State 4    3.95 0.13   3931230.61    1.23 0.60 4.290000e+00
# State 3 - State 4    0.57 0.00        11.20    2.25 0.34 1.316000e+01
# >
# >  drugageint
#                   HR.age0 2.5%       97.5% HR.age1 2.5%        97.5%
# State 1 - State 2    2.88 0.17 36777205.35    1.70 0.01 3.370700e+02
# State 1 - State 4    3.42 0.81       44.19    2.42 0.27 2.317994e+08
# State 2 - State 3    0.22 0.00  9721949.58    0.24 0.00 2.401281e+13
# State 2 - State 4    0.59 0.03        2.70    4.81 0.81 2.897000e+01
# State 3 - State 4    1.19 0.11        9.62    2.50 0.22 1.654000e+01
# >
#   >  drugmaleint
#                   HR.male0 2.5%     97.5% HR.male1 2.5%    97.5%
# State 1 - State 2     0.55 0.00 736684.83     9.80 1.18  2719.66
# State 1 - State 4     1.28 0.02     41.94     6.20 0.26    43.45
# State 2 - State 3     0.28 0.00      8.31     1.88 0.00 72133.95
# State 2 - State 4     0.85 0.15 209070.47     7.86 0.15    43.99
# State 3 - State 4     1.83 0.11  68480.71     2.10 0.16    13.93
# >


 # adjusted analysis for overall effect with bootsrap

 dobootstrap <- 1

 if(dobootstrap == 1){
   B <- 1000 # the bootstrap size (WARNING: this takes 10 minutes on my PC)

   ids <- unique(ncc3$ID) # we get bootstrap from re-sampling on the patients
   N   <- length(ids)

   set.seed(123)
   hm.boot <- NULL
   for(s in 1:B){
     bt.ids <- sample(ids, N, replace=T) # important to specify replace=T for bootstrap
     ncc.boot <- NULL

     for(i in 1:N){
       tt <- subset(ncc3, ID %in% bt.ids[i])
       tt$id2 <- paste(tt$newid, as.character(i), sep="-") #need to define a new id for msm
       ncc.boot <- rbind(ncc.boot, tt)
     }
     fit <- msm(STATUS ~ MONTH, subject=id2, data = ncc.boot,
                censor=99, censor.states=c(2,3,4),
                qmatrix = Q.ini,  covariates = ~ drug.x + pparen + ptrans)

     hm <- hazard.msm(fit)
     hm.boot <- rbind(hm.boot, hm$drug.x[, 1])
   }

   my.quantile <- function(x){
     quantile(x, probs = c(0.025, 0.975))
   }

   hm.median <- apply(hm.boot, 2, median);
   hm.CI    <- apply(hm.boot, 2, function(x) my.quantile(x))
   round(t(rbind(hm.median, hm.CI)),2)
 }
#                   hm.median 2.5% 97.5%
# State 1 - State 2      2.06 0.50 11.15
# State 1 - State 4      3.17 1.50 12.94
# State 2 - State 3      0.81 0.00 11.57
# State 2 - State 4      1.53 0.62  6.32
# State 3 - State 4      1.35 0.38  5.53


