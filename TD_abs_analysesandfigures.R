########Inbreeding depression and Va of LH and disease
library(dplyr)
library(brms)
library(broom)
library(tidyverse)
library(nadiv)
library(tidybayes)
library(data.table)

setwd("")

pd<-read.csv("phenotypicdata_females.csv")
gd<-read.csv("geneticdata_females.csv")

##model 1
priorP<-set_prior("normal(0,5)", class = "b") +
  set_prior("cauchy(0,3)", class = "sd")

fit_pcb1_hw<-brm(breed~agecat*(year+Headwidth+Weight)+dftd+
                   (1|year)+
                   (1|Microchip)+(1|month)+(1|TrapID),
                 prior=priorP,
                 data=pd,
                 control = list(adapt_delta = 0.99,max_treedepth=15),
                 chains=4,cores=4,
                 family=bernoulli(),
                 iter = 9000,
                 warmup= 2000)


plot(fit_pcb1_hw,N=5)

summary(fit_pcb1_hw)
conditional_effects(fit_pcb1_hw)

###with a smoothed effect of year 
fit_pcb1_hw2<-brm(breed~s(year,by=agecat,k=5)+agecat*(Headwidth+Weight)+dftd+
                  (1|year)+(1|Microchip)+(1|month)+(1|TrapID),
                  prior=priorP,
                  data=pd,
                  control = list(adapt_delta = 0.99,max_treedepth=15),
                  chains=4,cores=4,
                  family=bernoulli(),
                  iter = 9000,
                  warmup= 2000)

plot(fit_pcb1_hw2,N=5)
conditional_effects(fit_pcb1_hw2)

#######
##inbreeding depression model
str(gd)
fit_pcbG<-brm(breed~agecat*(Fhat3)+Headwidth+Weight+
                (1|year)+(1|Microchip)+(1|month)+(1|TrapID),
              prior=priorP,
              data=gd,
              control = list(adapt_delta = 0.99, max_treedepth = 15),
              chains=4,cores=4,
              family=bernoulli(),
              iter = 9000,
              warmup= 2000)

plot(fit_pcbG,N=5)
summary(fit_pcbG)
conditional_effects(fit_pcbG)
