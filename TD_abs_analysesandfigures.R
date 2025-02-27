########Inbreeding depression and Va of LH and disease
library(dplyr)
library(ggplot2)
library(wesanderson)
library(brms)
library(broom)
library(tidyverse)
library(nadiv)
library(tidybayes)
library(data.table)

setwd("")

pd<-read.csv("phenotypicdata_females.csv")##meta3_females2003
gd<-read.csv("geneticdata_females.csv")##meta3_females2003G

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

#############
#####plots
##Figure 1
########## Visualise raw data and check sample sizes across years
unique(pd$agecat)
breed_props<-pd %>%
  group_by(year, agecat,breed) %>%
  summarise(n = n()) %>%
  #mutate(se=sd(breed)/sqrt(n)) %>%
  mutate(freq = n / sum(n))
str(breed_props)
head(breed_props)

breed_props$se<-sqrt(breed_props$freq*(1-breed_props$freq)/breed_props$n)
  
bp_se<- pd %>%
  group_by(year,agecat) %>%
  mutate(sd=sd(breed),n=n(),se = sd/sqrt(n)) 

bp_se<-unique(bp_se[c("year","agecat","se")])

bp_se2 <- pd %>%
  group_by(year, agecat,breed) %>%
  summarise(n = n()) %>%
  #mutate(se=sd(breed)/sqrt(n)) %>%
  mutate(freq = n / sum(n))

breed_props<-left_join(bp_se2,bp_se,by=c("year","agecat"))

breed_props$N<-(1/breed_props$freq)*breed_props$n
breed_props[is.na(breed_props$se),]$se<-0

##very little data before 2003 for 1 year olds. Don't know why that might be...
breed_props<-subset(breed_props,breed_props$breed==1)
breed_props$lwr<-breed_props$freq-breed_props$se
breed_props$upr<-breed_props$freq+breed_props$se
breed_props$lwr<-ifelse(breed_props$lwr<=0,0,breed_props$lwr)
breed_props$upr<-ifelse(breed_props$upr>=1,1,breed_props$upr)
breed_props1year<-subset(breed_props,breed_props$agecat=="1year")
breed_propsadult<-subset(breed_props,breed_props$agecat=="Adult")

#Scale factor
sfactor <- max(breed_props1year$freq)/max(breed_props1year$N)
#Plot
pdf("Fig1A.pdf",width = 11)
ggplot(data = breed_props1year, aes(x = factor(year)))+
  geom_point(aes(y = freq),alpha=0.3)+
  geom_line(aes(y = freq), size=0.8,group = 1,alpha=0.75) +
  geom_errorbar(data = breed_props1year, 
                aes(x = factor(year), ymin = lwr,ymax = upr),
                color="darkgrey",width=.2,
                position=position_dodge(0.05),alpha=0.8)+
  geom_col(aes(y=N*sfactor),alpha=0.1)+
  scale_y_continuous(name='Probability of female breeding\n',
                     sec.axis = sec_axis(trans= ~. /sfactor, name = "Sample Size\n"))+
  xlab("\nYear")+
  theme_classic()+
  theme(text=element_text(size=16))
dev.off()

#Scale factor
sfactor_adult <- max(breed_propsadult$freq)/max(breed_propsadult$N)
pdf("Fig1B.pdf",width = 11)
ggplot(data = breed_propsadult, aes(x = factor(year)))+
  geom_point(aes(y = freq),alpha=0.3)+
  geom_line(aes(y = freq), size=0.8,group = 1,alpha=0.75) +
  geom_errorbar(data = breed_propsadult, 
                aes(x = factor(year), ymin = lwr,ymax = upr), 
                color="darkgrey",width=.2,
                position=position_dodge(0.05),alpha=0.8)+
  geom_col(aes(y=N*sfactor_adult),alpha=0.1)+
  scale_y_continuous(name='Probability of female breeding\n',
                     sec.axis = sec_axis(trans= ~. /sfactor_adult, name = "Sample Size\n"))+
  xlab("\nYear")+
  theme_classic()+
  theme(text=element_text(size=16))
dev.off()

##############
###Figure 2.
summary(fit_pcb1_hw2)
ceAY<-conditional_effects(fit_pcb1_hw2,effects="year:agecat")
str(ceAY)
ceAY<-as.data.frame(ceAY[[1]])

pdf("/Volumes/KS_work/ABS_TasDevs_MSinprep/ABS_Fig2.pdf",height=8,width=14)
ggplot(data=pd,aes(x=year,y=breed,colour=agecat))+
  geom_point(position = position_jitter(w = 0.3, h = 0))+
  geom_line(data=ceAY,aes(x=effect1__,y=estimate__,colour=effect2__),linewidth=1)+
  geom_line(data=ceAY,aes(x=effect1__,y=lower__,colour=effect2__),linetype="dashed")+
  geom_line(data=ceAY,aes(x=effect1__,y=upper__,colour=effect2__),linetype="dashed")+
  facet_grid(.~agecat)+
  theme_classic()+
  theme(text = element_text(size=16))+
  theme(text=element_text(size=14),
        legend.text.align = 0,
        legend.title.align = 0.5,
        legend.position = c(0.93,0.22),
        legend.justification="top",
        legend.box.background=element_rect(),
        legend.box.margin=margin(5,5,5,5))+
  theme(panel.margin=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1),
        strip.text.x = element_blank())+
  scale_colour_manual("Age",values=c("goldenrod1","darkgreen"),labels=c("1-year-olds","Adults"))+
  xlab("\nYear")+
  ylab("Annual breeding success\n")
dev.off()  

require(modelr)

##plot effect -- Weight
pdf("/Volumes/KS_work/ABS_TasDevs_MSinprep/ABS_Fig1a.pdf")
pd %>%
  #data_grid(Weight = seq_range(Weight, n = 101)) %>%
  data_grid(Weight=seq_range(Weight , n = 101), 
            agecat = agecat, 
            #year=floor(mean(year)),
            Headwidth=floor(mean(Headwidth)),
            dftd=0.5,
            month=sample(month,20),
            Microchip=sample(Microchip,20),
            TrapID=sample(TrapID,20),
            year=factor(year))%>%
  # NOTE: this shows the use of ndraws to subsample within add_epred_draws()
  # ONLY do this IF you are planning to make spaghetti plots, etc.
  # NEVER subsample to a small sample to plot intervals, densities, etc.
  add_epred_draws(fit_pcb1_hw,ndraws = 100,allow_new_levels=TRUE) %>%
  ggplot(aes(x = Weight, y = breed)) +
  geom_line(aes(y = .epred, group = .draw),stat="smooth",method="glm",method.args=list(family="binomial"),
            se=FALSE,alpha = 1/6, color = "#808080") +
  geom_line(aes(y = .epred),stat="smooth",method="glm",method.args=list(family="binomial"),se=FALSE,color = "#808080") +
  geom_point(data = meta3_females2003,color = "#808080",alpha = 1/3) +
  xlab("\nWeight (kg)")+
  ylab("Breeding success\n")+
  theme_classic()+
  theme(text = element_text(size=14))
dev.off()

pdf("/Volumes/KS_work/ABS_TasDevs_MSinprep/ABS_Fig1b.pdf")
meta3_females2003 %>%
  #data_grid(Weight = seq_range(Weight, n = 101)) %>%
  data_grid(Headwidth=seq_range(Headwidth , n = 101), 
            agecat = agecat, 
            year=floor(mean(year)),
            Weight=floor(mean(Weight)),
            dftd=0.5,
            Fmonth=sample(Fmonth,20),
            Microchip=sample(Microchip,20),
            TrapID=sample(TrapID,20),
            Fyear=factor(year))%>%
  # NOTE: this shows the use of ndraws to subsample within add_epred_draws()
  # ONLY do this IF you are planning to make spaghetti plots, etc.
  # NEVER subsample to a small sample to plot intervals, densities, etc.
  add_epred_draws(fit_pcb1_hw,ndraws = 100,allow_new_levels=TRUE) %>%
  ggplot(aes(x = Headwidth, y = breed)) +
  geom_line(aes(y = .epred, group = .draw),stat="smooth",method="glm",method.args=list(family="binomial"),
            se=FALSE,alpha = 1/6, color = "#808080") +
  geom_line(aes(y = .epred),stat="smooth",method="glm",method.args=list(family="binomial"),se=FALSE,color = "#808080") +
  geom_point(data = meta3_females2003,color = "#808080",alpha = 1/3) +
  xlab("\nHeadwidth (mm)")+
  ylab("Breeding success\n")+
  theme_classic()+
  theme(text = element_text(size=14))
dev.off()

##plot for adults and 1 year olds separately
ceAH<-conditional_effects(fit_pcb1_hw,effects="Headwidth:agecat")
str(ceAH)
ceAH<-as.data.frame(ceAH[[1]])

pdf("/Volumes/KS_work/ABS_TasDevs_MSinprep/ABS_FigS1.pdf",height=6,width=7)
ggplot(data=pd,aes(x=Headwidth,y=breed,colour=agecat))+
  geom_point(position = position_jitter(w = 0.3, h = 0))+
  geom_line(data=ceAH,aes(x=effect1__,y=estimate__,colour=effect2__),linewidth=1)+
  geom_line(data=ceAH,aes(x=effect1__,y=lower__,colour=effect2__),linetype="dashed")+
  geom_line(data=ceAH,aes(x=effect1__,y=upper__,colour=effect2__),linetype="dashed")+
  theme_classic()+
  theme(text = element_text(size=16))+
  theme(text=element_text(size=14),
        legend.text.align = 0,
        legend.title.align = 0.5,
        #legend.position = c(0.93,0.22),
        legend.justification="top",
        legend.box.background=element_rect(),
        legend.box.margin=margin(5,5,5,5))+
  scale_colour_manual("Age",values=c("goldenrod1","darkgreen"),labels=c("1-year-olds","Adults"))+
  xlab("\nHeadwidth")+
  ylab("Annual breeding success\n")
dev.off()  

ceAW<-conditional_effects(fit_pcb1_hw,effects="Weight:agecat")
str(ceAW)
ceAW<-as.data.frame(ceAW[[1]])

pdf("/Volumes/KS_work/ABS_TasDevs_MSinprep/ABS_FigS2.pdf",height=6,width=7)
ggplot(data=pd,aes(x=Weight,y=breed,colour=agecat))+
  geom_point(position = position_jitter(w = 0.3, h = 0))+
  geom_line(data=ceAW,aes(x=effect1__,y=estimate__,colour=effect2__),linewidth=1)+
  geom_line(data=ceAW,aes(x=effect1__,y=lower__,colour=effect2__),linetype="dashed")+
  geom_line(data=ceAW,aes(x=effect1__,y=upper__,colour=effect2__),linetype="dashed")+
  theme_classic()+
  theme(text = element_text(size=16))+
  theme(text=element_text(size=14),
        legend.text.align = 0,
        legend.title.align = 0.5,
        #legend.position = c(0.93,0.22),
        legend.justification="top",
        legend.box.background=element_rect(),
        legend.box.margin=margin(5,5,5,5))+
  scale_colour_manual("Age",values=c("goldenrod1","darkgreen"),labels=c("1-year-olds","Adults"))+
  xlab("\nWeight")+
  ylab("Annual breeding success\n")
dev.off()  

