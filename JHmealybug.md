---
title: R analyses for Vea et al. Differential juvenile hormone modulation establishes
  extreme sexual dimorphism, for submission to PloS ONE in scale insects
author: "Isabelle Vea"
date: "September 14, 2015"
output:
  html_document:
    fig_caption: yes
---

#Introduction
This file presents the analyses performed in R to obtain the figures presented in Vea et al. Differential juvenile hormone modulation establishes extreme sexual dimorphism in scale insects.

```{r,echo=FALSE}
library(ggplot2)
library(plyr)
library(tidyr)

##script source for summarizeSE: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

```



#Loading datasets
```{r}
#for details of variable in each file, cf. readme.md
dataA<-read.csv(file="expressionprofile.csv",header = TRUE)
dataB<-read.csv(file="Pyri5mM.csv",header = TRUE)

#reshaping data
dataA2<-gather(dataA,Gene,SDM,6:19)
#summarizing by gene, day afther hatching and sex
dataA3<-ddply(dataA2,~Day.after.hatching +Sex+Gene,summarise,meanrp49=mean(SDM.rp49.2),meanSDM=mean(SDM))

#normalize dataB with housekeeping gene
dataB$gene.normal<-dataB$SDM.gene/dataB$rp49.1

#making mean and error values
dataBs<-summarySE(dataB, measurevar="gene.normal", groupvars=c("Treatment","Stage.treated","Gene"))

```

#Expression profiles
##Figure 2: Expression profiles of PkJHAMT, PkMet, PkTai, PkKr-h1-common
```{r}
Figure2<-subset(dataA3,Gene=="SDM.JHAMT"|Gene=="SDM.Met"|Gene=="SDM.Tai"|Gene=="SDM.Pkkr.h1_26")

pFig1<-ggplot(Figure2, aes(x=Day.after.hatching,y=meanSDM/meanrp49,group=Sex)) +
  geom_point(aes(linetype=Sex),size=2)+
  geom_line(aes(linetype=Sex),size=0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Gene relative expression (/rpL32)") +xlab("days after oviposition") +
  theme_bw(15) + 
  theme(axis.title.x = element_text(colour = "#242424"))

pFig1+facet_grid(Gene~.,scales="free")
```

##Figure 3: Effect on Met and Kr-h1
###Figure 3B: PkMet
```{r}
PkMet<-subset(dataBs,dataBs$Gene=="B.Met")

ggplot(PkMet, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("PkMet relative expression/rp49") +
    scale_fill_manual(name="", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on PkMet") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.1)) +
    theme_bw(15)
```


###Figure 3C: PkKr-h1
```{r}
PkKrh1<-subset(dataBs,dataBs$Gene=="D.kr-h1")

ggplot(PkKrh1, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("PkKr-h1 relative expression/rp49") +
    scale_fill_manual(name="Treatment", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on PkKr-h1") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.1)) +
    theme_bw(15)
```

###Figure 3D: PkKr-h1 A
```{r}
PkKrh1A<-subset(dataBs,(dataBs$Gene=="E.Pkkr-h1A"))

ggplot(PkKrh1A, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("PkKr-h1 A relative expression/rp49") +
    scale_fill_manual(name="Treatment", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on PkKr-h1 A") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.02)) +
    theme_bw(15)
```

###Figure 3E: PkKr-h1 B
```{r}
PkKrh1B<-subset(dataBs,(dataBs$Gene=="F.Pkkr-h1B"))

ggplot(PkKrh1B, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("PkKr-h1 B relative expression/rp49") +
    scale_fill_manual(name="Treatment", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on PkKr-h1 B") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.0005)) +
    theme_bw(15)

```

###Statistical tests
```{r}
#subsetting original data for statistical tests prepupae
prepupatestMet<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="B.Met")
t.test(gene.normal~Treatment,data=prepupatestMet)

prepupatestPkkrh1<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="D.kr-h1")
t.test(gene.normal~Treatment,data=prepupatestPkkrh1)

prepupatestPkkrh1A<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="E.Pkkr-h1A")
t.test(gene.normal~Treatment,data=prepupatestPkkrh1A)

prepupatestPkkrh1B<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="F.Pkkr-h1B")
t.test(gene.normal~Treatment,data=prepupatestPkkrh1B)

#subsetting original data for statistical tests pupae
pupatestMet<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="B.Met")
t.test(gene.normal~Treatment,data=pupatestMet)

pupatestPkkrh1<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="D.kr-h1")
t.test(gene.normal~Treatment,data=pupatestPkkrh1)

pupatestPkkrh1A<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="E.Pkkr-h1A")
t.test(gene.normal~Treatment,data=pupatestPkkrh1A)

pupatestPkkrh1B<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="F.Pkkr-h1B")
t.test(gene.normal~Treatment,data=pupatestPkkrh1B)
```

##Figure 5: Broad
###Figure 5A: Expression profile of Pkbr1 and Pkbr2
```{r}
Figure5<-subset(dataA3,Gene=="SDM.Pkbr1"|Gene=="SDM.Pkbr3")
pFig5<-ggplot(Figure5, aes(x=Day.after.hatching,y=meanSDM/meanrp49,group=Sex)) +
   geom_point(aes(linetype=Sex))+
   geom_line(aes(linetype=Sex))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Pkbr copies relative expression (/rpL32)") +xlab("days after oviposition") +
  theme_bw(15) + 
  theme(axis.title.x = element_text(colour = "#242424"))

pFig5+facet_grid(Gene~.,scale="free")

```

###Figure5B
```{r}
Pkbr1<-subset(dataBs,dataBs$Gene=="G.Pkbr1")

ggplot(Pkbr1, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("Pkbr1 relative expression/rp49") +
    scale_fill_manual(name="Treatment", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on Pkbr1") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.02)) +
    theme_bw(15)
```

###Figure5C
```{r}
Pkbr1Z2<-subset(dataBs,dataBs$Gene=="H.Pkbr1.Z2")

ggplot(Pkbr1Z2, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("Pkbr1 Z2 relative expression/rp49") +
    scale_fill_manual(name="Treatment", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on Pkbr1 Z2") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.01)) +
    theme_bw(15)
```

###Figure5D
```{r}
Pkbr1Z4<-subset(dataBs,dataBs$Gene=="I.Pkbr1.Z4")

ggplot(Pkbr1Z4, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("Pkbr1 Z4 relative expression/rp49") +
    scale_fill_manual(name="Treatment", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on Pkbr1 Z4") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.01)) +
    theme_bw(15)
```


###Figure5E
```{r}
Pkbr3<-subset(dataBs,dataBs$Gene=="M.Pkbr3")

ggplot(Pkbr3, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("Pkbr3 relative expression/rp49") +
    scale_fill_manual(name="Treatment", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on Pkbr3") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.1)) +
    theme_bw(15)
```
##Figure5F
```{r}
Pkbr3Z2<-subset(dataBs,dataBs$Gene=="N.Pkbr3.Z2")

ggplot(Pkbr3Z2, aes(x=Stage.treated, y=gene.normal, fill=Treatment)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=1) +      # Thinner lines
    geom_errorbar(aes(ymin=gene.normal-se, ymax=gene.normal+se),
                  size=.8,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("Pkbr3 Z2 relative expression/rp49") +
    scale_fill_manual(name="Treatment", # Legend label, use darker colors
                   breaks=c("A.Methanol", "B.5mM"),
                   labels=c("methanol","5mM pyriproxyfen"),
                   values=c("white", "#888888"))+
    ggtitle("Effect of pyriproxyfen on Pkbr3 Z2") + 
    scale_y_continuous(breaks=seq(0,2,by = 0.1)) +
    theme_bw(15)
```

###Statistical tests
```{r}
#subsetting original data for statistical tests prepupae
prepupatestPkbr1<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="G.Pkbr1")
t.test(gene.normal~Treatment,data=prepupatestPkbr1)

prepupatestPkbr1z2<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="H.Pkbr1.Z2")
t.test(gene.normal~Treatment,data=prepupatestPkbr1z2)

prepupatestPkbr1z4<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="I.Pkbr1.Z4")
t.test(gene.normal~Treatment,data=prepupatestPkbr1z4)

prepupatestPkbr2<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="J.Pkbr2")
t.test(gene.normal~Treatment,data=prepupatestPkbr2)

prepupatestPkbr2z2<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="K.Pkbr2.Z2")
t.test(gene.normal~Treatment,data=prepupatestPkbr2z2)

prepupatestPkbr2z4<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="L.Pkbr2.Z4")
t.test(gene.normal~Treatment,data=prepupatestPkbr2z4)

prepupatestPkbr3<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="M.Pkbr3")
t.test(gene.normal~Treatment,data=prepupatestPkbr3)

prepupatestPkbr3z2<-subset(dataB,dataB$Stage.treated=="PreD1" & dataB$Gene=="N.Pkbr3.Z2")
t.test(gene.normal~Treatment,data=prepupatestPkbr3z2)
```

```{r}
#subsetting original data for statistical tests pupae
pupatestPkbr1<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="G.Pkbr1")
t.test(gene.normal~Treatment,data=pupatestPkbr1)

pupatestPkbr1z2<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="H.Pkbr1.Z2")
t.test(gene.normal~Treatment,data=pupatestPkbr1z2)

pupatestPkbr1z4<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="I.Pkbr1.Z4")
t.test(gene.normal~Treatment,data=pupatestPkbr1z4)

pupatestPkbr2<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="J.Pkbr2")
t.test(gene.normal~Treatment,data=pupatestPkbr2)

pupatestPkbr2z2<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="K.Pkbr2.Z2")
t.test(gene.normal~Treatment,data=pupatestPkbr2z2)

pupatestPkbr2z4<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="L.Pkbr2.Z4")
t.test(gene.normal~Treatment,data=pupatestPkbr2z4)

pupatestPkbr3<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="M.Pkbr3")
t.test(gene.normal~Treatment,data=pupatestPkbr3)

pupatestPkbr3z2<-subset(dataB,dataB$Stage.treated=="PuD0" & dataB$Gene=="N.Pkbr3.Z2")
t.test(gene.normal~Treatment,data=pupatestPkbr3z2)
```


