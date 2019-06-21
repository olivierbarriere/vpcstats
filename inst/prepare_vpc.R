#library(plyr)
#detach("package:plyr", unload = TRUE)
library(rlang)
library(vpc)
library(magrittr)
library(dplyr)
library(data.table)
library(vpcstats)
library(ggplot2)
library(ggrepel)
library(colourpicker)

#devtools::install_github("olivierbarriere/vpcstats")
#detach("package:vpcstats", unload = TRUE)
#source("vpcstats.r") #use locally edited version for debugging

covdata$Observed_Time_hours
covdata <- read.csv("obsdata.csv",skip = 0,header = TRUE)
covdata <- covdata[!is.na(covdata$exosome_ug_per_uL),]
covdata$IVAR <- covdata$Observed_Time_hours

exampleobsdata <- read.csv("Residuals.csv")
exampleobsdata$REPLICATE <- 0

exampleobsdata<- left_join(exampleobsdata,
                           covdata[,c("Animal_ID","Cell_line","IVAR")])

exampleobsdata$ID <- exampleobsdata$Animal_ID
exampleobsdata$TIME <- exampleobsdata$IVAR

ggquickeda::run_ggquickeda(covdata)


examplesimdata <- read.csv("PredCheckAll.csv")
examplesimdata$REPLICATE <- examplesimdata$Replicate
examplesimdata$Cell_line <- rep(exampleobsdata$Cell_line,1000)
examplesimdata$ID <- examplesimdata$Animal_ID
exampleobsdata$LLOQ <- -Inf
examplesimdata$LLOQ <- -Inf
examplesimdata$TAD <- rep(exampleobsdata$TAD,1000)
examplesimdata$TIME <- examplesimdata$IVAR


source("plotfunction.R")
source("vpc_app.r")
vpc_app(obsdata = exampleobsdata,
        simdata = examplesimdata,
        stratalist=c("","Cell_line") ) 


vpcstats(obsdata = exampleobsdata,
         simdata = examplesimdata,
         REPL = REPLICATE, NBINS = 7)

#filter qname zoom other enhancements