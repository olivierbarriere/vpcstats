source("R/Debug/computePI.bins.R")
source("R/Debug/compare_outputs.R")
source("R/computePI.R")
source("R/quantile_cens.R")
require(ggbeeswarm)
require(ggplot2)
require(vpc)
require(data.table)

theme_set(theme_bw())

# we go through examples to showcase use case scenarios and functionality
# at the same time a kind of a validation against vpc package
#and where we make it more general
exampleobs <- simple_data$obs
exampleobs <- exampleobs[exampleobs$MDV == 0, ]
examplesim <- simple_data$sim

A <- vpc(sim = examplesim, obs = exampleobs, pred_corr = FALSE, stratify = c("ISM"), bins = "none", lloq = 50)
B <- vpc_cens(sim = simple_data$sim, stratify = c("ISM"), obs = simple_data$obs, lloq = 50)
egg::ggarrange(A, B)
# nice plot including lloq

# current function is software agnostic
# so some granted nonmem outputs are not always there
# e.g. simulation table might not contain PRED

## for the computepi when pred corrected it needs PRED also it needs REP
## LLOQ need to be a column when requested
## LLOQ can vary by observation we will see later in the examples

#using the same data as above jsut making it compatible with the new function
exampleobs$PRED <- examplesim[1:nrow(exampleobs), "PRED"]
exampleobs$REP <- 0
exampleobs$LLOQ <- 50
examplesim <- examplesim[examplesim$MDV == 0, ]#
examplesim$LLOQ <- 50


# small benchmark removing plotting from vpc command
system.time(
  vpc(sim = examplesim, obs = exampleobs,
      pred_corr = FALSE, stratify = c("ISM"), bins = "none", lloq = 50,vpcdb=TRUE)
)


system.time(
  {
    compute.PI.bins(
      data = exampleobs, simdata = examplesim, stratify1 = ISM,
      nbins = NULL, LLOQ = LLOQ, logadditive = FALSE, predcorrection = FALSE
    )
  }
)

system.time(
  {
    compute.PI(
      obsdata = exampleobs, simdata = examplesim, stratify = ~ISM,
      NBINS = NULL, LLOQ = LLOQ, predcorrection_islogdv = FALSE, predcorrection = FALSE
    )
  }
)


# almost the same for vpc and the old compute.PI.bins (~1s), but almost twice as fast for the new computePI (~0.5s)


# first example no strata with lloqoof 10
TESTA0 <- 
  compute.PI.bins(
    data = exampleobs, simdata = examplesim, stratify1 = NULL,
    nbins = NULL, LLOQ = LLOQ, logadditive = FALSE, predcorrection = FALSE
  )

TESTA <- 
  compute.PI(
    obsdata = exampleobs, simdata = examplesim, stratify = NULL,
    NBINS = NULL, LLOQ = LLOQ, predcorrection_islogdv = FALSE, predcorrection = FALSE
  )

compare_outputs(TESTA0, TESTA)
#Same output (beside renaming)

A <- vpc(sim = examplesim, obs = exampleobs,
         pred_corr = FALSE, bins = "none", lloq = 50)

AAA <- ggplot(TESTA) +
  geom_ribbon(data = TESTA, aes(XMED,
                                ymin = (`DV2.5%CI`), ymax = (`DV97.5%CI`),
                                fill = DVQNAME, col = DVQNAME,
                                group = DVQNAME
  ), alpha = 0.1, col = NA) +
  geom_line(data = TESTA, aes(XMED,
                              y = `DV50%CI`,
                              col = DVQNAME, group = DVQNAME
  )) +
  geom_line(data = TESTA, aes(
    x = XMED, y = DVOBS, group = DVQNAME,
    linetype = DVQNAME
  ), size = 2) +
  geom_hline(yintercept = 50, col = "red") +
  geom_rug(data = TESTA, aes(x = XMIN), sides = "t") +
  geom_rug(data = TESTA, aes(x = XMAX), sides = "t") +
  scale_colour_manual(
    name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
    breaks = c("5%PI", "50%PI", "95%PI", "Percent BLQ"),
    values = c("red", "blue", "red", "black"),
    labels = c("5%", "50%", "95%", "Percent BLQ")
  ) +
  scale_fill_manual(
    name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
    breaks = c("5%PI", "50%PI", "95%PI", "Percent BLQ"),
    values = c("red", "blue", "red", "black"),
    labels = c("5%", "50%", "95%", "Percent BLQ")
  ) +
  scale_linetype_manual(
    name = "Observed Percentiles\n(black lines)",
    breaks = c("5%PI", "50%PI", "95%PI"),
    values = c("dotted", "solid", "dashed"),
    labels = c("5%", "50%", "95%")
  ) +
  guides(
    fill = guide_legend(order = 2),
    colour = guide_legend(order = 2),
    linetype = guide_legend(order = 1)
  ) +
  theme(
    legend.position = "top", legend.key.width = grid::unit(2, "cm"),
    axis.text.x = element_text(angle = 30), axis.title.x = element_blank()
  )
egg::ggarrange(AAA, A)

# for now keepign raw ggplot as it enables use to do the customization you want

###########################
#############

# second strata with lloq of 50 and STRATA by is male

exampleobs$LLOQ <- 50
examplesim$LLOQ <- 50

TESTA0 <- compute.PI.bins(
  data = exampleobs, simdata = examplesim, stratify1 = ISM,
  nbins = NULL, LLOQ = LLOQ
)

TESTA <- compute.PI(
  obsdata = exampleobs, simdata = examplesim, stratify = ~ISM,
  NBINS = NULL, LLOQ = LLOQ
)

compare_outputs(TESTA0, TESTA)

A <- vpc(
  sim = examplesim, obs = exampleobs, pred_corr = FALSE, stratify = c("ISM"), # n_bins=5
  lloq = 50
)

AAA <- ggplot(TESTA) +
  facet_wrap(~ISM, scales = "free_x", labeller = label_wrap_gen(multi_line = FALSE), ncol = 2) +
  geom_ribbon(data = TESTA, aes(XMED,
                                ymin = (`DV2.5%CI`), ymax = (`DV97.5%CI`),
                                fill = DVQNAME, col = DVQNAME,
                                group = DVQNAME
  ), alpha = 0.1, col = NA) +
  geom_line(data = TESTA, aes(XMED,
                              y = `DV50%CI`,
                              col = DVQNAME, group = DVQNAME
  )) +
  geom_line(data = TESTA, aes(
    x = XMED, y = DVOBS, group = DVQNAME,
    linetype = DVQNAME
  ), size = 2) +
  geom_hline(yintercept = 50, col = "red") +
  geom_rug(data = TESTA, aes(x = XMIN), sides = "t") +
  geom_rug(data = TESTA, aes(x = XMAX), sides = "t") +
  scale_colour_manual(
    name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
    breaks = c("5%PI", "50%PI", "95%PI", "Percent BLQ"),
    values = c("red", "blue", "red", "black"),
    labels = c("5%", "50%", "95%", "Percent BLQ")
  ) +
  scale_fill_manual(
    name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
    breaks = c("5%PI", "50%PI", "95%PI", "Percent BLQ"),
    values = c("red", "blue", "red", "black"),
    labels = c("5%", "50%", "95%", "Percent BLQ")
  ) +
  scale_linetype_manual(
    name = "Observed Percentiles\n(black lines)",
    breaks = c("5%PI", "50%PI", "95%PI"),
    values = c("dotted", "solid", "dashed"),
    labels = c("5%", "50%", "95%")
  ) +
  guides(
    fill = guide_legend(order = 2),
    colour = guide_legend(order = 2),
    linetype = guide_legend(order = 1)
  ) +
  theme(
    legend.position = "top", legend.key.width = grid::unit(2, "cm"),
    axis.text.x = element_text(angle = 30), axis.title.x = element_blank()
  ) +
  geom_hline(yintercept = 50)
egg::ggarrange(AAA, A)
#


BBB <- ggplot(TESTA[TESTA$DVQNAME == TESTA$DVQNAME[1], ]) +
  facet_wrap(~ISM, scales = "free_x", labeller = label_wrap_gen(multi_line = FALSE), ncol = 2) +
  geom_ribbon(data = TESTA, aes(XMED,
                                ymin = (`PCTBLQ2.5%CI`), ymax = (`PCTBLQ97.5%CI`),
                                fill = PCTBLQQNAME, col = PCTBLQQNAME,
                                group = PCTBLQQNAME
  ), alpha = 0.1, col = NA) +
  geom_line(data = TESTA, aes(XMED,
                              y = `PCTBLQ50%CI`,
                              col = PCTBLQQNAME, group = PCTBLQQNAME
  )) +
  geom_line(data = TESTA, aes(
    x = XMED, y = PCTBLQOBS, group = PCTBLQQNAME,
    linetype = PCTBLQQNAME
  ), size = 2) +
  geom_rug(data = TESTA, aes(x = XMIN), inherit.aes = FALSE, sides = "t") +
  geom_rug(data = TESTA, aes(x = XMAX), inherit.aes = FALSE, sides = "t") +
  scale_colour_manual(
    name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
    breaks = c("5%PI", "50%PI", "95%PI", "PercentBLQ"),
    values = c("red", "blue", "red", "black"),
    labels = c("5%", "50%", "95%", "Percent BLQ")
  ) +
  scale_fill_manual(
    name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
    breaks = c("5%PI", "50%PI", "95%PI", "PercentBLQ"),
    values = c("red", "blue", "red", "black"),
    labels = c("5%", "50%", "95%", "Percent BLQ")
  ) +
  scale_linetype_manual(
    name = "Observed Percentiles\n(black lines)",
    breaks = c("5%PI", "50%PI", "95%PI"),
    values = c("dotted", "solid", "dashed"),
    labels = c("5%", "50%", "95%")
  ) +
  guides(
    fill = guide_legend(order = 2),
    colour = guide_legend(order = 2),
    linetype = guide_legend(order = 1)
  ) +
  theme(
    legend.position = "bottom", legend.key.width = grid::unit(2, "cm"),
    axis.text.x = element_text(angle = 30), axis.title.x = element_blank()
  )

B <- vpc_cens(sim = simple_data$sim, stratify = c("ISM"), obs = simple_data$obs, lloq = 50)

egg::ggarrange(BBB, B)
egg::ggarrange(AAA+
                 geom_hline(yintercept=c(25,100)),
               BBB)


### now examples with diferences and things that vpc cannot handle

exampleobs$LLOQ <- ifelse(exampleobs$ISM == 0, 100, 25)
examplesim$LLOQ <- ifelse(examplesim$ISM == 0, 100, 25)

TESTA0 <- compute.PI.bins(
  data = exampleobs, simdata = examplesim, stratify1 = ISM,
  nbins = NULL, LLOQ = LLOQ
)

TESTA <- compute.PI(
  obsdata = exampleobs, simdata = examplesim, stratify = ~ISM,
  NBINS = NULL, LLOQ = LLOQ
)

compare_outputs(TESTA0, TESTA)


A <- vpc(
  sim = examplesim, obs = exampleobs, pred_corr = FALSE, stratify = c("ISM"), # n_bins=5
  lloq = LLOQ
) # does not work with a column name

A <- vpc(
  sim = examplesim, obs = exampleobs, pred_corr = FALSE, stratify = c("ISM"), # n_bins=5
  lloq = 100
) # does with a scalar but 100 making something not work

A <- vpc(
  sim = examplesim, obs = exampleobs, pred_corr = FALSE, stratify = c("ISM"), # n_bins=5
  lloq = 50
)



# rerun AAA and BBB as is
egg::ggarrange(AAA, A)
egg::ggarrange(BBB, B)
egg::ggarrange(AAA, BBB)

# now with more than one strata

IDCOV <- sample(unique(exampleobs$ID), 25, replace = FALSE)
exampleobs$ISF <- ifelse(exampleobs$ID %in% IDCOV, 1, 0)
examplesim$ISF <- rep(exampleobs$ISF, 100)

TESTA0 <- compute.PI.bins(
  data = exampleobs, simdata = examplesim,
  nbins = NULL, LLOQ = LLOQ, stratify1 = ISM, stratify2 = ISF
)

TESTA <- compute.PI(
  obsdata = exampleobs, simdata = examplesim, stratify = ~ISM+ISF,
  NBINS = NULL, LLOQ = LLOQ
)

compare_outputs(TESTA0, TESTA)

A <- vpc(
  sim = examplesim, obs = exampleobs, pred_corr = FALSE, stratify = c("ISM", "ISF"), # n_bins=5
  lloq = 50
) # does not work


ggplot(TESTA) +
  facet_wrap(ISF ~ ISM, scales = "free_x", labeller = label_wrap_gen(multi_line = FALSE), ncol = 2) +
  geom_ribbon(data = TESTA, aes(XMED,
                                ymin = (`DV2.5%CI`), ymax = (`DV97.5%CI`),
                                fill = DVQNAME, col = DVQNAME,
                                group = DVQNAME
  ), alpha = 0.1, col = NA) +
  geom_line(data = TESTA, aes(XMED,
                              y = `DV50%CI`,
                              col = DVQNAME, group = DVQNAME
  )) +
  geom_line(data = TESTA, aes(
    x = XMED, y = DVOBS, group = DVQNAME,
    linetype = DVQNAME
  ), inherit.aes = FALSE, size = 2)


# Binning strategy global (vpc) or by group (computePI)
# truncate data by ISM
exampleobs1=simple_data$obs %>% filter(MDV==0 & ((ISM==0 & TIME>=2) | (ISM==1 & TIME<=6)))
examplesim1=simple_data$sim %>% filter(MDV==0 & ((ISM==0 & TIME>=2) | (ISM==1 & TIME<=6)))
n_bins = 5
A <- vpc(sim = examplesim1, obs = exampleobs1, stratify="ISM", bins="jenks", n_bins=n_bins) 
AA <- vpc(sim = examplesim1, obs = exampleobs1, stratify="ISM", bins="jenks", n_bins=n_bins, vpcdb=T) 
AA$vpc_dat
AA$aggr_obs
vpc:::auto_bin(exampleobs1 %>% mutate(idv=TIME), type="jenks", n_bins=n_bins)
vpc::bin_data
#uses cut(right=F), so I updated computePI to use the same.

#n global bins
B <- compute.PI(obsdata=exampleobs1, simdata = examplesim1, stratify=~ISM, NBINS=n_bins, style="jenks", bin_by_strata=F)
pvpc <- function(data) {
  ggplot(data, aes(x=XMID))+
  geom_ribbon(aes(ymin=`DV2.5%CI`, ymax=`DV97.5%CI`, fill=DVQNAME))+
  geom_line(aes(y=DVOBS, linetype=DVQNAME))+
  geom_rug(aes(x=XLEFT), sides="bt")+
  geom_rug(aes(x=XRIGHT), sides="bt")+
  facet_wrap(~ISM)
}
egg::ggarrange(A, pvpc(B))  
#Same

#nbins per strata
C <- compute.PI(obsdata=exampleobs1, simdata = examplesim1, stratify=~ISM, NBINS=n_bins, style="jenks", bin_by_strata=T)
egg::ggarrange(A, pvpc(B), pvpc(C))  


# to do list make the plotting a function that retrun ggplot object ready to be ggrranged

# ASCPT challenge parent metabolite pd

obsdata <- read.csv("data/obsdata.csv")
simdata <- read.csv("data/simdata.csv")


obsdata <- obsdata %>%
  arrange(DVTYPE, ID, TIME)
simdata <- simdata %>%
  arrange(REPL, DVTYPE, ID, TIME)

head(simdata)
head(obsdata)
obsdata <- obsdata[obsdata$TIME > 120, ]
simdata <- simdata[simdata$REPL < 20, ]
simdata <- simdata[simdata$TIME > 120, ]

# simdata<- simdata[simdata$DVTYPE=="EObs",]
# obsdata<- obsdata[obsdata$DVTYPE=="EObs",]

# source("computePI.method.strat.R")

PKPDVPC0 <- compute.PI.bins(
  data = obsdata, simdata = simdata, TIME = TIME, REPL = REPL,
  nbins = NULL, LLOQ = NULL, stratify1 = DVTYPE, stratify2 = NULL, predcorrection = FALSE
)

PKPDVPC <- compute.PI(
  obsdata = obsdata, simdata = simdata, TIME = TIME, REPL = REPL,
  NBINS = NULL, LLOQ = NULL, stratify = ~DVTYPE, predcorrection = FALSE
)

compare_outputs(PKPDVPC0, PKPDVPC)


AAA <- vpc(
  sim = simdata, obs = obsdata,
  pred_corr = FALSE, stratify = c("DVTYPE"), bins = "none", facet = "wrap"
)

AAA <- AAA + facet_wrap(~DVTYPE, scales = "free_y")


BBB <- ggplot(PKPDVPC) +
  facet_wrap(~DVTYPE, labeller = label_wrap_gen(multi_line = FALSE), ncol = 3, scales = "free_y") +
  geom_ribbon(data = PKPDVPC, aes(XMID,
                                  ymin = (`DV2.5%CI`), ymax = (`DV97.5%CI`),
                                  fill = DVQNAME, col = DVQNAME,
                                  group = DVQNAME
  ), alpha = 0.1, col = NA) +
  geom_line(data = PKPDVPC, aes(XMID,
                                y = `DV50%CI`,
                                col = DVQNAME, group = DVQNAME
  )) +
  geom_line(data = PKPDVPC, aes(
    x = XMID, y = DVOBS, group = DVQNAME,
    linetype = DVQNAME
  ), inherit.aes = FALSE, size = 2)
egg::ggarrange(AAA, BBB)

# vpc package is wrong did not investigate why and how nor did I submit github issue


# report ready VPC PLOT

ggplot(PKPDVPC) +
  geom_ribbon(aes(XMID,
                  ymin = (`DV2.5%CI`), ymax = (`DV97.5%CI`),
                  fill = DVQNAME, col = DVQNAME,
                  group = DVQNAME
  ), alpha = 0.1, col = NA) +
  facet_wrap(~DVTYPE, labeller = label_wrap_gen(multi_line = FALSE), ncol = 3, scales = "free_y") +
  geom_line(aes(XMID, y = `DV50%CI`, col = DVQNAME, group = DVQNAME)) +
  geom_line(aes(
    x = XMID, y = DVOBS, group = DVQNAME,
    linetype = "Observed"
  ), inherit.aes = FALSE, size = 1.5) +
  scale_linetype_manual(name = "", breaks = c("Observed"), values = c("dashed")) +
  scale_fill_manual(
    name = "Prediction Intervals 95% CI (ribbons)\nMedian Predictions (solid lines)",
    breaks = c("5%PI", "50%PI", "95%PI", "Obs"),
    labels = c("5%", "50%", "95%", "Obs"),
    values = c("#2c7fb8", "#df65b0", "#2c7fb8")
  ) +
  scale_colour_manual(
    name = "Prediction Intervals 95% CI (ribbons)\nMedian Predictions (solid lines)",
    breaks = c("5%PI", "50%PI", "95%PI"),
    labels = c("5%", "50%", "95%"),
    values = c("#2c7fb8", "#df65b0", "#2c7fb8")
  ) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30)) +
  ylab("Obs/Simulated") +
  xlab("Time (h)") +
  theme(axis.title = element_text(size = 20)) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(strip.background = element_rect(fill = "white", colour = "black")) +
  theme(panel.grid.minor = element_line(colour = "gray", linetype = "dotted")) +
  theme(panel.grid.major = element_line(colour = "gray", linetype = "solid")) +
  theme(strip.text = element_text(size = 20, colour = "black"), aspect.ratio = 1) +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  theme(legend.key.width = grid::unit(2, "cm")) +
  scale_x_continuous(breaks = seq(144, 168, 6))




# Bonus for reading that far how to produce shaded PI ribbons with legends etc.
# Benjamin Guiastrennec has a better function instead of the manual process
# still it is a good learning on data wrangling

simdata <- simdata[simdata$TIME > 120, ]

simquantiles <- simdata %>%
  group_by(DVTYPE, TIME) %>%
  filter(!is.na(DV)) %>%
  summarize(
    P50 = median(DV, na.rm = TRUE),
    P05 = quantile(DV, probs = 0.05, type = 7),
    P95 = quantile(DV, probs = 0.95, type = 7),
    P10 = quantile(DV, probs = 0.10, type = 7),
    P90 = quantile(DV, probs = 0.90, type = 7),
    P15 = quantile(DV, probs = 0.15, type = 7),
    P85 = quantile(DV, probs = 0.85, type = 7),
    P25 = quantile(DV, probs = 0.25, type = 7),
    P75 = quantile(DV, probs = 0.75, type = 7),
    P375 = quantile(DV, probs = 0.375, type = 7),
    P625 = quantile(DV, probs = 0.625, type = 7)
  )
require(tidyr)
VPCSTATPLOT2 <- simquantiles %>%
  gather(DVQNAME, quantile, -DVTYPE, -TIME)



VPCSTATPLOT2 <- VPCSTATPLOT2 %>%
  mutate(PI = if_else(DVQNAME == "P05", "90%PI",
                      if_else(DVQNAME == "P10", "80%PI",
                              if_else(DVQNAME == "P15", "70%PI",
                                      if_else(DVQNAME == "P25", "50%PI",
                                              if_else(DVQNAME == "P375", "25%PI",
                                                      
                                                      if_else(DVQNAME == "P95", "90%PI",
                                                              if_else(DVQNAME == "P90", "80%PI",
                                                                      if_else(DVQNAME == "P85", "70%PI",
                                                                              if_else(DVQNAME == "P625", "25%PI",
                                                                                      
                                                                                      if_else(DVQNAME == "P75", "50%PI",
                                                                                              "Median"
                                                                                      )
                                                                              )
                                                                      )
                                                              )
                                                      )
                                              )
                                      )
                              )
                      )
  ))


VPCSTATPLOT2 <- VPCSTATPLOT2 %>%
  mutate(PILIMIT = if_else(DVQNAME == "P05", "PILOW",
                           if_else(DVQNAME == "P10", "PILOW",
                                   if_else(DVQNAME == "P15", "PILOW",
                                           if_else(DVQNAME == "P25", "PILOW",
                                                   if_else(DVQNAME == "P375", "PILOW",
                                                           if_else(DVQNAME == "P625", "PIUP",
                                                                   if_else(DVQNAME == "P95", "PIUP",
                                                                           if_else(DVQNAME == "P90", "PIUP",
                                                                                   if_else(DVQNAME == "P85", "PIUP",
                                                                                           if_else(DVQNAME == "P75", "PIUP", "Median")
                                                                                   )
                                                                           )
                                                                   )
                                                           )
                                                   )
                                           )
                                   )
                           )
  ))



VPCSTATPLOTWIDE <- VPCSTATPLOT2 %>%
  select(-DVQNAME) %>%
  filter(PILIMIT != "Median")
VPCSTATPLOTWIDE <- VPCSTATPLOTWIDE %>%
  tidyr::spread(PILIMIT, quantile)


VPCSTATPLOTWIDEMedian <- VPCSTATPLOT2 %>%
  select(-DVQNAME) %>%
  filter(PILIMIT == "Median")
VPCSTATPLOTWIDEMedian <- VPCSTATPLOTWIDEMedian %>%
  tidyr::spread(PILIMIT, quantile)

ggplot(VPCSTATPLOTWIDE, aes(TIME,
                            ymin = PILOW, ymax = PIUP,
                            fill = PI, col = PI, group = PI
)) +
  facet_wrap(~DVTYPE, scales = "free_y", labeller = label_both) +
  geom_ribbon(alpha = 0.2, col = NA) +
  geom_line(
    data = VPCSTATPLOTWIDEMedian,
    aes(x = TIME, y = Median, group = PI, linetype = "Median"), inherit.aes = FALSE, color = "#009ACDFF",
    size = 2
  ) +
  geom_quasirandom(data = obsdata, aes(TIME, DV), inherit.aes = FALSE, alpha = 0.05, size = 2, width = 0.5) +
  scale_fill_manual(values = c("#009ACDFF", "#009ACDFF", "#009ACDFF", "#009ACDFF", "#009ACDFF")) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  ) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(strip.background = element_rect(fill = "white", colour = "black")) +
  theme(panel.grid.minor = element_line(colour = "gray", linetype = "dotted")) +
  theme(panel.grid.major = element_line(colour = "gray", linetype = "solid")) +
  theme(strip.text = element_text(size = 16, colour = "black")) +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  guides(fill = guide_legend(
    reverse = TRUE,
    title = "",
    override.aes = list(alpha = c(0.2, 0.4, 0.6, 0.8, 1))
  )) +
  guides(linetype = guide_legend(
    title = ""
  )) +
  labs(y = "Prediction Intervals") +
  scale_x_continuous(breaks = seq(144, 168, 6))


#################
# enough simulation or unreal data let us look at real project examples



