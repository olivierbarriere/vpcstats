require(ggplot2)
require(dplyr)
require(rlang)
require(vpc) # to compare test and augment on it
require(classInt)


compute.PI <- function(obsdata,
                       simdata,
                       stratify = NULL, #formula
                       TIME = TIME,
                       DV = DV,
                       REPL = REP,
                       nbins = NULL,
                       style = "ntile", # "ntile", ("fixed", "sd", "equal", "pretty", "quantile", "kmeans", "hclust", "bclust", "fisher", "jenks")-> ClassInt, "breaks"
                       breaks = NULL, #either a vector or a data.frame with a column "breaks"
                       LLOQ = NULL,
                       filterblq = FALSE, # remove the BLQ just before computing the quantiles, but after "merging" with rep(, sim) 
                       predcorrection = FALSE,
                       predcorrection_islogdv = FALSE,
                       predcorrection_lowerbnd = 0, #not used yet
                       PI = c(0.05, 0.5, 0.95),
                       CIPI = c(0.025, 0.5, 0.975),
                       bootstrapobsdata = FALSE) {
  
  if (!is.null(breaks) & is.null(nbins)) {
    nbins <- "breaks"
  }
  nbins <- enquo(nbins)
  if (quo_is_null(nbins)) {
    message("nbins is null no binning done")
  }
  TIME <- enquo(TIME)
  DV <- enquo(DV)
  REPL <- enquo(REPL)
  LLOQ <- enquo(LLOQ)
  if (quo_is_null(LLOQ)) {
    message("LLOQ is null")
  }
  if (filterblq) {
    message("Equal censoring of observed and simulated data")
  }
  
  if(bootstrapobsdata){
    warning("warning: bootstrapping observed data ignored not yet implemented")
  }
  
  originaldatabins <- obsdata
  simdatabins <- simdata
  
  NREP <- nrow(simdatabins) / nrow(originaldatabins)
  message(paste("detected", NREP, "simulations"))
  
  if (NREP %%1 != 0) {
    stop("Error: Sim dataset is not a multiple of obs dataset")
  }
  
  if (!all.equal(
    originaldatabins %>% pull(ID),
    simdatabins %>% filter(!!REPL == min(!!REPL)) %>% pull(ID))) {
    stop("Error: ID's of your Obs and Sim data are not identical")
  }
  
  if (!all.equal(
    originaldatabins %>% mutate(TIME = !!TIME) %>% pull(TIME),
    simdatabins %>% filter(!!REPL == min(!!REPL)) %>% mutate(TIME = !!TIME) %>% pull(TIME),
    check.attributes = FALSE)) {
    stop("Error: Time columns of Obs and Sim data are not identical or not in the same order")
  }
  
  originaldatabins <- originaldatabins %>%
    group_by_at(all.vars(stratify))
  simdatabins <- simdatabins %>%
    group_by_at(all.vars(stratify))
  
  if (!quo_is_null(LLOQ)) {
    originaldatabins <- originaldatabins %>%
      mutate(LLOQFL = ifelse(!!DV < !!LLOQ, 1, 0))
    
    simdatabins <- simdatabins %>%
      mutate(LLOQFL = ifelse(!!DV < !!LLOQ, 1, 0))
  }
  
  if (quo_is_null(nbins)) {
    originaldatabins <- originaldatabins %>%
      mutate(BIN = !!TIME)
  } else {
    if (!is.null(breaks)) {
      if (is.vector(breaks)) {
        breaks <- data.frame(breaks=breaks)
        
        originaldatabins <- originaldatabins %>%
          mutate(BIN = as.numeric(cut(!!TIME, breaks$breaks, include.lowest=T)))
        
      } else {
        cutbreaks <- function(data, breaks) {
          subsetbreaks = semi_join(breaks, data)
          data %>% mutate(BIN = as.numeric(cut(!!TIME, subsetbreaks$breaks, include.lowest=T)))
        }
        
        originaldatabins <- originaldatabins %>% 
          left_join(originaldatabins %>% do(cutbreaks(., breaks))) #join to preserver original order
      }
      
      breaks <- breaks %>% 
        group_by_at(intersect(group_vars(originaldatabins), names(breaks %>% select(-breaks)))) %>%
        mutate(BIN = seq_len(n()),
               XLEFT = breaks,
               XRIGHT = c(breaks[-1],NA)) %>%
        select(-breaks)
      
    } else { 
      if (style == "ntile") {
        originaldatabins <- originaldatabins %>%
          mutate(BIN = ntile(!!TIME, !!nbins))
      } else if (style %in% c("fixed", "sd", "equal", "pretty", "quantile", "kmeans", "hclust", "bclust", "fisher", "jenks")) {
        originaldatabins <- originaldatabins %>%
          mutate(BIN = as.numeric(cut(!!TIME, classIntervals(!!TIME, n=!!nbins, style=style)$brks, include.lowest=T)))
      } else {
        stop("Error: Unknown style")
      }
    }
  }
  
  originaldatabins <- originaldatabins %>%
    group_by(BIN, add = TRUE) %>%
    mutate(XMIN = min(!!TIME),
           XMAX = max(!!TIME),
           XMED = median(!!TIME),
           XMID = median(c(XMIN, XMAX)),
           NOBS = length(DV))
  
  if (is.null(breaks)) {
    breaks <- originaldatabins %>%
      group_by_at(all.vars(stratify)) %>%
      distinct(BIN, XMIN, XMAX) %>%
      arrange_at(c(group_vars(.),"BIN")) %>%
      #mutate(XLEFT = (XMIN+c(-Inf,XMAX[-length(XMAX)]))/2, #Infinite left and right boundaries
      #       XRIGHT = (XMAX+c(XMIN[-1],Inf))/2) %>%
      mutate(XLEFT = (XMIN+c(XMIN[1],XMAX[-length(XMAX)]))/2, #Min and max left and right boundaries
             XRIGHT = (XMAX+c(XMIN[-1],XMAX[length(XMAX)]))/2)
  }
  
  if (quo_is_null(LLOQ)) {
    bins <- originaldatabins %>%
      distinct(BIN, XMIN, XMAX, XMED, XMID, NOBS)
  } else {
    bins <- originaldatabins %>%
      mutate(LLOQ = !!LLOQ) %>%
      distinct(BIN, XMIN, XMAX, XMED, XMID, NOBS, LLOQ) #Warning, can produce unexpected results if LLOQ is not unique within each strata
  }
  bins <- bins %>% 
    left_join(breaks)
  
  simdatabins$BIN <- rep(originaldatabins$BIN, time = NREP)
  
  if (filterblq & !quo_is_null(LLOQ)) {
    originaldatabins <- originaldatabins %>% filter(LLOQFL==0)
    simdatabins <- simdatabins %>% filter(LLOQFL==0)
  }
  
  if (predcorrection) {
    originaldatabins <- originaldatabins %>%
      dplyr::mutate(MEDPRED = median(PRED))
    
    if (!predcorrection_islogdv) {
      originaldatabins <- originaldatabins %>%
        dplyr::mutate(DVC = !!DV * MEDPRED / PRED)
    }
    if (predcorrection_islogdv) {
      originaldatabins <- originaldatabins %>%
        dplyr::mutate(DVC = !!DV + (MEDPRED - PRED))
    }
  } else {
    originaldatabins <- originaldatabins %>%
      dplyr::mutate(DVC = !!DV)
  }
  
  if (quo_is_null(LLOQ)) {
    PIobs <- originaldatabins %>%
      dplyr::summarize(
        PLPI = quantile(DVC, probs = PI[1]),
        PMPI = quantile(DVC, probs = PI[2]),
        PUPI = quantile(DVC, probs = PI[3])) %>%
      tidyr::gather(QNAME, QVALUE,
                    PLPI, PMPI, PUPI)
  } else {
    PIobs <- originaldatabins %>%
      dplyr::mutate(PCTBLQ_VALUE = 100 * mean(LLOQFL)) %>%
      group_by(PCTBLQ_VALUE, add = TRUE) %>%
      #mutate(LLOQ = !!LLOQ) %>%
      #group_by(LLOQ, add = TRUE) %>%
      dplyr::summarize(
        PLPI = quantile_cens(DVC, p = PI[1], limit = !!LLOQ),
        PMPI = quantile_cens(DVC, p = PI[2], limit = !!LLOQ),
        PUPI = quantile_cens(DVC, p = PI[3], limit = !!LLOQ)) %>%
      tidyr::gather(QNAME, QVALUE,
                    PLPI, PMPI, PUPI)
  }
  
  
  if (predcorrection) {
    simdatabins$MEDPRED <- rep(originaldatabins$MEDPRED, time = NREP)
    if (!predcorrection_islogdv) {
      simdatabins <- simdatabins %>%
        dplyr::mutate(DVC = !!DV * MEDPRED / PRED)
    }
    if (predcorrection_islogdv) {
      simdatabins <- simdatabins %>%
        dplyr::mutate(DVC = !!DV + (MEDPRED - PRED))
    }
  } else {
    simdatabins <- simdatabins %>%
      dplyr::mutate(DVC = !!DV)
  }
  
  SIMPIPI <- simdatabins %>%
    group_by(BIN, !!REPL, add = TRUE) %>%
    dplyr::summarize(
      PLPI = quantile(DVC, probs = PI[1]),
      PMPI = quantile(DVC, probs = PI[2]),
      PUPI = quantile(DVC, probs = PI[3])) %>%
    tidyr::gather(QNAME, QVALUE,
                  PLPI, PMPI, PUPI)
  
  VPCSTAT <- SIMPIPI %>%
    group_by(QNAME, add = TRUE) %>% #summarise will drop the last group (here !!REPL) https://github.com/tidyverse/dplyr/issues/862
    dplyr::summarize(
      QLCI = quantile(QVALUE, probs = CIPI[1]),
      QMCI = quantile(QVALUE, probs = CIPI[2]),
      QUCI = quantile(QVALUE, probs = CIPI[3]))
  
  
  if (!quo_is_null(LLOQ)) {
    percentblqsimSTAT <- simdatabins %>%
      group_by(BIN, !!REPL, add = TRUE) %>%
      dplyr::summarize(PERCENTBLQ = 100 * mean(LLOQFL)) %>% #summarise will drop the last group (here !!REPL)
      dplyr::summarize(PCTBLQ_LCI = quantile(PERCENTBLQ, probs = CIPI[1]),
                       PCTBLQ_MCI = quantile(PERCENTBLQ, probs = CIPI[2]),
                       PCTBLQ_UCI = quantile(PERCENTBLQ, probs = CIPI[3]))
    
    VPCSTAT <- left_join(VPCSTAT, percentblqsimSTAT)
  }
  
  VPCSTAT <- left_join(VPCSTAT, PIobs)
  VPCSTAT <- left_join(VPCSTAT, bins)
  
  VPCSTAT %>% ungroup()
}
