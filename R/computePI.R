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
                       breaks = NULL,
                       LLOQ = NULL,
                       filterblq = FALSE, # remove the BLQ just before computing the quantiles, but after "merging" with rep(, sim) 
                       predcorrection = FALSE,
                       predcorrection_islogdv = FALSE,
                       predcorrection_lowerbnd = 0, #not used yet
                       PI = c(0.05, 0.5, 0.95),
                       CIPI = c(0.025, 0.5, 0.975),
                       bootstrapobsdata = FALSE) {
  
  
  if (!is.null(breaks) & is.null(nbins)) nbins = "breaks"
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
  
  originaldatabins <- obsdata
  simdatabins <- simdata
  
  NREP <- nrow(simdatabins) / nrow(originaldatabins)
  message(paste("detected", NREP, "simulations"))
  
  if(bootstrapobsdata){
    warning("warning: bootstrapping observed data ignored not yet implemented")
  }
  
  if (!all.equal(
    originaldatabins$ID,
    simdatabins$ID[1:(nrow(simdatabins) / NREP)]
  )) {
    stop("Error: ID's of your Obs and Sim data are not identical ")
  }
  if (!all.equal(
    as.vector(originaldatabins %>%
              mutate(TIME = !!TIME) %>% select(TIME)),
    as.vector(simdatabins %>%
              filter(!!REPL == min(!!REPL)) %>%
              mutate(TIME = !!TIME) %>% select(TIME))
    ,
    check.attributes = FALSE
  )
  ) {
    stop("Error: Time columns of Obs and Sim data are not identical or not in the same order ")
  }
  
  if (!quo_is_null(LLOQ)) {
    originaldatabins <- originaldatabins %>%
      mutate(LLOQFL = ifelse(!!DV < !!LLOQ, 1, 0))
    
    simdatabins <- simdatabins %>%
      mutate(LLOQFL = ifelse(!!DV < !!LLOQ, 1, 0))
  }
  
  
  
  if (quo_is_null(nbins)) {
    originaldatabins <- originaldatabins %>%
      mutate(BINS = !!TIME)
    # originaldatabins <- originaldatabins %>%
    #   group_by(BINS) %>%
    #   mutate(XMIN = min(!!TIME), XMAX = max(!!TIME)) %>%
    #   mutate(XMED = median(!!TIME)) %>%
    #   mutate(XMID = median(c(XMIN, XMAX))) %>%
    #   mutate(BINSC = paste("(", XMIN, "-", XMAX, "]", sep = "")) %>%
    #   mutate(NOBS = length(DV))
    # 
    # originaldatabins = originaldatabins %>%
    #   mutate(XLEFT = XMIN, XRIGHT = XMAX) #May be improved (more right/left than min/max) if choosing a direction and the value (min or max) of the current bin and the same value in the previous/next one
    
  } else {
    originaldatabins <- originaldatabins %>%
      group_by_at(c(group_vars(.), all.vars(stratify)))
    
    if (!is.null(breaks)) {
      if (is.vector(breaks)) {
        breaks = data.frame(breaks=breaks)
        
        originaldatabins <- originaldatabins %>%
          mutate(BINS = as.numeric(cut(!!TIME, breaks$breaks, include.lowest=T)))
        
      } else {
        cutbreaks <- function(data, breaks) {
          subsetbreaks = semi_join(breaks, data)
          data %>% mutate(BINS = as.numeric(cut(!!TIME, subsetbreaks$breaks, include.lowest=T)))
        }
        
        originaldatabins <- originaldatabins %>% 
          left_join(originaldatabins %>% do(cutbreaks(., breaks))) #join to preserver original order
      }
      
      breaks = breaks %>% 
        group_by_at(intersect(group_vars(originaldatabins), names(breaks %>% select(-breaks)))) %>%
        mutate(BINS = seq_len(n()),
               XLEFT = breaks,
               XRIGHT = c(breaks[-1],NA))
      
    } else { 
      if (is.null(style)) {
        originaldatabins <- originaldatabins %>%
          mutate(BINS = ntile(!!TIME, !!nbins))
      } else {
        originaldatabins <- originaldatabins %>%
          mutate(BINS = as.numeric(cut(!!TIME, classIntervals(!!TIME, n=!!nbins, style=style)$brks, include.lowest=T)))
      }
    }
    
    originaldatabins <- originaldatabins %>%
      group_by(BINS, add = TRUE) %>%
      mutate(XMIN = min(!!TIME), XMAX = max(!!TIME)) %>%
      mutate(XMED = median(!!TIME)) %>%
      mutate(XMID = median(c(XMIN, XMAX))) %>%
      mutate(BINSC = paste("(", XMIN, "-", XMAX, "]", sep = "")) %>%
      mutate(NOBS = length(DV))
    
    if (!is.null(breaks)) {
      originaldatabins = originaldatabins %>%
        left_join(breaks %>% select(-breaks)) 
    } else {
      originaldatabins = originaldatabins %>%
        mutate(XLEFT = XMIN, XRIGHT = XMAX) #May be improved (more right/left than min/max) if choosing a direction and the value (min or max) of the current bin and the same value in the previous/next one
    }
  }
  
  simdatabins$XMIN <- rep(originaldatabins$XMIN, time = NREP)
  simdatabins$XMED <- rep(originaldatabins$XMED, time = NREP)
  simdatabins$XMID <- rep(originaldatabins$XMID, time = NREP)
  simdatabins$XMAX <- rep(originaldatabins$XMAX, time = NREP)
  simdatabins$XLEFT <- rep(originaldatabins$XLEFT, time = NREP)
  simdatabins$XRIGHT <- rep(originaldatabins$XRIGHT, time = NREP)
  simdatabins$BINS <- rep(originaldatabins$BINS, time = NREP)
  simdatabins$BINSC <- rep(originaldatabins$BINSC, time = NREP)
  simdatabins$NOBS <- rep(originaldatabins$NOBS, time = NREP)
  
  if (filterblq & !quo_is_null(LLOQ)) {
    originaldatabins = originaldatabins %>% filter(LLOQFL==0)
    simdatabins = simdatabins %>% filter(LLOQFL==0)
  }
  
  if (quo_is_null(LLOQ)) {
    percentblqobs <- originaldatabins %>%
      dplyr::mutate(percentblq = NA)
    originaldatabins <- left_join(originaldatabins, percentblqobs)
  }
  if (!quo_is_null(LLOQ)) {
    percentblqobs <- originaldatabins %>%
      group_by(BINSC, XMIN, XMAX, XMED, XMID, XLEFT, XRIGHT, NOBS)
    
    
    percentblqobs <- percentblqobs %>%
      group_by_at(c(group_vars(.), all.vars(stratify)))
    
    
    percentblqobs <- percentblqobs %>%
      dplyr::mutate(percentblq = 100 * mean(LLOQFL))
    originaldatabins <- left_join(originaldatabins, percentblqobs)
  }
  
  if (predcorrection) {
    originaldatabins <- originaldatabins %>%
      dplyr::group_by(BINSC)
    
    
    originaldatabins <- originaldatabins %>%
      group_by_at(c(group_vars(.), all.vars(stratify)))
    
    
    originaldatabins <- originaldatabins %>%
      dplyr::mutate(mp = median(PRED))
    
    if (!predcorrection_islogdv) {
      originaldatabins <- originaldatabins %>%
        dplyr::mutate(dvc = !!DV * mp / PRED)
    }
    if (predcorrection_islogdv) {
      originaldatabins <- originaldatabins %>%
        dplyr::mutate(dvc = !!DV + (mp - PRED))
    }
  }
  
  if (!predcorrection) {
    originaldatabins <- originaldatabins %>%
      dplyr::mutate(dvc = !!DV)
  }
  
  if (quo_is_null(LLOQ)) {
    PIobs <- originaldatabins %>%
      group_by(BINS, BINSC, percentblq, XMIN, XMAX, XMED, XMID, XLEFT, XRIGHT, NOBS)
    
    
    PIobs <- PIobs %>%
      group_by_at(c(group_vars(.), all.vars(stratify)))
    
    
    PIobs <- PIobs %>%
      dplyr::summarize(
        PLPI = quantile(dvc, probs = PI[1]),
        PMPI = quantile(dvc, probs = PI[2]),
        PUPI = quantile(dvc, probs = PI[3])
      )
    PIobs <- tidyr::gather(
      PIobs, quantilename, quantilevalue,
      PLPI, PMPI, PUPI
    )
  }
  if (!quo_is_null(LLOQ)) {
    PIobs <- originaldatabins %>%
      mutate(LLOQ = !!LLOQ) %>%
      group_by(BINS, BINSC, percentblq, XMIN, XMAX, XMED, XMID, XLEFT, XRIGHT, NOBS, LLOQ)
    
    PIobs <- PIobs %>%
      group_by_at(c(group_vars(.), all.vars(stratify)))
    
    
    PIobs <- PIobs %>%
      dplyr::summarize(
        PLPI = quantile_cens(dvc, p = PI[1], limit = LLOQ),
        PMPI = quantile_cens(dvc, p = PI[2], limit = LLOQ),
        PUPI = quantile_cens(dvc, p = PI[3], limit = LLOQ)
      )
    PIobs <- tidyr::gather(
      PIobs, quantilename, quantilevalue,
      PLPI, PMPI, PUPI
    )
  }
  
  PIobs <- PIobs
  
  if (predcorrection) {
    simdatabins$mp <- rep(originaldatabins$mp, time = NREP)
    if (!predcorrection_islogdv) {
      simdatabins <- simdatabins %>%
        dplyr::mutate(dvc = !!DV * mp / PRED)
    }
    if (predcorrection_islogdv) {
      simdatabins <- simdatabins %>%
        dplyr::mutate(dvc = !!DV + (mp - PRED))
    }
  }
  
  if (!predcorrection) {
    simdatabins <- simdatabins %>%
      dplyr::mutate(dvc = !!DV)
  }
  
  PICIsim <- simdatabins %>%
    group_by(!!REPL, BINS, BINSC, XMIN, XMAX, XMED, XMID, XLEFT, XRIGHT, NOBS)
  
  
  PICIsim <- PICIsim %>%
    group_by_at(c(group_vars(.), all.vars(stratify)))
  
  
  PICIsim <- PICIsim %>%
    dplyr::summarize(
      PLPI = quantile(dvc, probs = PI[1]),
      PMPI = quantile(dvc, probs = PI[2]),
      PUPI = quantile(dvc, probs = PI[3])
    )
  SIMPIPI <- tidyr::gather(PICIsim, quantilename, quantilevalue, PLPI, PMPI, PUPI)
  
  VPCSTAT <- SIMPIPI %>%
    group_by(quantilename, BINS, BINSC, XMIN, XMAX, XMED, XMID, XLEFT, XRIGHT, NOBS)
  
  
  VPCSTAT <- VPCSTAT %>%
    group_by_at(c(group_vars(.), all.vars(stratify)))
  
  
  VPCSTAT <- VPCSTAT %>%
    dplyr::summarize(
      QELCI = quantile(quantilevalue, probs = CIPI[1]),
      QEMCI = quantile(quantilevalue, probs = CIPI[2]),
      QEUCI = quantile(quantilevalue, probs = CIPI[3])
    )
  
  
  if (!quo_is_null(LLOQ)) {
    percentblqsim <- simdatabins %>%
      group_by(!!REPL, BINSC, XMIN, XMAX, XMED, XMID, XLEFT, XRIGHT, NOBS)
    
    
    percentblqsim <- percentblqsim %>%
      group_by_at(c(group_vars(.), all.vars(stratify)))
    
    
    percentblqsim <- percentblqsim %>%
      dplyr::summarize(percentblq = mean(LLOQFL))
    percentblqsim <- tidyr::gather(percentblqsim, blqquantilename, blqquantilevalue, percentblq)
    
    percentblqsimSTAT <- percentblqsim %>%
      group_by(blqquantilename, BINSC, XMIN, XMAX, XMED, XMID, XLEFT, XRIGHT, NOBS)
    
    
    percentblqsimSTAT <- percentblqsimSTAT %>%
      group_by_at(c(group_vars(.), all.vars(stratify)))
    
    
    percentblqsimSTAT <- percentblqsimSTAT %>%
      dplyr::summarize(
        percentblqQELCI = 100 * quantile(blqquantilevalue, probs = CIPI[1]),
        percentblqQEMCI = 100 * quantile(blqquantilevalue, probs = CIPI[2]),
        percentblqQEUCI = 100 * quantile(blqquantilevalue, probs = CIPI[3])
      )
    percentblqsimSTAT <- percentblqsimSTAT
    VPCSTAT <- left_join(VPCSTAT, percentblqsimSTAT)
  }
  VPCSTAT <- left_join(VPCSTAT, PIobs)
  
  VPCSTAT %>% ungroup()
}
