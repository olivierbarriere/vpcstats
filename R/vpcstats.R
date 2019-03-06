#' vpcstats
#'
#' @param obsdata Obs dataset
#' @param simdata Sim dataset
#' @param stratify Stratify formula: right hand side only
#' @param TIME Time column: a name or an expression
#' @param DV Dependent Variable column: a name or an expression
#' @param REPL Replication Index column: a name
#' @param LLOQ LLOQ column: a name, an expression, or a fixed value
#' @param NBINS NBINS column: a name, an expression, or a fixed value
#' @param bin_by_strata Flag to indicate if the binning is by strata (default) or global
#' @param bin_style Style of the binning procedure: "ntile" (default), {"fixed", "sd", "equal", "pretty", "quantile", "kmeans", "hclust", "bclust", "fisher", or "jenks"} (ClassInterval styles), "breaks" (user defined breaks)
#' @param breaks User defined breaks: either a vector or a wide data.frame with a column "breaks", and other columns corresponding to the stratify variables
#' @param cut_right Flag to indicate if the intervals should be closed on the right (and open on the left) or vice versa. Default to False, unlike the cut function.
#' @param quantile_type Quantile algorithm: an integer between 4 and 9 selecting one of the nine quantile algorithms
#' @param filterblq Flag to remove the BLQ values just before computing the quantiles, but after "merging" with rep(, sim) 
#' @param predcorrection Flag to indicate if the VPC has to be pred corrected or not
#' @param predcorrection_islogdv Flag to indicate if the data was log transformed
#' @param predcorrection_lowerbnd Not used yet
#' @param PI Prediction Intervals: vector of any number of  values
#' @param CI Confidence Intervals around the Prediction Intervals : vector of any number of value
#' @param bootstrapobsdata Not used yet
#' 
#' @rawNamespace importFrom("stats", "as.formula", "median", "quantile")
#' @rawNamespace importFrom("rlang", ".data")
#' @rawNamespace import(dplyr, except = c(last,between,first))
#' @rawNamespace import(data.table, except = c(last,between,first))
#' @export
#' @examples
#' 
#' library(vpc)
#' library(vpcstats)
#' exampleobs <- simple_data$obs
#' exampleobs <- exampleobs[exampleobs$MDV == 0, ]
#' examplesim <- simple_data$sim
#' exampleobs$PRED <- examplesim[1:nrow(exampleobs), "PRED"]
#' exampleobs$REP <- 0
#' examplesim <- examplesim[examplesim$MDV == 0, ]#
#' exampleobs$LLOQ <- ifelse(exampleobs$ISM == 0, 100, 25)
#' examplesim$LLOQ <- ifelse(examplesim$ISM == 0, 100, 25)
#' VPCDATA<- vpcstats(
#' obsdata = exampleobs, simdata = examplesim, stratify = ~ISM,
#' NBINS = NULL, LLOQ = LLOQ)

vpcstats <- function(obsdata = NULL,
                      simdata,
                      stratify = NULL,
                      TIME = TIME,
                      DV = DV,
                      REPL = REP,
                      LLOQ = NULL,
                      NBINS = NULL,
                      bin_by_strata = TRUE,
                      bin_style = "ntile",
                      breaks = NULL,
                      cut_right = FALSE,
                      quantile_type = 7,
                      filterblq = FALSE,
                      predcorrection = FALSE,
                      predcorrection_islogdv = FALSE,
                      predcorrection_lowerbnd = 0,
                      PI = c(0.05, 0.5, 0.95),
                      CI = c(0.025, 0.5, 0.975),
                      bootstrapobsdata = FALSE) {
REP = ID = BIN = LLOQFL = PRED = MEDPRED = DVC = SIM = NOBS = NULL
XMIN = XMAX = XMED = XMEAN = XLEFT= XRIGHT = NULL

  NBINS <- rlang::enquo(NBINS)
  if ( rlang::quo_is_null(NBINS)) {
    message("No binning done")
  }
  TIME <-  rlang::enquo(TIME)
  DV <-  rlang::enquo(DV)
  REPL <-  rlang::enquo(REPL)
  LLOQ <-  rlang::enquo(LLOQ)
  if ( rlang::quo_is_null(LLOQ)) {
    message("LLOQ is not defined")
  }
  if (filterblq) {
    message("Equal censoring of observed and simulated data")
  }
  
  if(bootstrapobsdata){
    warning("Warning: bootstrapping observed data ignored not yet implemented")
  }
  
  simdatabins <- simdata
  if (!is.null(obsdata)) {
    obsdatabins <- obsdata
  } else {
    obsdatabins <-  simdatabins %>%
      dplyr::filter(!!REPL==min(!!REPL)) %>%
      dplyr::mutate(!!rlang::quo_name(DV) := NA)
  }
  
  NREP <- nrow(simdatabins) / nrow(obsdatabins)
  message(paste("Detected", NREP, "simulations"))
  
  if (NREP %%1 != 0) {
    stop("Error: Sim dataset length is not a multiple of obs dataset length")
  }
  
  if ("ID" %in% names(obsdatabins) & "ID" %in% names(simdatabins)) {
    if (!isTRUE(all.equal(
      obsdatabins %>% pull(ID),
      simdatabins %>% filter(!!REPL == min(!!REPL)) %>% pull(ID)))) {
      stop("Error: ID's of your Obs and Sim data are not identical")
    }
  }
  
  if (!isTRUE(all.equal(
    obsdatabins %>% mutate(TIME = !!TIME) %>% pull(TIME),
    simdatabins %>% filter(!!REPL == min(!!REPL)) %>% mutate(TIME = !!TIME) %>% pull(TIME),
    check.attributes = FALSE))) {
    stop("Error: Time columns of Obs and Sim data are not identical or not in the same order")
  }
  
  stratifyvars <- all.vars(stratify)
  
  if (!rlang::quo_is_null(LLOQ)) {
    obsdatabins <- obsdatabins %>%
      mutate(LLOQFL = ifelse(!!DV < !!LLOQ, 1, 0))
    
    simdatabins <- simdatabins %>%
      mutate(LLOQFL = ifelse(!!DV < !!LLOQ, 1, 0))
  }
  
  if (bin_by_strata) {
    obsdatabins <- obsdatabins %>%
      group_by_at(stratifyvars)
  }
  simdatabins <- simdatabins %>%
    group_by_at(stratifyvars)
  
  if (rlang::quo_is_null(NBINS) & is.null(breaks)) {
    obsdatabins <- obsdatabins %>%
      mutate(BIN = !!TIME)
  } else {
    if (!is.null(breaks)) {
      if (is.vector(breaks)) {
        breaks <- data.frame(breaks=breaks)
        
        obsdatabins <- obsdatabins %>%
          mutate(BIN = as.numeric(cut(!!TIME, breaks$breaks, right=cut_right, include.lowest=T)))
        
      } else {
        cutbreaks <- function(data, breaks) {
          subsetbreaks = semi_join(breaks, data, by=intersect(names(breaks), names(data)))
          data %>% mutate(BIN = as.numeric(cut(!!TIME, subsetbreaks$breaks, right=cut_right, include.lowest=T)))
        }
        
        obsdatabins <- obsdatabins %>% 
          left_join(obsdatabins %>% do(cutbreaks(.data, breaks)), #join to preserver original order
                    by=names(obsdatabins)) 
      }
      
    } else { 
      if (!is.null(bin_style) && bin_style == "ntile") {
        obsdatabins <- obsdatabins %>%
          mutate(BIN = ntile(!!TIME, !!NBINS))
      } else if (!is.null(bin_style) && bin_style %in% c("fixed", "sd", "equal", "pretty", "quantile", "kmeans", "hclust", "bclust", "fisher", "jenks")) {
        obsdatabins <- obsdatabins %>%
          mutate(BIN = as.numeric(cut(!!TIME, classInt::classIntervals(!!TIME, n=!!NBINS, style=bin_style)$brks, right=cut_right, include.lowest=T)))
        
        breaks <- obsdatabins %>%
          mutate(TIME = !!TIME,
                 NBINS = !!NBINS) %>%
          select_at(c(stratifyvars, "TIME","NBINS")) %>%
          do(data.frame(breaks = classInt::classIntervals(.data$TIME, n=unique(.data$NBINS), style=bin_style)$brks))
      } else {
        stop("Error: Unknown binning style")
      }
    }
  }
  
  obsdatabins <- obsdatabins %>%
    group_by_at(stratifyvars) %>% #Now that the binning is done, unconditionnaly group by strata (even if bin_by_strata is false)
    group_by(BIN, add = TRUE) %>%
    mutate(XMIN = min(!!TIME),
           XMAX = max(!!TIME),
           XMED = median(!!TIME),
           XMEAN = mean(!!TIME),
           NOBS = length(DV))
  
  if (is.null(breaks)) {
    breaks <- obsdatabins %>%
      group_by_at(stratifyvars) %>%
      distinct(BIN, XMIN, XMAX)
      groupvars<- group_vars(breaks)
      breaks <- breaks %>%
      arrange_at(c(groupvars,"BIN")) %>%
      #mutate(XLEFT = (XMIN+c(-Inf,XMAX[-length(XMAX)]))/2, #Infinite left and right boundaries
      #       XRIGHT = (XMAX+c(XMIN[-1],Inf))/2) %>%
      mutate(XLEFT = (XMIN+c(XMIN[1],XMAX[-length(XMAX)]))/2, #Min and max left and right boundaries
             XRIGHT = (XMAX+c(XMIN[-1],XMAX[length(XMAX)]))/2) %>%
      select(-XMIN, -XMAX)
  } else {
    breaks <- breaks %>% 
      group_by_at(intersect(stratifyvars, names(breaks %>% select(-breaks)))) %>%
      mutate(BIN = seq_len(n()),
             XLEFT = breaks,
             XRIGHT = c(breaks[-1],NA)) %>%
      select(-breaks)
  }
  
  if (rlang::quo_is_null(LLOQ)) {
    BINS <- obsdatabins %>%
      distinct(BIN, XMIN, XMAX, XMED, XMEAN, NOBS)
  } else {
    BINS <- obsdatabins %>%
      mutate(LLOQ = !!LLOQ) %>%
      distinct(BIN, XMIN, XMAX, XMED, XMEAN, NOBS, LLOQ) #Warning, can produce unexpected results if LLOQ is not unique within each strata
  }
  
  BINS <- BINS %>%
    left_join(breaks, by=intersect(names(BINS), names(breaks)))
  
  BINS <- BINS %>%
    mutate(XMID = median(c(XMIN, XMAX)),
           XCENTER = median(c(XLEFT, XRIGHT)))
  
  breaks <- breaks %>% filter(BIN %in% unique(BINS$BIN)) 
  stratifyvarsbreaks <- intersect(stratifyvars, names(breaks))
  breaks %>%
    group_by_at(stratifyvarsbreaks) %>%
    do({
      if (length(stratifyvarsbreaks)>0) {
        uv <- t(unique(.data[,stratifyvarsbreaks]))
        msg1 <- paste0(paste(apply(uv, 2, function(x) paste(rownames(uv), x, sep="=")), collapse=", "), " ")
      } else {
        msg1 <- ""
      }
      msg2 <- paste0("[",length(.data$XRIGHT),"] ")
      msg3 <- paste(signif(sort(c(min(.data$XLEFT), .data$XRIGHT)),3), collapse=" < ")
      message(paste0(msg1, msg2, msg3))
      
      data.frame()
    })
  
  simdatabins$BIN <- rep(obsdatabins$BIN, time = NREP)
  
  if (filterblq & !rlang::quo_is_null(LLOQ)) {
    obsdatabins <- obsdatabins %>% filter(LLOQFL==0)
    simdatabins <- simdatabins %>% filter(LLOQFL==0)
  }
  
  if (predcorrection) {
    obsdatabins <- obsdatabins %>%
      dplyr::mutate(MEDPRED = median(PRED))
    
    if (!predcorrection_islogdv) {
      obsdatabins <- obsdatabins %>%
        dplyr::mutate(DVC = !!DV * MEDPRED / PRED)
    }
    if (predcorrection_islogdv) {
      obsdatabins <- obsdatabins %>%
        dplyr::mutate(DVC = !!DV + (MEDPRED - PRED))
    }
  } else {
    obsdatabins <- obsdatabins %>%
      dplyr::mutate(DVC = !!DV)
  }
  
  if (rlang::quo_is_null(LLOQ)) {
    PIobs <- data.table::as.data.table(obsdatabins)[, 
                                                    list(QNAME = paste0("",100*PI,"%PI"),
                                                         RAWOBS = quantile(DVC, probs = PI, type = quantile_type, na.rm = T)),
                                                    by = c(stratifyvars,"BIN")]    
    
  } else {
    PIobs <- data.table::as.data.table(obsdatabins %>% 
                                         mutate(LLOQ = !!LLOQ))[, #Don't know how to use !! inside data.table  
                                                                list(QNAME=paste0("",100*PI,"%PI"),
                                                                     RAWOBS=as.numeric(quantile_cens(DVC, p = PI, limit = LLOQ, type = quantile_type, na.rm = T))),
                                                                by = c(stratifyvars,"BIN")]
    PCTBLQobs <- as.data.table(obsdatabins)[, 
                                            list(QNAME = "PercentBLQ",
                                                 RAWOBS = 100 * mean(LLOQFL)),
                                            by = c(stratifyvars,"BIN")]
  }
  
  if (predcorrection) {
    simdatabins$MEDPRED <- rep(obsdatabins$MEDPRED, time = NREP)
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
  
  PIsim <- as.data.table(simdatabins)[, 
                                      list(QNAME=paste0("",100*PI,"%PI"),
                                           SIM=quantile(DVC, probs = PI, type = quantile_type, na.rm = T)),
                                      by = c(stratifyvars,"BIN",quo_name(REPL))]
  
  PI <- PIsim[, 
              list(CI=paste0("SIM",100*CI,"%CI"),
                   Q=quantile(SIM, probs = CI, type = quantile_type, na.rm = T)),
              by = c(stratifyvars,"BIN","QNAME")]
  
  PI <- dcast(PI, as.formula(paste0(paste(c(stratifyvars, "BIN", "QNAME"), collapse="+"), "~CI")), value.var = "Q")
  
  PI <- merge(PI, PIobs, all.x=TRUE, by=c(stratifyvars,"BIN","QNAME"))
  
  if (!rlang::quo_is_null(LLOQ)) {
    PCTBLQsim <- as.data.table(simdatabins)[, 
                                            list(QNAME = "PercentBLQ",
                                                 SIM = 100 * mean(LLOQFL)),
                                            by = c(stratifyvars,"BIN",quo_name(REPL))]
    
    PCTBLQ <- PCTBLQsim[, 
                        list(CI=paste0("SIM",100*CI,"%CI"),
                             Q=quantile(SIM, probs = CI, type = quantile_type, na.rm = T)),
                        by = c(stratifyvars,"BIN","QNAME")]
    
    PCTBLQ <- dcast(PCTBLQ, as.formula(paste0(paste(c(stratifyvars, "BIN","QNAME"), collapse="+"), "~CI")), value.var = "Q")
    
    PCTBLQ <- merge(PCTBLQ, PCTBLQobs, all.x=TRUE, by=c(stratifyvars,"BIN","QNAME"))
    
  } else {
    PCTBLQ <- NULL 
  }
  
  PI <- as.data.frame(merge(PI, BINS, all.x=TRUE, by=c(stratifyvars,"BIN")))
  PCTBLQ <- if (!is.null(PCTBLQ)) as.data.frame(merge(PCTBLQ, BINS, all.x=TRUE, by=c(stratifyvars,"BIN"))) else NULL
  BINS <- as.data.frame(BINS)
  
  list(PI = PI,
       PCTBLQ = PCTBLQ,
       BINS = BINS)
}