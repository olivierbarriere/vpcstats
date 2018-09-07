require(dplyr)
require(data.table)
require(rlang)
require(classInt)

#' Compute PI
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
#' @param style Style of the binning procedure: "ntile" (default), {"fixed", "sd", "equal", "pretty", "quantile", "kmeans", "hclust", "bclust", "fisher", or "jenks"} (ClassInterval styles), "breaks" (user defined breaks)
#' @param breaks User defined breaks: either a vector or a wide data.frame with a column "breaks", and other colums corresponding to the stratify variables
#' @param cut_right Flag to indicate if the intervals should be closed on the right (and open on the left) or vice versa. Default to False, unlike the cut function.
#' @param filterblq Flag to remove the BLQ values just before computing the quantiles, but after "merging" with rep(, sim) 
#' @param predcorrection Flag to indicate if the VPC has to be pred corrected or not
#' @param predcorrection_islogdv Flag to indicate if the data was log transformed
#' @param predcorrection_lowerbnd Not used yet
#' @param PI Prediction Intervals: vector of any number of  values
#' @param CI Confidence Intervals around the Prediction Intervals : vector of any number of value
#' @param bootstrapobsdata Not used yet
#' @param sort_output Output sorted by strata, bin, and quantile of dv
#'
#' @return
#' @export
#'
#' @examples


compute.PI <- function(obsdata = NULL,
                       simdata,
                       stratify = NULL,
                       TIME = TIME,
                       DV = DV,
                       REPL = REP,
                       LLOQ = NULL,
                       NBINS = NULL,
                       bin_by_strata = TRUE,
                       style = "ntile",
                       breaks = NULL,
                       cut_right = FALSE,
                       filterblq = FALSE,
                       predcorrection = FALSE,
                       predcorrection_islogdv = FALSE,
                       predcorrection_lowerbnd = 0,
                       PI = c(0.05, 0.5, 0.95),
                       CI = c(0.025, 0.5, 0.975),
                       bootstrapobsdata = FALSE,
                       sort_output = TRUE) {
  
  NBINS <- enquo(NBINS)
  if (quo_is_null(NBINS)) {
    message("No binning done")
  }
  TIME <- enquo(TIME)
  DV <- enquo(DV)
  REPL <- enquo(REPL)
  LLOQ <- enquo(LLOQ)
  if (quo_is_null(LLOQ)) {
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
      filter(!!REPL==min(!!REPL)) %>%
      mutate(!!quo_name(DV) := NA)
  }
  
  NREP <- nrow(simdatabins) / nrow(obsdatabins)
  message(paste("Detected", NREP, "simulations"))
  
  if (NREP %%1 != 0) {
    stop("Error: Sim dataset length is not a multiple of obs dataset length")
  }
  
  if ("ID" %in% names(obsdatabins) & "ID" %in% names(simdatabins)) {
    if (!all.equal(
      obsdatabins %>% pull(ID),
      simdatabins %>% filter(!!REPL == min(!!REPL)) %>% pull(ID))) {
      stop("Error: ID's of your Obs and Sim data are not identical")
    }
  }
  
  if (!all.equal(
    obsdatabins %>% mutate(TIME = !!TIME) %>% pull(TIME),
    simdatabins %>% filter(!!REPL == min(!!REPL)) %>% mutate(TIME = !!TIME) %>% pull(TIME),
    check.attributes = FALSE)) {
    stop("Error: Time columns of Obs and Sim data are not identical or not in the same order")
  }
  
  stratifyvars <- all.vars(stratify)
  
  if (!quo_is_null(LLOQ)) {
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
  
  if (quo_is_null(NBINS) & is.null(breaks)) {
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
          left_join(obsdatabins %>% do(cutbreaks(., breaks)), #join to preserver original order
                    by=names(obsdatabins)) 
      }
      
    } else { 
      if (!is.null(style) && style == "ntile") {
        obsdatabins <- obsdatabins %>%
          mutate(BIN = ntile(!!TIME, !!NBINS))
      } else if (!is.null(style) && style %in% c("fixed", "sd", "equal", "pretty", "quantile", "kmeans", "hclust", "bclust", "fisher", "jenks")) {
        obsdatabins <- obsdatabins %>%
          mutate(BIN = as.numeric(cut(!!TIME, classIntervals(!!TIME, n=!!NBINS, style=style)$brks, right=cut_right, include.lowest=T)))
        
        breaks <- obsdatabins %>%
          mutate(TIME = !!TIME,
                 NBINS = !!NBINS) %>%
          select_at(c(stratifyvars, "TIME","NBINS")) %>%
          do(data.frame(breaks = classIntervals(.$TIME, n=unique(.$NBINS), style=style)$brks))
      } else {
        stop("Error: Unknown style")
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
      distinct(BIN, XMIN, XMAX) %>%
      arrange_at(c(group_vars(.),"BIN")) %>%
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
  
  if (quo_is_null(LLOQ)) {
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
        uv <- t(unique(.[,stratifyvarsbreaks]))
        msg1 <- paste0(paste(apply(uv, 2, function(x) paste(rownames(uv), x, sep="=")), collapse=", "), " ")
      } else {
        msg1 <- ""
      }
      msg2 <- paste0("[",length(.$XRIGHT),"] ")
      msg3 <- paste(signif(sort(c(min(.$XLEFT), .$XRIGHT)),3), collapse=" < ")
      message(paste0(msg1, msg2, msg3))
      
      data.frame()
    })
  
  simdatabins$BIN <- rep(obsdatabins$BIN, time = NREP)
  
  if (filterblq & !quo_is_null(LLOQ)) {
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
  
  if (quo_is_null(LLOQ)) {
    PIobs <- as.data.table(obsdatabins)[, 
                                        .(DVQNAME=paste0("",100*PI,"%PI"),
                                          DVOBS=quantile(DVC, probs = PI, na.rm=T)),
                                        by = c(stratifyvars,"BIN")]    
    
  } else {
    PIobs <- as.data.table(obsdatabins %>% 
                             mutate(LLOQ = !!LLOQ))[, #Don't know how to use !! inside data.table  
                                                    .(DVQNAME=paste0("",100*PI,"%PI"),
                                                      DVOBS=as.numeric(quantile_cens(DVC, p = PI, limit = LLOQ, na.rm=T))),
                                                    by = c(stratifyvars,"BIN")]
    PCTBLQobs <- as.data.table(obsdatabins)[, 
                                            .(PCTBLQQNAME = "PercentBLQ",
                                              PCTBLQOBS = 100 * mean(LLOQFL)),
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
                                      .(DVQNAME=paste0("",100*PI,"%PI"),
                                        DVSIM=quantile(DVC, probs = PI, na.rm=T)),
                                      by = c(stratifyvars,"BIN",quo_name(REPL))]
  
  PI <- PIsim[, 
              .(CI=paste0("DV",100*CI,"%CI"),
                Q=quantile(DVSIM, probs = CI, na.rm=T)),
              by = c(stratifyvars,"BIN","DVQNAME")]
  
  PI <- dcast(PI, as.formula(paste0(paste(c(stratifyvars, "BIN", "DVQNAME"), collapse="+"), "~CI")), value.var = "Q")
  
  PI <- merge(PI, PIobs, all.x=TRUE, by=c(stratifyvars,"BIN","DVQNAME"))
  
  if (!quo_is_null(LLOQ)) {
    PCTBLQsim <- as.data.table(simdatabins)[, 
                                            .(PCTBLQQNAME = "PercentBLQ",
                                              PCTBLQSIM = 100 * mean(LLOQFL)),
                                            by = c(stratifyvars,"BIN",quo_name(REPL))]
    
    PCTBLQ <- PCTBLQsim[, 
                        .(CI=paste0("PCTBLQ",100*CI,"%CI"),
                          Q=quantile(PCTBLQSIM, probs = CI, na.rm=T)),
                        by = c(stratifyvars,"BIN","PCTBLQQNAME")]
    
    PCTBLQ <- dcast(PCTBLQ, as.formula(paste0(paste(c(stratifyvars, "BIN", "PCTBLQQNAME"), collapse="+"), "~CI")), value.var = "Q")
    
    PCTBLQ <- merge(PCTBLQ, PCTBLQobs, all.x=TRUE, by=c(stratifyvars,"BIN","PCTBLQQNAME"))
    
  } else {
    PCTBLQ <- NULL 
  }
  
  # Return a list, as in pi_plot::gg_pi or vpc(..., vpcdb=T)
  # Advantage, the "DV" and "PCTBLQ" prefixes could be removed! Even PCTBLQQNAME is not very useful (kept for compatibility)
  #list(PI, PCTBLQ, BINS)
  
  VPC <- merge(PI,
                   if (!quo_is_null(LLOQ)) merge(PCTBLQ, BINS, all.x=TRUE, by=c(stratifyvars,"BIN")) else BINS,
                   all.x=TRUE, by=c(stratifyvars,"BIN"))
  
  if (sort_output) {
    VPC <- VPC %>% mutate(DVQNAMEN=as.numeric(gsub("%PI","",DVQNAME))) %>% arrange_at(c(stratifyvars,"BIN", "DVQNAMEN")) %>% select(-DVQNAMEN)
  }
  
  as.data.frame(VPC)
}

