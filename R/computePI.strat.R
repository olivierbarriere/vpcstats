require(ggplot2)
require(dplyr)
require(rlang)
require(vpc) # to compare test and augment on it

quantile_cens <- function(x, p = 0.5, limit = 1, cens = "left") {
  if (cens %in% c("left", "lower", "bloq", "loq", "lloq")) {
    x[is.na(x)] <- -Inf
    x[x < limit] <- -Inf
  }
  else {
    x[is.na(x)] <- Inf
    x[x > limit] <- Inf
  }
  q <- quantile(x, p)
  ifelse(q %in% c(Inf, -Inf), NA, q)
}

compute.PI.bins <- function(data = exampleobs,
                              simdata = examplesim,
                              stratify1 = NULL,
                              stratify2 = NULL,
                              stratify3 = NULL,
                              TIME = TIME,
                              DV = DV,
                              REPL = REP,
                              nbins = NULL,
                              LLOQ = NULL,
                              logadditive = FALSE,
                              predcorrection = FALSE,
                              PI = c(0.05, 0.5, 0.95),
                              CIPI = c(0.025, 0.5, 0.975),
                              bootstrapobsdata=FALSE) {
  stratify1 <- enquo(stratify1)
  if (quo_is_null(stratify1)) {
    print("stratify1 is null")
  }
  stratify2 <- enquo(stratify2)
  if (quo_is_null(stratify2)) {
    print("stratify2 is null")
  }
  stratify3 <- enquo(stratify3)
  if (quo_is_null(stratify3)) {
    print("stratify3 is null")
  }

  nbins <- enquo(nbins)
  if (quo_is_null(nbins)) {
    print("nbins is null no binning done")
  }

  TIME <- enquo(TIME)
  DV <- enquo(DV)
  REPL <- enquo(REPL)
  LLOQ <- enquo(LLOQ)
  if (quo_is_null(LLOQ)) {
    print("LLOQ is null")
  }
  originaldatabins <- data
  simdatabins <- simdata


  NREP <- nrow(simdatabins) / nrow(originaldatabins)
  print(paste("detected", NREP, "simulations"))

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
              select(!!TIME)),
    as.vector(simdatabins %>%
              filter(!!REPL == 1) %>%
              select(!!TIME))
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
    originaldatabins <- originaldatabins %>%
      group_by(BINS) %>%
      mutate(XMIN = min(!!TIME), XMAX = max(!!TIME)) %>%
      mutate(XMED = median(!!TIME)) %>%
      mutate(XMID = median(c(XMIN, XMAX))) %>%
      mutate(BINSC = paste("(", XMIN, "-", XMAX, "]", sep = ""))
  }
  if (!quo_is_null(nbins)) {
    if (!quo_is_null(stratify1)) {
      originaldatabins <- originaldatabins %>%
        group_by(!!stratify1, add = TRUE)
    }
    if (!quo_is_null(stratify2)) {
      originaldatabins <- originaldatabins %>%
        group_by(!!stratify2, add = TRUE)
    }
    if (!quo_is_null(stratify3)) {
      originaldatabins <- originaldatabins %>%
        group_by(!!stratify3, add = TRUE)
    }

    originaldatabins <- originaldatabins %>%
      mutate(BINS = ntile(!!TIME, !!nbins))

    originaldatabins <- originaldatabins %>%
      group_by(BINS, add = TRUE) %>%
      mutate(XMIN = min(!!TIME), XMAX = max(!!TIME)) %>%
      mutate(XMED = median(!!TIME)) %>%
      mutate(XMID = median(c(XMIN, XMAX))) %>%
      mutate(BINSC = paste("(", XMIN, "-", XMAX, "]", sep = ""))
  }


  if (quo_is_null(LLOQ)) {
    percentblqobs <- originaldatabins %>%
      dplyr::mutate(percentblq = NA)
    originaldatabins <- left_join(originaldatabins, percentblqobs)
  }
  if (!quo_is_null(LLOQ)) {
    percentblqobs <- originaldatabins %>%
      group_by(BINSC, XMIN, XMAX, XMED, XMID)

    if (!quo_is_null(stratify1)) {
      percentblqobs <- percentblqobs %>%
        group_by(!!stratify1, add = TRUE)
    }
    if (!quo_is_null(stratify2)) {
      percentblqobs <- percentblqobs %>%
        group_by(!!stratify2, add = TRUE)
    }
    if (!quo_is_null(stratify3)) {
      percentblqobs <- percentblqobs %>%
        group_by(!!stratify3, add = TRUE)
    }

    percentblqobs <- percentblqobs %>%
      dplyr::mutate(percentblq = 100 * mean(LLOQFL))
    originaldatabins <- left_join(originaldatabins, percentblqobs)
  }

  if (predcorrection) {
    originaldatabins <- originaldatabins %>%
      dplyr::group_by(BINSC)

    if (!quo_is_null(stratify1)) {
      originaldatabins <- originaldatabins %>%
        group_by(!!stratify1, add = TRUE)
    }
    if (!quo_is_null(stratify2)) {
      originaldatabins <- originaldatabins %>%
        group_by(!!stratify2, add = TRUE)
    }
    if (!quo_is_null(stratify3)) {
      originaldatabins <- originaldatabins %>%
        group_by(!!stratify3, add = TRUE)
    }

    originaldatabins <- originaldatabins %>%
      dplyr::mutate(mp = median(PRED))

    if (!logadditive) {
      originaldatabins <- originaldatabins %>%
        dplyr::mutate(dvc = !!DV * mp / PRED)
    }
    if (logadditive) {
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
      group_by(BINS, BINSC, percentblq, XMIN, XMAX, XMED, XMID)

    if (!quo_is_null(stratify1)) {
      PIobs <- PIobs %>%
        group_by(!!stratify1, add = TRUE)
    }
    if (!quo_is_null(stratify2)) {
      PIobs <- PIobs %>%
        group_by(!!stratify2, add = TRUE)
    }
    if (!quo_is_null(stratify3)) {
      PIobs <- PIobs %>%
        group_by(!!stratify3, add = TRUE)
    }
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
      group_by(BINS, BINSC, percentblq, XMIN, XMAX, XMED, XMID)

    if (!quo_is_null(stratify1)) {
      PIobs <- PIobs %>%
        group_by(!!stratify1, add = TRUE)
    }
    if (!quo_is_null(stratify2)) {
      PIobs <- PIobs %>%
        group_by(!!stratify2, add = TRUE)
    }
    if (!quo_is_null(stratify3)) {
      PIobs <- PIobs %>%
        group_by(!!stratify3, add = TRUE)
    }

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

  simdatabins$XMIN <- rep(originaldatabins$XMIN, time = NREP)
  simdatabins$XMED <- rep(originaldatabins$XMED, time = NREP)
  simdatabins$XMID <- rep(originaldatabins$XMID, time = NREP)
  simdatabins$XMAX <- rep(originaldatabins$XMAX, time = NREP)
  simdatabins$BINS <- rep(originaldatabins$BINS, time = NREP)
  simdatabins$BINSC <- rep(originaldatabins$BINSC, time = NREP)

  if (predcorrection) {
    simdatabins$mp <- rep(originaldatabins$mp, time = NREP)
    if (!logadditive) {
      simdatabins <- simdatabins %>%
        dplyr::mutate(dvc = !!DV * mp / PRED)
    }
    if (logadditive) {
      simdatabins <- simdatabins %>%
        dplyr::mutate(dvc = !!DV + (mp - PRED))
    }
  }

  if (!predcorrection) {
    simdatabins <- simdatabins %>%
      dplyr::mutate(dvc = !!DV)
  }


  PICIsim <- simdatabins %>%
    group_by(!!REPL, BINS, BINSC, XMIN, XMAX, XMED, XMID)

  if (!quo_is_null(stratify1)) {
    PICIsim <- PICIsim %>%
      group_by(!!stratify1, add = TRUE)
  }


  if (!quo_is_null(stratify2)) {
    PICIsim <- PICIsim %>%
      group_by(!!stratify2, add = TRUE)
  }
  if (!quo_is_null(stratify3)) {
    PICIsim <- PICIsim %>%
      group_by(!!stratify3, add = TRUE)
  }

  PICIsim <- PICIsim %>%
    dplyr::summarize(
      PLPI = quantile(dvc, probs = PI[1]),
      PMPI = quantile(dvc, probs = PI[2]),
      PUPI = quantile(dvc, probs = PI[3])
    )
  SIMPIPI <- tidyr::gather(PICIsim, quantilename, quantilevalue, PLPI, PMPI, PUPI)
  VPCSTAT <- SIMPIPI %>%
    group_by(quantilename, BINS, BINSC, XMIN, XMAX, XMED, XMID)
  if (!quo_is_null(stratify1)) {
    VPCSTAT <- VPCSTAT %>%
      group_by(!!stratify1, add = TRUE)
  }
  if (!quo_is_null(stratify2)) {
    VPCSTAT <- VPCSTAT %>%
      group_by(!!stratify2, add = TRUE)
  }
  if (!quo_is_null(stratify3)) {
    VPCSTAT <- VPCSTAT %>%
      group_by(!!stratify3, add = TRUE)
  }

  VPCSTAT <- VPCSTAT %>%
    dplyr::summarize(
      QELCI = quantile(quantilevalue, probs = CIPI[1]),
      QEMCI = quantile(quantilevalue, probs = CIPI[2]),
      QEUCI = quantile(quantilevalue, probs = CIPI[3])
    )


  if (!quo_is_null(LLOQ)) {
    percentblqsim <- simdatabins %>%
      group_by(!!REPL, BINSC, XMIN, XMAX, XMED, XMID)

    if (!quo_is_null(stratify1)) {
      percentblqsim <- percentblqsim %>%
        group_by(!!stratify1, add = TRUE)
    }
    if (!quo_is_null(stratify2)) {
      percentblqsim <- percentblqsim %>%
        group_by(!!stratify2, add = TRUE)
    }
    if (!quo_is_null(stratify3)) {
      percentblqsim <- percentblqsim %>%
        group_by(!!stratify3, add = TRUE)
    }

    percentblqsim <- percentblqsim %>%
      dplyr::summarize(percentblq = mean(LLOQFL))
    percentblqsim <- tidyr::gather(percentblqsim, blqquantilename, blqquantilevalue, percentblq)

    percentblqsimSTAT <- percentblqsim %>%
      group_by(blqquantilename, BINSC, XMIN, XMAX, XMED)
    if (!quo_is_null(stratify1)) {
      percentblqsimSTAT <- percentblqsimSTAT %>%
        group_by(!!stratify1, add = TRUE)
    }
    if (!quo_is_null(stratify2)) {
      percentblqsimSTAT <- percentblqsimSTAT %>%
        group_by(!!stratify2, add = TRUE)
    }
    if (!quo_is_null(stratify3)) {
      percentblqsimSTAT <- percentblqsimSTAT %>%
        group_by(!!stratify3, add = TRUE)
    }

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

  VPCSTAT
}

