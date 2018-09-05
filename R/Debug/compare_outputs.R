compare_outputs = function(df1, df2) {
  arrangevars = names(df2)[1:which(names(df2)=="DVPI")]

  df1 = df1 %>%
    ungroup() %>%
    rename(BIN=BINS,
           DVPI=quantilename,
           `DV2.5%CI`=QELCI,
           `DV50%CI`=QEMCI,
           `DV97.5%CI`=QEUCI,
           DVOBS=quantilevalue) %>%
    arrange_at(arrangevars)
  
  if ("percentblqQELCI" %in% names(df1)) df1 = df1 %>% rename(`PCTBLQ2.5%CI`=percentblqQELCI,
                                                              `PCTBLQ50%CI`=percentblqQEMCI,
                                                              `PCTBLQ97.5%CI`=percentblqQEUCI, 
                                                              PCTBLQOBS=percentblq)
  
  df1$DVPI=ifelse(df1$DVPI=="PLPI", "5%", ifelse(df1$DVPI=="PMPI", "50%", ifelse(df1$DVPI=="PUPI", "95%", NA)))
  
  df2 = df2 %>%
    arrange_at(arrangevars)
  
  sapply(intersect(names(df2),names(df1)), function(n) all.equal(df1 %>% pull(n), df2 %>% pull(n)))
}