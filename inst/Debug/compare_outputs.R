compare_outputs = function(df1, listdf2) {

  df2PI=listdf2$PI %>% rename(DVQNAME=QNAME, DVOBS=OBS)
  names(df2PI)[grepl("%CI$",names(df2PI))] = paste0("DV",names(df2PI)[grepl("%CI$",names(df2PI))])
  
  df2PCTBLQ=listdf2$PCTBLQ
  if (!is.null(df2PCTBLQ)) {
    df2PCTBLQ=df2PCTBLQ %>% rename(PCTBLQOBS=OBS)
    names(df2PCTBLQ)[grepl("%CI$",names(df2PCTBLQ))] = paste0("PCTBLQ",names(df2PCTBLQ)[grepl("%CI$",names(df2PCTBLQ))])
  }
  
  df2=left_join(df2PI,df2PCTBLQ, by=intersect(names(df2PI), names(df2PCTBLQ)))
  
  arrangevars = names(df2)[1:which(names(df2)=="DVQNAME")]
  
  df1 = df1 %>%
    ungroup() %>%
    rename(BIN=BINS,
           DVQNAME=quantilename,
           `DV2.5%CI`=QELCI,
           `DV50%CI`=QEMCI,
           `DV97.5%CI`=QEUCI,
           DVOBS=quantilevalue) %>%
    arrange_at(arrangevars)
  
  df1$DVQNAME=ifelse(df1$DVQNAME=="PLPI", "5%PI", ifelse(df1$DVQNAME=="PMPI", "50%PI", ifelse(df1$DVQNAME=="PUPI", "95%PI", NA)))
  
  if ("blqquantilename" %in% names(df1)) {
    df1 = df1 %>% rename(PCTBLQQNAME=blqquantilename,
                         `PCTBLQ2.5%CI`=percentblqQELCI,
                         `PCTBLQ50%CI`=percentblqQEMCI,
                         `PCTBLQ97.5%CI`=percentblqQEUCI, 
                         PCTBLQOBS=percentblq)
    df1$PCTBLQQNAME=ifelse(df1$PCTBLQQNAME=="percentblq", "PercentBLQ", NA)
  }  
  
  df2 = df2 %>%
    arrange_at(arrangevars)
  
  sapply(intersect(names(df2),names(df1)), function(n) all.equal(df1 %>% pull(n), df2 %>% pull(n)))
}