is.formula <- function(x){
  inherits(x,"formula")
}
# fix grouping by all specified strata
vpcplotfunction <- function (data=vpcdara, facet_formula = "",facetscales="free_x",
                             stratvars="",
                             xlab="Time",ylab="VPC",shownobs=FALSE,showrug=TRUE,
                             pointrange=FALSE,pointrangenudge= 0,smoothlines=FALSE,rectangles=TRUE,
                             detachedrectangles=FALSE,
                             bin_mid= "XMID",
                             PItransparency=0.2,PIcolors =  rep(c("red", "blue"),10),
                             PIlinetypes= rep(c( "dashed", "solid"),10), obslinesizes=1,
                             log_y=FALSE,log_x=FALSE,exp_y=FALSE,log_y_min=-Inf,
                             facetwrap = FALSE,showLLOQ=FALSE,
                             legend.position =c("right")
                             ){
  
  # PIdata<- VPCDATA$PI
  # BINSdata<- VPCDATA$BINS ; PItransparency=0.2;PIcolors =  rep(c("red", "blue"),10);obslinesizes=1;bin_mid= "XMID"
  # need to  clearly tell the plot strat variables since facetting might be different
  PIdata<- data$PI
  BINSdata<- data$BINS
  shownobs <- shownobs
  QNAMES <- unique(PIdata$QNAME)
  fix_names <- function(x) gsub("%", "percent", x)
  names(PIdata)<- fix_names(  names(PIdata))
  CInames <- tidyselect::vars_select(names(PIdata), starts_with("SIM"))
  CInames<- as.vector(CInames)
  CILEVELS <-  readr::parse_number(CInames)
  CILEVEL <- max(CILEVELS) - min(CILEVELS)
  
  PIdata$QNAME <- factor(PIdata$QNAME)
  PIdata$QNAME <- reorder(  PIdata$QNAME,PIdata$RAWOBS,mean,na.rm=TRUE)
  QNAMES <- levels(PIdata$QNAME)
  
  if(exp_y){
    for (i in c(CInames,"RAWOBS")  ){
      PIdata[i] <- exp(PIdata[i])
          }
    
  }
  if (log_y) {
    for (i in c(CInames,"RAWOBS")  ){
      PIdata[i][PIdata[i] < log_y_min] <- log_y_min
    }
  }

  p <- ggplot(PIdata)
  if (smoothlines) {
    p <- p +
      geom_ribbon(
        aes_string(
          x = bin_mid,
          ymin =  CInames[1] ,
          ymax =  CInames[3],
          fill = "QNAME",
          group = "QNAME"
        ),
        alpha = PItransparency,
        linetype = 0
      )+
      geom_line(aes_string(
        x=bin_mid,
        y = CInames[2],
        col = "QNAME",
        group =  "QNAME"
      )) +
      geom_line(aes_string(
        x=bin_mid,
        y = "RAWOBS",
        group = "QNAME",
        linetype =  "QNAME"
      ), size = obslinesizes) 
  }

  
  if (rectangles) {
    if(!detachedrectangles){
      p <- p +
      geom_rect(
        aes_string(
          xmin="XLEFT",
          xmax="XRIGHT",
          ymin =  CInames[1] ,
          ymax =  CInames[3],
          fill = "QNAME",
          group = "QNAME"
        ),
        alpha = PItransparency,
        col = NA,
        inherit.aes = FALSE
      )+
      geom_segment(aes_string(
        x="XLEFT",
        xend="XRIGHT",
        y =  CInames[2],
        yend =  CInames[2],
        col = "QNAME",
        group = "QNAME"
      )) +
      geom_segment(aes_string(
        x="XLEFT",
        xend="XRIGHT",
        y = "RAWOBS",
        yend = "RAWOBS",
        group = "QNAME",
        linetype = "QNAME"
      ), size = obslinesizes) 
  }
    if(detachedrectangles){
      p <- p +
        geom_rect(
          aes_string(
            xmin="XMIN",
            xmax="XMAX",
            ymin =  CInames[1] ,
            ymax =  CInames[3],
            fill = "QNAME",
            group = "QNAME"
          ),
          alpha = PItransparency,
          col = NA,
          inherit.aes = FALSE
        )+
        geom_segment(aes_string(
          x="XMIN",
          xend="XMAX",
          y =  CInames[2],
          yend =  CInames[2],
          col =  "QNAME",
          group =  "QNAME"
        )) +
        geom_segment(aes_string(
          x="XMIN",
          xend="XMAX",
          y = "RAWOBS",
          yend = "RAWOBS",
          group = "QNAME",
          linetype = "QNAME"
        ), size = obslinesizes) 
    }
  }
  
     if (pointrange) {
    p <- p +
      geom_pointrange(aes_string(
        x=paste0(bin_mid,"-pointrangenudge"),
        y= CInames[2],
        ymin = CInames[1],
        ymax = CInames[3],
        col =  "QNAME"
      ))+
      geom_point(aes_string(
        x=paste0(bin_mid,"+pointrangenudge"),
        y= "RAWOBS",shape=shQuote("obs")
      ),col = "black")
     }
  p <- p +
    scale_colour_manual(
      name = paste0("Simulated\nMedian (lines) ",CILEVEL," % CI (areas)"),
      breaks = QNAMES,
      values = PIcolors,
      labels = gsub('.{2}$', '', QNAMES)
    ) +
    scale_fill_manual(
      name = paste0("Simulated\nMedian (lines) ",CILEVEL," % CI (areas)"),
      breaks = QNAMES,
      values = PIcolors,
      labels =  gsub('.{2}$', '', QNAMES)
    ) +
    scale_linetype_manual(
      name = "Observed (black lines)",
      breaks = QNAMES,
      values = PIlinetypes,
      labels =  gsub('.{2}$', '', QNAMES)
    ) +
    guides(
      fill = guide_legend(order = 2, reverse = TRUE),
      colour = guide_legend(order = 2, reverse = TRUE),
      linetype = guide_legend(order = 1, reverse = TRUE)
    ) +
    theme_bw(base_size = 18)+
    theme(
      legend.key.width = grid::unit(2, "cm"),
      axis.text.x = element_text(angle = 0))+
      labs(x=xlab,y=ylab, shape="")
  

  if (shownobs&& !log_y&& !facetwrap) {
    p <- p+ geom_text_repel(data=BINSdata,
      aes(x = XMID,y=-Inf,label=NOBS), direction="x" )
  } 
  
  if (shownobs&& !log_y&& facetwrap) {
    p <- p+ geom_text(data=BINSdata,
                            aes(x = XMID,y=0 ,label=NOBS))
  } 
  
  if (showrug) {
    p <-  p+
      geom_rug(data=BINSdata, aes(x=XLEFT ), size =1,sides = "t")+
      geom_rug(data=BINSdata, aes(x=XRIGHT ),size =1,sides = "t")
  }
  if (showLLOQ) {
  p <-  p+
    geom_hline(data = BINSdata,
               aes_string(
                 yintercept = "LLOQ"
               ), linetype = "dotted", size = 1) 
  
  }
  
  if (log_y) {
    p <-  p+ 
      scale_y_log10()
    if (shownobs) {
      p <- p+ geom_text_repel(data=BINSdata,
                              aes(x=XMED,y=0,label=NOBS), direction="x",vjust=1 )
    }
  } 
  
  if (exp_y & !log_y) {
    if (shownobs) {
      p <- p+ geom_text_repel(data=BINSdata,
                              aes(x=XMED,y=1,label=NOBS), direction="x",vjust=1 )
    }
  } 
  
  if (log_x) {
    p <-  p+ 
      scale_x_log10()
    
  } 
  if ( facet_formula  != "")  {
    facet_formula <- stats::as.formula(facet_formula)
    if(is.formula(facet_formula)&&!facetwrap)    p <- p+ facet_grid(facet_formula,scales=facetscales,as.table = FALSE)
    if(is.formula(facet_formula)&&facetwrap)    p <- p+ facet_wrap(facet_formula,scales=facetscales,as.table = FALSE)
    
  }
  p <- p +theme(legend.position = legend.position)

  p
  
    }


pctblqplotfunction <- function (data=vpcdara, facet_formula = "",facetscales="free_x",
                                stratvars="",
                                xlab="Time",ylab="VPC",shownobs=FALSE,showrug=TRUE,
                                pointrange=FALSE,pointrangenudge= 0,smoothlines=FALSE,rectangles=TRUE,
                                detachedrectangles=FALSE,
                                bin_mid= "XMID",
                                PItransparency=0.2,PIcolors =  rep(c("red", "blue"),10),
                                PIlinetypes= rep(c( "dashed", "solid"),10),obslinesizes=1,
                                log_y=FALSE,log_x=FALSE,
                                facetwrap = FALSE,
                                legend.position =c("right")
                                
){
  PIdata<- data$PI
  BINSdata<- data$BINS
  PCTBLQdata<- data$PCTBLQ

  shownobs <- shownobs
  QNAMES <- unique(PIdata$QNAME)
  fix_names <- function(x) gsub("%", "percent", x)
  names(PIdata)<- fix_names(  names(PIdata))
  names(PCTBLQdata)<- fix_names(  names(PCTBLQdata))
  
  CInames <- tidyselect::vars_select(names(PIdata), starts_with("SIM"))
  CInames<- as.vector(CInames)
  CILEVELS <-  readr::parse_number(CInames)
  CILEVEL <- max(CILEVELS) - min(CILEVELS)
  
  PIdata$QNAME <- factor(PIdata$QNAME)
  PIdata$QNAME <- reorder(  PIdata$QNAME,PIdata$RAWOBS,mean,na.rm=TRUE)
  PIdata$QNAME <- factor(PCTBLQdata$QNAME)
  PCTBLQdata$QNAME<- "BLQ"
  QNAMES <- levels(factor(PCTBLQdata$QNAME))
  p <- ggplot(PCTBLQdata)
  if (smoothlines) {
    p <- p +
      geom_ribbon(
        aes_string(
          x = bin_mid,
          ymin =  CInames[1] ,
          ymax =  CInames[3],
          fill = "QNAME",
          group = "QNAME"
        ),
        alpha = PItransparency,
        linetype = 0
      )+
      geom_line(aes_string(
        x=bin_mid,
        y = CInames[2],
        col = "QNAME",
        group =  "QNAME"
      )) +
      geom_line(aes_string(
        x=bin_mid,
        y = "RAWOBS",
        group = "QNAME",
        linetype =  "QNAME"
      ), size = obslinesizes) 
  }

  if (rectangles) {
    if(!detachedrectangles){
      p <- p +
        geom_rect(
          aes_string(
            xmin="XLEFT",
            xmax="XRIGHT",
            ymin =  CInames[1] ,
            ymax =  CInames[3],
            fill = "QNAME",
            group = "QNAME"
          ),
          alpha = PItransparency,
          col = NA,
          inherit.aes = FALSE
        )+
        geom_segment(aes_string(
          x="XLEFT",
          xend="XRIGHT",
          y =  CInames[2],
          yend =  CInames[2],
          col = "QNAME",
          group = "QNAME"
        )) +
        geom_segment(aes_string(
          x="XLEFT",
          xend="XRIGHT",
          y = "RAWOBS",
          yend = "RAWOBS",
          group = "QNAME",
          linetype = "QNAME"
        ), size = obslinesizes) 
    }
    if(detachedrectangles){
      p <- p +
        geom_rect(
          aes_string(
            xmin="XMIN",
            xmax="XMAX",
            ymin =  CInames[1] ,
            ymax =  CInames[3],
            fill = "QNAME",
            group = "QNAME"
          ),
          alpha = PItransparency,
          col = NA,
          inherit.aes = FALSE
        )+
        geom_segment(aes_string(
          x="XMIN",
          xend="XMAX",
          y =  CInames[2],
          yend =  CInames[2],
          col =  "QNAME",
          group =  "QNAME"
        )) +
        geom_segment(aes_string(
          x="XMIN",
          xend="XMAX",
          y = "RAWOBS",
          yend = "RAWOBS",
          group = "QNAME",
          linetype = "QNAME"
        ), size = obslinesizes) 
    }
  }
  
  if (pointrange) {
    p <- p +
      geom_pointrange(aes_string(
        x=paste0(bin_mid,"-pointrangenudge"),
        y= CInames[2],
        ymin = CInames[1],
        ymax = CInames[3],
        col =  "QNAME"
      ))+
      geom_point(aes_string(
        x=paste0(bin_mid,"+pointrangenudge"),
        y= "RAWOBS",shape=shQuote("obs")
      ),col = "black")
  }
  p <- p +
    scale_colour_manual(
      name = paste0("Simulated\nMedian (lines) ",CILEVEL," % CI (areas)"),
      breaks = QNAMES,
      values = PIcolors,
      labels = QNAMES
    ) +
    scale_fill_manual(
      name = paste0("Simulated\nMedian (lines) ",CILEVEL," % CI (areas)"),
      breaks = QNAMES,
      values = PIcolors,
      labels =  QNAMES
    ) +
    scale_linetype_manual(
      name = "Observed (black lines)",
      breaks = QNAMES,
      values = PIlinetypes,
      labels =  QNAMES
    ) +
    guides(
      fill = guide_legend(order = 2, reverse = TRUE),
      colour = guide_legend(order = 2, reverse = TRUE),
      linetype = guide_legend(order = 1, reverse = TRUE)
    ) +
    theme_bw(base_size = 18)+
    theme(
      legend.key.width = grid::unit(2, "cm"),
      axis.text.x = element_text(angle = 0))+
    labs(x=xlab,y=ylab, shape="")


  if (shownobs&& !log_y&& !facetwrap) {
    p <- p+ geom_text_repel(data=BINSdata,
                            aes(x = XMID,y=-Inf,label=NOBS), direction="x" )
  } 
  
  if (shownobs&& !log_y&& facetwrap) {
    p <- p+ geom_text(data=BINSdata,
                      aes(x = XMID,y=0 ,label=NOBS))
  } 
  
  if (showrug) {
    p <-  p+
      geom_rug(data=BINSdata, aes(x=XLEFT ), size =1,sides = "t")+
      geom_rug(data=BINSdata, aes(x=XRIGHT ),size =1,sides = "t")
  }
  if (log_y) {
    p <-  p+ 
      scale_y_log10()
    if (shownobs) {
      p <- p+ geom_text_repel(data=BINSdata,
                              aes(x=XMED,y=0,label=NOBS), direction="x",vjust=1 )
    }
  } 
  if (log_x) {
    p <-  p+ 
      scale_x_log10()
    
  } 
  if ( facet_formula  != "")  {
    facet_formula <- stats::as.formula(facet_formula)
    if(is.formula(facet_formula)&&!facetwrap)    p <- p+ facet_grid(facet_formula,scales=facetscales,as.table = FALSE)
    if(is.formula(facet_formula)&&facetwrap)    p <- p+ facet_wrap(facet_formula,scales=facetscales,as.table = FALSE)
    
  }
  p <- p + theme(legend.position = legend.position)
  
  p
}