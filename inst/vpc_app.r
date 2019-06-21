vpc_app <-  function(obsdata = obsdata, simdata = simdata, stratalist = stratalist ) {
  
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(shiny) 
  library(scales)  
  library(colourpicker)  
  library(plotly)
  library(shinyjs)
  shinyApp(
    ui = fluidPage(
      useShinyjs(),
      fluidRow(
        plotOutput('vpcplot', height = "100%")  
      ),
      fluidRow(style = "padding-bottom: 0px;",
               column(2,
                      textInput("PI",label = "Prediction Intervals:",value= c("0.05,0.5,0.95")),
                      sliderInput("CI", "Confidence Intervals:",min = 0, max = 0.99, value=c(0.95),step=0.01),
                      sliderInput("quantile_type", "Quantile Method:",min=6, max=7, 7,step=1),
                      checkboxInput("predcorrection",label = "Predcorrection",value=FALSE),
                      conditionalPanel(condition = 'input.predcorrection',
                      checkboxInput("predcorrection_islogdv",label = "was dv log transformed ?",value=FALSE) ),
                      radioButtons("TIME", "Time Variable:",c("TIME" = "TIME","TAD" = "TAD"), inline = TRUE)
                      
               )
               ,
               column(2,  radioButtons("bin_style", "Binning Style:",
                                       c("ntile" = "ntile",
                                         "sd" = "sd",
                                         "equal" = "equal",
                                         "pretty" = "pretty",
                                         "quantile" = "quantile",
                                         "kmeans" = "kmeans",
                                         "jenks" = "jenks",
                                         "breaks" = "breaks"),
                                       
               ),
               checkboxInput("nobinning",label = "No Binning",value=FALSE),
               checkboxInput("bin_by_strata",label = "Different Binning by strata ?",value=TRUE),
                numericInput('nbins', 'N Bins:', 4, min = 2, max = 20),
                uiOutput('breaksstratanostrata')
               ),
               column(2,
                      radioButtons("LLOQ", "LLOQ:",c("None" = "None","LLOQ" = "LLOQ"), inline = TRUE),
                      checkboxInput("filterblq",label = "Filter <LLOQ ?",value=FALSE),
                      conditionalPanel(condition = 'input.LLOQ!="None" ', 
                                       checkboxInput("plotblq",label = "Plot % BLQ ?",value=FALSE)
                                       )
                      
               ),
               column(2,
                      selectInput('vpcformula' ,'VPC stratify variables:',c(stratalist ), multiple = TRUE,selected= stratalist[1] ),
                      uiOutput('checkboxesforselectedstratvars') 
                      
               ),
               column(2,
                      checkboxInput("showobspoints",label = "Show Obs Points ?",value=TRUE),
                      checkboxInput("rectangles",label = "Rectangles:",value=FALSE),
                      conditionalPanel(condition = 'input.rectangles',
                                       checkboxInput("detachedrectangle",label = "Disconnected Rect?",value=FALSE) ),
                      checkboxInput("smoothlines",label = "Smooth Lines ?",value=TRUE),
                      checkboxInput("pointrange",label = "Point Range ?",value=FALSE),
                      conditionalPanel(condition = 'input.smoothlines|input.pointrange',
                                       radioButtons("bin_mid", "Midpoint Style:",
                                                    c("XMID" = "XMID",
                                                      "XMED" = "XMED",
                                                      "XMEAN" = "XMEAN",
                                                      "XCENTER" = "XCENTER")) 
                      ),
                      checkboxInput("showrug",label = "Show Bin Delimeters ?",value=TRUE),
                      checkboxInput("shownobs",label = "Show N OBS ?",value=TRUE),
                      checkboxInput("log_y",label = "log_y",value=FALSE),
                      conditionalPanel(condition = 'input.log_y ',
                      checkboxInput("log_y_min_set",label = "Specify Min Value for y on log scale ?",value=FALSE)),
                      conditionalPanel(condition = 'input.log_y&&input.log_y_min_set ',
                                       numericInput('log_y_min', 'Min Value for Y with logged scale', 0.0001)),
                      checkboxInput("log_x",label = "log_x",value=FALSE),
                      checkboxInput("exp_y",label = "Exponentiate y first",value=FALSE),
                      checkboxInput("legend.position",label = "Hide Legend ?",value=FALSE)
                      
                      
               ),
               column(2,
                      uiOutput("facetvarsrows"),uiOutput("facetvarscols"),
                      uiOutput("facetvarsrows2"),uiOutput("facetvarscols2"),
                      selectInput('facetscales' ,'Facet Scales:',c("fixed","free_x","free_y","free"),selected =c("free_y") ),
                      checkboxInput("facetwrap", label = "Facet Wrap ?", value = FALSE) 
               )
               
      ),
      fluidRow(style = "padding-bottom: 0px;",
               column(3,
                      sliderInput("PItransparency", "PI Opacity:", min = 0, max = 1, value = c(0.2), step=0.01),
                      checkboxInput("custompicolors", label = "User Colors for PIs ?", value = FALSE) ,
                      conditionalPanel(condition = 'input.custompicolors', uiOutput('userdefinedcolor'))
               ),
               column(3,
                      sliderInput("obslinesizes", "Obs line size:", min = 0, max = 3, value = c(1), step = 0.025),
                      checkboxInput("custompilines", label = "User Lines for PIs ?", value = FALSE) ,
                      conditionalPanel(condition = 'input.custompilines', uiOutput('userdefinedlinetypes'))
               ),
               column(3,
                      sliderInput("blqtransparency", "Blq Opacity:", min = 0, max = 1, value = c(0.2), step=0.01),
                      checkboxInput("customblqcolors", label = "User Colors for %BLQ ?", value = FALSE) ,
                      conditionalPanel(condition = 'input.customblqcolors', uiOutput('userdefinedblqcolor'))
               ),
               column(3,
                      sliderInput("obsblqlinesize", "Obs Blq line size:", min = 0, max = 3, value = c(1), step = 0.025),
                      checkboxInput("customblqlines", label = "User Lines for %BLQ ?", value = FALSE) ,
                      conditionalPanel(condition = 'input.customblqlines', uiOutput('userdefinedblqlinetype'))
               )
      )#fluidrow
    ),#fluidpage ui start
    
    server = function(input, output, session) {
      
      tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
                     "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
      
      tableau20 <- c("#1F77B4","#AEC7E8", "#FF7F0E","#FFBB78"  ,"#2CA02C",
                     "#98DF8A" ,"#D62728","#FF9896" ,"#9467BD","#C5B0D5" ,
                     "#8C564B","#C49C94" ,"#E377C2","#F7B6D2" ,"#7F7F7F",
                     "#C7C7C7" ,"#BCBD22","#DBDB8D" ,"#17BECF","#9EDAE5")
      
      PIcolors <- rep(c("red", "blue"),10)
      PIlinetypes <- rep(c( "dashed", "solid"),10)
      blqcolor <-"red"
      blqlinetype <- "dashed"

      output$facetvarsrows <- renderUI({
                        selectizeInput("facetrows",label="Facet Rows",
                         choices = c(None='.',"QNAME") ,
                         multiple = FALSE,
                         selected = "QNAME") })
      
      output$facetvarsrows2 <- renderUI({
        selectizeInput("facetrows2",label="Facet Rows(2)",
                       choices = c(None='.',"QNAME") ,
                       multiple = FALSE,
                       selected = "QNAME") })
      
      output$facetvarscols <- renderUI({
                        selectizeInput("facetcols",label="Facet Cols",
                       choices = c(None='.',"QNAME") ,
                       multiple = FALSE,
                       selected = ".") })
      
      output$facetvarscols2 <- renderUI({
        selectizeInput("facetcols2",label="Facet Cols(2)",
                       choices = c(None='.',"QNAME") ,
                       multiple = FALSE,
                       selected = ".") })
      
      
      facetformulaconstructed = reactive({
        req(input$facetrows)
        req(input$facetrows2)
        req(input$facetcols)
        req(input$facetcols2)
        allfacetsvariables<- c(input$facetrows,input$facetrows2,input$facetcols,input$facetcols2)
        allfacetsvariables[which(duplicated(allfacetsvariables))]<- "." # make it not fail
        
        facets <- paste(allfacetsvariables[1],'+',allfacetsvariables[2],'~',allfacetsvariables[3],'+',allfacetsvariables[4])
        facets <-ifelse(facets%in%c('. + . ~ . + .'),"",facets)
        facets
        
      })
      vpcformulaupdate = reactive({
        if(length(input$vpcformula) == 0)   vpcformula <- NULL
        if(length(input$vpcformula)>= 1)  vpcformula <- as.formula(paste(" ~ ", paste(c(input$vpcformula), collapse = " + "), sep = " "))
        if(!is.formula(vpcformula)) vpcformula <- ""
        vpcformula
      })

 

      
      
      observe({
        vpcformulaupdate()
        stratvars<- all.vars( vpcformulaupdate() )
        choicesupdate<- c( c(None='.',"QNAME") , stratvars)
        updateSelectInput(session, "facetcols", choices = choicesupdate,
                          selected = ifelse(length(input$vpcformula)>= 1,stratvars[1],'.' ) )
      })

      observe({
        vpcformulaupdate()
        stratvars<- all.vars( vpcformulaupdate() )
        choicesupdate<- c( c(None='.',"QNAME") , stratvars)
        updateSelectInput(session, "facetrows", choices = choicesupdate,selected = "QNAME")
      })
      observe({
        vpcformulaupdate()
        stratvars<- all.vars( vpcformulaupdate() )
        choicesupdate<- c( c(None='.',"QNAME") , stratvars)
        updateSelectInput(session, "facetrows2", choices = choicesupdate,selected = '.')
      })
      
      output$checkboxesforselectedstratvars <- renderUI({ 
        req(input$vpcformula)
        if(length(input$vpcformula) == 0)  return()
        if(length(input$vpcformula) >= 1) {
          lev <- 1:length(input$vpcformula)
          lapply(seq_along(lev), function(i) {
            checkboxGroupInput(inputId = paste0('checkboxgrp', lev[i]),
                               label = paste0("What to Include in Plot for ", input$vpcformula[i]),
                               choices =unique(levels(as.factor( obsdata[,input$vpcformula[i]] )) ),
                               selected = unique(obsdata[,input$vpcformula[i]] ), inline = TRUE
            )
            
            
            
          })
          
        }
        
      })
      
      
      output$breaksstratanostrata <- renderUI({ 
          textInput(inputId = paste0('userbreakvaluesnostrata'),
                    label = paste0("What breaks to use:"),
                    value = c("0,10,100"))
      })

      
      observe({
        if (input$bin_style!="breaks" ) {
          shinyjs::enable("bin_by_strata")
        }
        if (input$bin_style=="breaks" ) {
          updateCheckboxInput(session, "bin_by_strata", value = FALSE)
          shinyjs::disable("bin_by_strata")
        }
      })
      
      
      vpcData = reactive({
        if (input$LLOQ == "None") LLOQ<- NULL
        if (input$LLOQ != "None") LLOQ<- rlang::sym(input$LLOQ)

        if(input$bin_style!="breaks"&&!input$nobinning) {
          vpcpercentiles<- vpcstats(
            obsdata = obsdata,
            simdata = simdata,
            stratify = vpcformulaupdate(),
            TIME = !!rlang::sym(input$TIME),
            REPL = REPLICATE,
            LLOQ = !!LLOQ,
            NBINS = input$nbins,
            bin_by_strata = input$bin_by_strata,
            bin_style = input$bin_style,
            breaks = NULL ,
            cut_right = FALSE,
            quantile_type = input$quantile_type,
            filterblq = input$filterblq,
            predcorrection = input$predcorrection,
            predcorrection_islogdv = input$predcorrection_islogdv,
            predcorrection_lowerbnd = 0,
            PI = as.numeric(unlist(strsplit(input$PI, split=","))),
            CI = c( (1-input$CI)/2,0.5  , ( 1-  (1-input$CI)/2 ) )
          )
          
        }
        if(input$bin_style=="breaks"&&!input$nobinning) {
          vpcpercentiles<- vpcstats(
            obsdata = obsdata,
            simdata = simdata,
            stratify = vpcformulaupdate(),
            TIME = !!rlang::sym(input$TIME),
            REPL = REPLICATE,
            LLOQ = !!LLOQ,
            NBINS = NULL,
            bin_by_strata = input$bin_by_strata,
            bin_style = input$bin_style,
            breaks=  as.numeric(unlist(strsplit(input$userbreakvaluesnostrata, split=","))) ,
            cut_right = FALSE,
            quantile_type = input$quantile_type,
            filterblq = input$filterblq,
            predcorrection = input$predcorrection,
            predcorrection_islogdv = input$predcorrection_islogdv,
            predcorrection_lowerbnd = 0,
            PI = as.numeric(unlist(strsplit(input$PI, split=","))),
            CI = c( (1-input$CI)/2,0.5  , ( 1-  (1-input$CI)/2 ) )
          )
        }
        if(input$nobinning) {
          vpcpercentiles<- vpcstats(
            obsdata = obsdata,
            simdata = simdata,
            stratify = vpcformulaupdate(),
            TIME = !!rlang::sym(input$TIME),
            REPL = REPLICATE,
            LLOQ = !!LLOQ,
            NBINS = NULL,
            bin_by_strata = input$bin_by_strata,
            bin_style = NULL,
            breaks=  NULL,
            cut_right = FALSE,
            quantile_type = input$quantile_type,
            filterblq = input$filterblq,
            predcorrection = input$predcorrection,
            predcorrection_islogdv = input$predcorrection_islogdv,
            predcorrection_lowerbnd = 0,
            PI = as.numeric(unlist(strsplit(input$PI, split=","))),
            CI = c( (1-input$CI)/2,0.5  , ( 1-  (1-input$CI)/2 ) )
          )
        }
        
        
        vpcpercentiles
      })
      
      output$userdefinedcolor <- renderUI({ 
        req(input$PI)
        lev <- 1:length(as.numeric(unlist(strsplit(input$PI, split=","))))
        if( length(lev) > 10 ){
          cols <- c(tableau20)
        }
        if( length(lev) <= 10 ){
          cols <- c(tableau10)
        }
        lapply(seq_along(lev), function(i) {
          colourpicker::colourInput(inputId = paste0("col", lev[i]),label = paste0("Color for PI", lev[i]), value = cols[i])        
        })
      })
      
      output$userdefinedlinetypes <- renderUI({ 
        req(input$PI)
        lev <- 1:length(as.numeric(unlist(strsplit(input$PI, split=","))))
        linetypechoices <-  c("solid","dashed", "dotted", "dotdash", "longdash", "twodash","blank","solid","dashed", "dotted")
        lapply(seq_along(lev), function(i) {
          selectizeInput(
            inputId = paste0('linetypes', lev[i]),
            label = paste0("Lines for PI", lev[i]),
            
            choices = linetypechoices,
            options = list(
              placeholder = 'Please select an option below',
              onInitialize = I('function() { this.setValue("solid"); }')
            )
          )
        })
        
      })
      
      output$userdefinedblqcolor <- renderUI({ 
        if(!input$plotblq) {
          return()        }
        if(input$plotblq) {
          colourpicker::colourInput(inputId = paste0("col","blq"),label = paste0("Color for %BLQ"), value = tableau20[6]) 
        }       
        
      }) 
      output$userdefinedblqlinetype <- renderUI({ 
        if(!input$plotblq) {
          return()        }
        if(input$plotblq) {
          linetypechoices <-  c("solid","dashed", "dotted", "dotdash", "longdash", "twodash","blank","solid","dashed", "dotted")
          selectizeInput(
            inputId = paste0('linetype',"blq"),
            label = paste0("Lines for %BLQ"),
            
            choices = linetypechoices,
            options = list(
              placeholder = 'Please select an option below',
              onInitialize = I('function() { this.setValue("solid"); }')
            )
          )        }       
        
      })
      
      
      
      vpcplotdata <- reactive({
        req(vpcData())
        checkboxfiltereddata<- vpcData()
        
        if(length(input$vpcformula) >= 1){
          PIsub<- vpcData()$PI
          BINSsub <-vpcData()$BINS
          DVCsub <- vpcData()$DVC
          PCTBLQsub <- vpcData()$PCTBLQ

          
          for (i in 1:(length(input$vpcformula) ) ){
            currentvar<- input$vpcformula[i]
            choices <-  eval(parse(text =  paste0("input$checkboxgrp",i))) 
            PIsub<-  PIsub [ PIsub[, currentvar] %in% choices  , ]
            BINSsub <-  BINSsub [ BINSsub[, currentvar] %in% choices  , ]
            DVCsub <- DVCsub [ DVCsub[, currentvar] %in% choices  , ]
            PCTBLQsub <- PCTBLQsub [ PCTBLQsub[, currentvar] %in% choices  , ]
            
          }
          checkboxfiltereddata$PI <- PIsub
          checkboxfiltereddata$BINS <- BINSsub
          checkboxfiltereddata$DVC <- DVCsub
          checkboxfiltereddata$PCTBLQ <- PCTBLQsub
          
        }
        
        if(length(input$vpcformula) == 0){
          checkboxfiltereddata
        }
        checkboxfiltereddata
      })

      
      obsplotdata <- reactive({
        checkboxfilteredobsdata<- obsdata
        if(length(input$vpcformula) >= 1){
          obsdatasub  <- obsdata
          for (i in 1:(length(input$vpcformula) ) ){
            currentvar<- input$vpcformula[i]
            choices <-  eval(parse(text =  paste0("input$checkboxgrp",i))) 
            obsdatasub <-   obsdatasub [ obsdatasub[, currentvar] %in% choices  , ]
          }
          checkboxfilteredobsdata <- obsdatasub
        }
        
        if(length(input$vpcformula) == 0){
          checkboxfilteredobsdata
        }
        checkboxfilteredobsdata
      })
      
      
      output$vpcplot <- renderPlot(height = 400, { # renderPlotly
         req(vpcplotdata())
        if(input$custompicolors){
          PIcolors<-  paste0("c(", paste0("input$col", 1:length(as.numeric(unlist(strsplit(input$PI, split=",")))), collapse = ", "), ")")
          PIcolors <- eval(parse(text = PIcolors))
        }
         if(input$custompilines){
           PIlinetypes<-  paste0("c(", paste0("input$linetypes", 1:length(as.numeric(unlist(strsplit(input$PI, split=",")))),
                                              collapse = ", "), ")")
           PIlinetypes <- eval(parse(text = PIlinetypes))
           }

     
        vpcplot <-  vpcplotfunction(vpcplotdata(), rectangles = input$rectangles,
                                    detachedrectangles = input$detachedrectangle,
                                    pointrange = input$pointrange,
                                    smoothlines = input$smoothlines,showrug = input$showrug, shownobs = input$shownobs,
                                    PItransparency=input$PItransparency,
                                    PIcolors =  PIcolors,
                                    PIlinetypes = PIlinetypes,
                                    obslinesizes=input$obslinesizes,
                                    bin_mid= input$bin_mid,
                                    facet_formula = facetformulaconstructed(),
                                    facetscales=input$facetscales,
                                    log_y= input$log_y,
                                    log_y_min=ifelse(input$log_y_min_set ,  input$log_y_min, -Inf),
                                    log_x= input$log_x,
                                    exp_y=input$exp_y,
                                    facetwrap = input$facetwrap,
                                    showLLOQ= ifelse(input$LLOQ == "None",FALSE,TRUE) ,
                                    legend.position = ifelse(input$legend.position,"none","right")
                                    )
        
        if(input$showobspoints&& !input$predcorrection&&!input$exp_y){
          vpcplot <- vpcplot +
            geom_point(data=obsplotdata(),aes_string(x= input$TIME, y ="DV"),alpha=0.2,size=3,shape=16)
        }
        
        if(input$showobspoints&& input$predcorrection&&!input$exp_y){
          vpcplot <- vpcplot +
            geom_point(data=vpcplotdata()$DVC,aes_string(x= input$TIME, y ="DVC"),alpha=0.2,size=3,shape=16)
        }
        
        if(input$showobspoints&& !input$predcorrection&&input$exp_y){
          vpcplot <- vpcplot +
            geom_point(data=obsplotdata(),aes_string(x= input$TIME, y ="exp(DV)"),alpha=0.2,size=3,shape=16)
        }
        
        if(input$showobspoints&& input$predcorrection&&input$exp_y){
          vpcplot <- vpcplot +
            geom_point(data=vpcplotdata()$DVC,aes_string(x= input$TIME, y ="exp(DVC)"),alpha=0.2,size=3,shape=16)
        }
        
        if(!input$plotblq) {
          outplot<- vpcplot
        }
        if(input$plotblq) {
          if(input$customblqcolors&&input$plotblq){
            blqcolor <- input$colblq
          }
          if(input$customblqlines&&input$plotblq){
            blqlinetype <- input$linetypeblq
          }
          blqplot <- pctblqplotfunction(vpcplotdata(), rectangles = input$rectangles,
                                     detachedrectangles = input$detachedrectangle,
                                     pointrange = input$pointrange,
                                     smoothlines = input$smoothlines,showrug = input$showrug, shownobs = input$shownobs,
                                     PItransparency=input$blqtransparency,
                                     PIcolors =  blqcolor,
                                     PIlinetypes = blqlinetype,
                                     obslinesizes=input$obsblqlinesize,
                                     bin_mid= input$bin_mid,
                                     facet_formula = facetformulaconstructed(),
                                     facetscales=input$facetscales,
                                     log_y= input$log_y,
                                     log_x= input$log_x,
                                     facetwrap = input$facetwrap,
                                     legend.position = ifelse(input$legend.position,"none","right")
          )
          outplot<- egg::ggarrange(vpcplot,blqplot)
          
        }
          
        outplot
        
        
      })
    },
    
    options = list(height = 500)
  )
}
