library(raster)
library(gstat)
library(sp)
library(dismo)
library(rgeos)
library(shiny)
library(shinyBS)
library(shinydashboard)
library(stringr)
library(shinyjs)
library(DT)
library(ggplot2)
library(deldir)


#setwd("/media/jberry/Extra Drive 1/Danforth/16S_Bacteria/Field/Shiny")

#*************************************************************************************************
# Getting data
#*************************************************************************************************
data <- read.table("data/all_2016_soil_data2.tsv",header = T,stringsAsFactors = F,sep = "\t")
data <- data[!stringr::str_detect(data$unique_id,"BK"),]
data <- data[-which(is.na(data$alpha_div)),]
data$coord_x <- as.numeric(substr(data$unique_id,3,4))*2 - (1*as.numeric(substr(data$unique_id,5,6))%%2)
data$coord_y <- ceiling(as.numeric(substr(data$unique_id,5,6))/2)
meta <- data[,c("treatment","sample_date")]
load("data/field_beta_diversity.RData")

ui <- dashboardPage(skin="purple",
                    dashboardHeader(
                      title = "Schachtman Fields",
                      titleWidth = 450
                    ),
                    dashboardSidebar(
                      disable = T
                    ),
                    dashboardBody(
                      tags$style(HTML("
                                      .tabbable > .nav > li > a[data-value='About'] {background-color: #676767;   color:white}
                                      .tabbable > .nav > li > a[data-value='Plot'] {background-color: #676767;   color:white}
                                      .tabbable > .nav > li > a[data-value='Tissue-Specific Heatmap'] {background-color: #676767;  color:white}
                                      .tabbable > .nav > li > a[data-value='Gene Finder'] {background-color: #676767; color:white}
                                      .tabbable > .nav > li > a[data-value='Contact Us'] {background-color: #676767; color:white}
                                      .tabbable > .nav > li[class=active]    > a {background-color: #444444; color:white}
                                      .multicol{
                                      -webkit-column-count: 4; /* Chrome, Safari, Opera */
                                      -moz-column-count: 4; /* Firefox */
                                      column-count: 4;
                                      }
                                      .twocol{
                                      -webkit-column-count: 2; /* Chrome, Safari, Opera */
                                      -moz-column-count: 2; /* Firefox */
                                      column-count: 2;
                                      }
                                      .warning { 
                                      color: red;
                                      }"
          )),
          uiOutput("login"),
          uiOutput("body")
                    )
)

server <- function(input, output) {
  output$login <- renderUI({
    showModal(modalDialog(easyClose = TRUE,footer = NULL,
                          passwordInput("password", "Password:"),
                          actionButton("go", "Go")
    ))
  })
  
  cc_auth <- reactiveValues(data="NO")
  
  observeEvent(input$go,{
    isolate({
      pass <- input$password
    })
    if(pass == "test"){
      cc_auth$data <- "YES"
      #      obs1$suspend()
      removeModal()
    }else{
      cc_auth$data <- "NO"
    }
  })
  
  output$body <- renderUI({
    if(cc_auth$data == "YES"){
      fluidRow(
        tabsetPanel(
          tabPanel("Field Soil",
                   fluidRow(
                     box(width=3,title = "Subsetting Data",solidHeader = T,status = 'success',collapsible = TRUE,
                         selectInput("which_field", width = 300,
                                     label = "Which Field: ",
                                     choices = c("Low N","Full N","WS","WW"),
                                     selected = "Low N"),
                         uiOutput("ui_which_date"),
                         hr(),
                         actionButton("use_config","Use these samples")
                     )
                   ),
                   tabsetPanel(
                     tabPanel("Measurement Associations",
                              uiOutput("corrplot"),
                              uiOutput("scatterplot")
                     ),
                     tabPanel("Spatial Analysis",
                              uiOutput("plots")
                     ),
                     tabPanel("Beta Diversity",
                              uiOutput("bd_ui_click"),
                              uiOutput("bd_ui_raster")
                     )
                   )
          )
        ) 
      )
    }
  })
  v <- reactiveValues(data=NULL)
  vca <- reactiveValues(data=NULL)
  r <- reactiveValues(data=NULL)
  vr <- reactiveValues(data=NULL)
  dsp <- reactiveValues(data=NULL)
  gs <- reactiveValues(data=NULL)
  nn <- reactiveValues(data=NULL)
  sub <- reactiveValues(data=NULL)
  sub1 <- reactiveValues(data=NULL)

  output$ui_which_date <- renderUI({
    sub <- data[data$treatment == input$which_field,]
    selectInput("which_date", width = 300,
                label = "Which Date: ",
                choices = c("all",unique(sub$sample_date)),
                selected = "all")
  })
  
  observeEvent(input$use_config,{
    if(input$which_date == "all"){
      temp <- data[data$treatment == input$which_field,]
    }else{
      temp <- data[data$treatment == input$which_field & data$sample_date == input$which_date,] 
    }
    sub$data <- temp
    sub1$data <- data.frame(sub$data,stringsAsFactors = F)
    sub1$data <- sub1$data[,c("Soil_pH","S_Salts","Nitrate_ppm","Nitrate_lbs","Potassium_ppm","Calcium_ppm","Magnesium_ppm","Sodium_ppm","Sulfate_ppm","CEC","H_Sat","K_Sat","Ca_Sat","Mg_Sat","Na_Sat","Mehlich.3_P","alpha_div")]
    
    sub1$data <- sub1$data[,as.logical(sapply(colnames(sub1$data), function(i) {
      if(length(na.omit(sub1$data[,i]))>0){
        if(var(sub1$data[,i],na.rm = T)>0){
          TRUE
        }else{
          FALSE
        }
      }else{
        FALSE
      }
      }))]

    dsp$data <- SpatialPoints(cbind(temp$coord_x,temp$coord_y))
    dsp$data <- SpatialPointsDataFrame(dsp$data,temp)

    v$data <- dismo::voronoi(dsp$data)
    vca$data <- intersect(v$data, dsp$data)
    r$data <- raster(dsp$data)
    
    gs$data <- gstat(formula=eval(parse(text = input$plot_by))~1, locations=vca$data, maxdist = 1.42)
  })
  
  output$bd_ui_click <- renderUI({
    if(!is.null(sub$data)){
      box(style = "overflow-y:scroll",width=5,title = "Click a sample",solidHeader = T,status = 'success',collapsible = TRUE,
          fluidRow(
            plotOutput("bd_click",click = "bd_clicked")
          )
      )
    }
  })

  selected <- reactiveValues(data = NULL)
  observe({
    selected$data <- nearPoints(sub$data, input$bd_clicked)
  })

  output$bd_click <- renderPlot({
    sub <- sub$data
    ggplot(sub, aes(coord_x,coord_y))+
      geom_point(size=3)+
      theme_light()
  })
  
  output$bd_ui_raster <- renderUI({
    if(!is.null(sub$data)){
      box(style = "overflow-y:scroll",width=5,title = "Raster",solidHeader = T,status = 'success',collapsible = TRUE,
          fluidRow(
            plotOutput("bd_raster_plot")
          )
      )
    }
  })
  
  output$bd_raster_plot <- renderPlot({
    sub <- sub$data
    sub <- dplyr::left_join(sub,m[m$c1 == as.character(selected$data$merge_id) & m$distance < 0.8,],by=c("merge_id"="c2"))
    sub <- sub[!is.na(sub$distance),]
    dsp <- SpatialPoints(cbind(sub$coord_x,sub$coord_y))
    dsp <- SpatialPointsDataFrame(dsp,sub)
    v <- dismo::voronoi(dsp)
    vca <- intersect(v, dsp)
    r <- raster(dsp)
    vr <- rasterize(vca, r, 'distance')
    gs <- gstat(formula=distance~1, locations=vca,maxdist=1.42)
    nn <- interpolate(r, gs)
    nnmsk <- mask(nn, vr)
    mmin <- min(sub$distance[sub$merge_id %in% m$c2]) 
    mmax <- max(sub$distance[sub$merge_id %in% m$c2])
    spacer <- (mmax-mmin)/8
    breakpoints <- round(seq(mmin-spacer,mmax+spacer,((mmax-spacer)-(mmin+spacer))/5),2)
    plot(nnmsk,xaxt = "n",yaxt = "n",main="Bray-Curtis Dissimilarity",col=rev(heat.colors(20)))
    axis(side = 1, at = 0:20, las = 2, hadj = 0.9)
    axis(side = 2, at = 0:15, las = 2)
  })
  
  output$meta_config <- renderDataTable({
    datatable(meta,rownames = F,selection = "none",filter = 'top',options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 5, autoWidth = TRUE))
  })
  
  output$corrplot <- renderUI({
    if(!is.null(dsp$data)){
      box(width=5,title = "Correlations",solidHeader = T,status = 'success',collapsible = TRUE,
          plotOutput("corrplot_plot")
      )
    }
  })
  
  output$scatterplot <- renderUI({
    if(!is.null(dsp$data)){
      box(style = "overflow-y:scroll",width=5,title = "Scatterplot",solidHeader = T,status = 'success',collapsible = TRUE,
        fluidRow(
          column(width=6,
                 selectInput("scat_x", width = 200,
                             label = "X-Axis Variable",
                             choices = colnames(sub1$data),
                             selected = "Nitrate_ppm")
          ),
          column(width=6,
                 selectInput("scat_y", width = 200,
                             label = "Y-Axis Variable",
                             choices = colnames(sub1$data),
                             selected = "Potassium_ppm")
          )  
        ),
        fluidRow(
          plotOutput("scatter_out")
        )
      )
    }
  })
  
  output$scatter_out <- renderPlot({
    sub <- sub1$data
    x <- input$scat_x
    y <- input$scat_y
    ggplot(sub,aes(eval(parse(text=x)),eval(parse(text=y))))+
      geom_point()+
      theme_light()+
      theme(plot.margin=unit(c(1,1,1,1), "cm"))+
      theme(axis.text = element_text(size = 14),
            axis.title= element_text(size = 18))+
      xlab(x)+
      ylab(y)
  })
  
  output$corrplot_plot <- renderPlot({
    mat <- cor(sub1$data)
    corrplot::corrplot(mat,type = "upper",order = "hclust",tl.col = "black")
  })
  
  output$plots <- renderUI({
    if(!is.null(dsp$data)){
      fluidRow(
        box(width=2,solidHeader = T,status = 'success',
          selectInput("plot_by", width = 300,
                      label = "Plot by: ",
                      choices = colnames(sub1$data),
                      selected = "Nitrate")
        ),
        box(width=5,title = "Raster Field",solidHeader = T,status = 'success',collapsible = TRUE,
            plotOutput("Raster")   
        ),
        box(width=5,title = "Variogram",solidHeader = T,status = 'success',collapsible = TRUE,
            plotOutput("Variogram")
        )
      ) 
    }
  })
  
  output$Raster <- renderPlot({
    vr <- rasterize(vca$data, r$data, input$plot_by)
    nn <- interpolate(r$data, gs$data)
    nnmsk <- mask(nn, vr)
    mmin <- min(nnmsk@data@values,na.rm = T) 
    mmax <- max(nnmsk@data@values,na.rm = T)
    spacer <- (mmax-mmin)/8
    breakpoints <- round(seq(mmin-spacer,mmax+spacer,((mmax-spacer)-(mmin+spacer))/5),2)
    plot(nnmsk,breaks=breakpoints,col=rev(terrain.colors(length(breakpoints))),xaxt = "n",yaxt = "n")
    axis(side = 1, at = 0:20, las = 2, hadj = 0.9)
    axis(side = 2, at = 0:15, las = 2)
  })
  
  output$Variogram <- renderPlot({
    v <- variogram(gs$data)
    plot(v)
  })
  
}


shinyApp(ui, server)