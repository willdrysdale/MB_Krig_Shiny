ui <- fluidPage(

   sidebarLayout(
      sidebarPanel(width = 3,
                   fileInput("file1","Choose File",
                             accept = c("text/csv","text/comma-separated-values,text/plain",".csv")
                   ),
                   conditionalPanel("input.data_ok == false",
                                    uiOutput("krig_param_ui"),
                                    conditionalPanel(condition = "output.mix_to_mgm3_option == true",
                                      checkboxInput("mix_to_mgm3","Change mixing ratio to mass volume"),
                                      conditionalPanel(condition = "input.mix_to_mgm3 == true",
                                                       selectInput("mix_to_mgm3_unit",choice = c("ppm","ppb","ppt"),
                                                                   selected = "ppt",label = "Mixing ratio input unit"),
                                                       numericInput("mix_to_mgm3_mass",min = 0,step = 1,
                                                                    label = "Mixing ratio species mass mg/m3",value = 46005)
                                                       )
                                    ),
                                    uiOutput("x_dim_ui"),
                                    uiOutput("y_dim_ui"),
                                    checkboxInput("x_lat_to_m","Convert Lat to meters",value = F),
                                    checkboxInput("resc_axis","Rescale Axis Range",value = F),
                                    tags$h5("Note: Fine grids will take a long time to krig"),
                                    uiOutput("x_grid_dim_ui"),
                                    uiOutput("y_grid_dim_ui"),
                                    uiOutput("anis_x_slider"),
                                    uiOutput("anis_y_slider"),
                                    uiOutput("check_data_ui")
                         ),
                   uiOutput("check_data_ok_ui"),
                   conditionalPanel(condition = "input.data_ok == true",
                                    actionButton(inputId = "vario_click",label = "Fit Variogram"),
                                    checkboxInput(inputId = "vario_ok",label = "Variogram Correct?",value = F)
                   ),
                   conditionalPanel(condition = "input.data_ok == true & input.vario_ok == true",
                                    actionButton(inputId = "krige",label = "Perform Kriging")
                   ),
                   actionButton(inputId = "update_plots",label = "Update Plots")
                   
      ),
      

      mainPanel(
        tabsetPanel(
          tabPanel("Inspect Data",
                   dataTableOutput("input_data")
                   ),
          tabPanel("View Data",
                   plotOutput("timeseries"),
                   uiOutput("ts_colour_ui"),
                   plotOutput("timeseries_processed")
                   ),
          tabPanel("Variogram",
                   plotOutput("vario_plot"),
                   verbatimTextOutput("vario_info")
                   ),
          tabPanel("Anisotropy",
                   plotOutput("anis_plot")
          ),
          tabPanel("Kriging",
                   plotOutput("krige_plot"),
                   plotOutput("krige_plot_interp"),
                   tags$h5("Grid Box Size"),
                   textOutput("grid_box_size"),
                   downloadButton(outputId = "export_data", "Export Data")
                   ),
          tabPanel("Kriging Variance",
                   plotOutput("krige_plot_var"),
                   plotOutput("krige_plot_var_interp")
          )
        )
      )
   )
)


server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2) 
  
  #Load Packages
  library(ggplot2)
  library(magrittr)
  library(lubridate)
  library(openair)
  library(viridisLite)
  library(dplyr)
  library(shiny)
  library(gstat)
  library(sp)
  library(wsdmiscr)
  
  #Format Reactive Data
  raw_input = reactive({
    req(input$file1)
    file_input = read.csv(input$file1$datapath)
    choices = names(file_input)
    
    raw_input = list(
      file_input = file_input,
      choices = choices
    )
    
    raw_input
  })
  #Allow option to convert mixing ratio to mgm3 automatically if the required columns exist in the input dataset
  output$mix_to_mgm3_option = reactive({
    all(c("ps_rvsm","tat_di_r") %in% raw_input()$choices)
    })
  outputOptions(output, "mix_to_mgm3_option", suspendWhenHidden = FALSE)
    
  
  data_list = reactiveValues()
  observeEvent(input$create_data,{
    withProgress(message = "Creating Data",value = 0,{
      incProgress(0.33)
      data_file_all = isolate({raw_input()$file_input})
      data_file = data_file_all[,c(isolate({input$x_dim}),isolate({input$y_dim}),isolate({input$krig_param}))]
      names(data_file) = c("x_dim","y_dim","krig_param")
      orig_dim = data.frame(x_orig = data_file$x_dim,y_orig = data_file$y_dim)
      #convert units of krig_param
      if(input$mix_to_mgm3){
        data_file$krig_param = wsdmiscr::faam_mixing_ratio_to_mgm3(df = data_file_all,
                                                                   unit = isolate({input$mix_to_mgm3_unit}),
                                                                   pollutant = isolate({input$krig_param}),
                                                                   pollutant_mass = isolate({input$mix_to_mgm3_mass}))
      }
      #Convert Latitude to meters
      if(isolate({input$x_lat_to_m})){
        mean_lat = mean(data_file$x_dim,na.rm = T)
        data_file$x_dim = data_file$x_dim - min(data_file$x_dim)
        data_file$x_dim = data_file$x_dim * wsdmiscr::length_of_lat(mean_lat)
        }
      #Equate Axis Lengths
      if(isolate({input$resc_axis})){
        x_range = max(data_file$x_dim,na.rm = T)-min(data_file$x_dim,na.rm = T)
        y_range = max(data_file$y_dim,na.rm = T)-min(data_file$y_dim,na.rm = T)
        
        if(x_range > y_range){
          ratio = y_range/x_range
          data_file$x_dim = data_file$x_dim * ratio
          data_file$x_dim = data_file$x_dim - min(data_file$x_dim,na.rm = T)
        }
        if(y_range > x_range){
          ratio = x_range/y_range
          data_file$y_dim = data_file$y_dim * ratio
          data_file$y_dim = data_file$y_dim - min(data_file$y_dim,na.rm = T)
        }
      }
      data_file$x_dim = data_file$x_dim * isolate({input$anis_x_scale})
      data_file$y_dim = data_file$y_dim * isolate({input$anis_y_scale})
      
      
      #Create SPDF of input data
      coords = data_file[,c("x_dim","y_dim")]
      data_spdf = sp::SpatialPointsDataFrame(coords = coords[complete.cases(data_file),],data = data_file[complete.cases(data_file),])
      incProgress(0.33)
      #Create grid to krig over
      data_bbox = data_spdf@bbox
      x_dim_vec = seq(floor(data_bbox[1]),ceiling(data_bbox[3]),isolate({input$x_grid_dim}))
      y_dim_vec = seq(floor(data_bbox[2]),ceiling(data_bbox[4]),isolate({input$y_grid_dim}))
      
      for(i in 1:length(x_dim_vec)){
        chunk = data.frame(x_dim = rep(x_dim_vec[i],length(y_dim_vec)),
                           y_dim = y_dim_vec)
        if(i == 1)
          grid_coords = chunk
        else
          grid_coords = rbind(grid_coords,chunk)
      }
      incProgress(0.33)
      grid = grid_coords
      grid$data = NA
      
      grid_spdf = sp::SpatialPointsDataFrame(grid_coords,grid)
      incProgress(0.34)
      data_list$data = data_file
      data_list$data_all = data_file_all
      data_list$data_spdf = data_spdf
      data_list$coords = coords
      data_list$grid = grid_spdf
      #Calculate grid box size
      data_list$grid_box_size = (x_range/length(x_dim_vec))*(y_range/length(y_dim_vec))
    })
  })
  
  vario_list = reactiveValues()
  
  observeEvent(input$vario_click,{
    withProgress(message = "Fitting Variogram",value = 0,{
      vario_model = paste("krig_param","1",sep = "~") %>% as.formula()
      vario_1d = gstat::variogram(vario_model,data_list$data_spdf)
      incProgress(0.5)
      vario_radial = gstat::variogram(vario_model,data_list$data_spdf,alpha = seq(0,360,10))
      incProgress(0.75)
      v.fit = gstat::fit.variogram(vario_radial, gstat::vgm(c("Sph","Gau","Exp","Mat","Ste","Cir","Lin","Bes","Pen","Wav")))

      vario_list$vario_1d = vario_1d
      vario_list$vario_radial = vario_radial
      vario_list$vario_fit = v.fit
      vario_list$vario_model = vario_model

      incProgress(0.25)

    })
  })
  

  
  #Perform Kriging
  kriged_data = reactiveValues()
  observeEvent(input$krige,{
    withProgress(message = "Kriging",value = 0.3,{
      x_kriged = gstat::krige(isolate({vario_list$vario_model}),
                              isolate({data_list$data_spdf}),
                              isolate({data_list$grid}),
                              model = isolate({vario_list$vario_fit}),
                              debug.level = -1)
      incProgress(0.7)
    
    
      kriged_data$df =  as.data.frame(x_kriged)
      
    })
  })
  #Create Export/Download Items
  
  output$export_data = downloadHandler(
    filename = function() {
      paste(isolate({input$krig_param}),
            isolate({input$x_grid_dim}),
            isolate({input$y_grid_dim}),"data.csv",sep ="_")
      },
    content = function(file) {
      write.csv(isolate({kriged_data$df}),file,row.names = F)}
    )

  
  #Create Reactive UI
  output$krig_param_ui = renderUI({
    selectInput(inputId = "krig_param", label = "Kriging Param:",
                choices = raw_input()$choices, selected = raw_input()$choices[3])
  })
  
  output$x_dim_ui = renderUI({
    selectInput(inputId = "x_dim", label = "x Dimention:",
                choices = raw_input()$choices, selected = raw_input()$choices[2])
  })
  
  output$y_dim_ui = renderUI({
    selectInput(inputId = "y_dim", label = "y Dimention:",
                choices = raw_input()$choices, selected = raw_input()$choices[1])
  })
  
  output$x_grid_dim_ui = renderUI({
    req(input$file1)
    numericInput(inputId = "x_grid_dim", label = "x Gridding Dimention :", value = 50, min = 1,step = 1)
  })
  
  output$y_grid_dim_ui = renderUI({
    req(input$file1)
    numericInput(inputId = "y_grid_dim", label = "y Gridding Dimention :", value = 50, min = 1,step = 1)
  })
  
  output$anis_x_slider = renderUI({
    req(input$file1)
    sliderInput(inputId = "anis_x_scale",label = "Anisotropy x scaler",min = 1,max = 10,value = 1)
  })
  output$anis_y_slider = renderUI({
    req(input$file1)
    sliderInput(inputId = "anis_y_scale",label = "Anisotropy y scaler",min = 1,max = 10,value = 1)
  })
  
  output$check_data_ui = renderUI({
    req(input$file1)
    actionButton(inputId = "create_data",label = "Create Data")
  })
  output$check_data_ok_ui = renderUI({
    req(input$file1)
    checkboxInput(inputId = "data_ok",label = "Data Correct?",value = F)
  })
  
  output$ts_colour_ui = renderUI({
    selectInput(inputId = "ts_colour", label = "Colour by: ",
                choices = c("none",raw_input()$choices),
                selected = isolate({input$krig_param}))
  })
  
  output$grid_box_size = renderText({data_list$grid_box_size
    })
  
  #Create Visual output
  
  fill_col = "grey97"
  line_col = "black"
  gen_theme = theme(
    plot.background = element_rect(fill = fill_col, colour = "grey92"),
    panel.background = element_rect(fill = fill_col, colour = fill_col),
    panel.grid.major = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(colour = line_col,size = 25,hjust = 0.5),
    axis.text = element_text(colour = line_col,size = 20,face = "bold"),
    axis.title = element_text(colour = line_col,size = 25),
    axis.ticks = element_line(colour = line_col),
    axis.line = element_line(colour = line_col),
    panel.border = element_rect(colour = line_col,fill = NA)
  )
  
  output$input_data = renderDataTable({
    data_list$data
  },options = list(pageLength = 10))
  
  output$timeseries = renderPlot({
    input$update_plots
    if(input$ts_colour == "none"){
      ggplot(
        isolate({data_list$data_all}))+
        geom_point(aes_string(isolate({input$x_dim}),isolate({input$y_dim})))+
        ylab(isolate({input$y_dim}))+
        xlab(isolate({input$x_dim}))+
        gen_theme
    }else{
      ggplot(
        isolate({data_list$data_all}))+
        geom_point(aes_string(isolate({input$x_dim}),
                              isolate({input$y_dim}),
                              col = isolate({input$ts_colour})
                              )
                   )+
        ylab(isolate({input$y_dim}))+
        xlab(isolate({input$x_dim}))+
        gen_theme+
        scale_colour_gradientn(colours = viridis(200))
    }
  },height = 400, width = 600)
  
  output$timeseries_processed = renderPlot({
    input$update_plots
    ggplot(
      isolate({data_list$data}))+
      geom_point(aes_string("x_dim","y_dim",colour = "krig_param"))+
      ylab(isolate({input$y_dim}))+
      xlab(isolate({input$x_dim}))+
      scale_colour_gradientn(colours = viridis(200))+
      gen_theme
  },height = 400, width = 600)
  
  output$vario_plot = renderPlot({
    input$update_plots
    plot(vario_list$vario_1d,vario_list$vario_fit,lw = 2,pch = 16,cex = 1.1,col = "#440154FF")
  },height = 400, width = 600)
  
  output$vario_info = renderPrint({
    vario_list$vario_fit
  }, width = 600)
  
  output$anis_plot = renderPlot({
    input$update_plots
    polarPlot(isolate({vario_list$vario_radial}),
              pol = "gamma",x = "dist",wd = "dir.hor",col = inferno(200))
  },height = 400, width = 600)
  
  
  output$krige_plot = renderPlot({
    input$update_plots
      ggplot(
        isolate({kriged_data$df}),
        aes(x=x_dim, y=y_dim))+
        ylab(isolate({input$y_dim}))+
        xlab(isolate({input$x_dim}))+
        ggtitle(isolate({input$krig_param}))+
        geom_raster(aes(fill=var1.pred))+
        scale_fill_gradientn(colours = viridis(200))+
        gen_theme
  },height = 400, width = 600)
  
  output$krige_plot_interp = renderPlot({
    input$update_plots
    ggplot(
      isolate({kriged_data$df}),
      aes(x=x_dim, y=y_dim))+
      ylab(isolate({input$y_dim}))+
      xlab(isolate({input$x_dim}))+
      geom_raster(aes(fill=var1.pred),interpolate = T)+
      stat_contour(aes(z = var1.pred),col = "white")+
      ggtitle(isolate({input$krig_param}))+
      scale_fill_gradientn(colours = viridis(200))+
      gen_theme
  },height = 400, width = 600)
  
  output$krige_plot_var = renderPlot({
    input$update_plots
    ggplot(
      isolate({kriged_data$df}),
      aes(x=x_dim, y=y_dim))+
      ylab(isolate({input$y_dim}))+
      xlab(isolate({input$x_dim}))+
      geom_raster(aes(fill=var1.var))+
      ggtitle(paste0(isolate({input$krig_param})," variance"))+
      scale_fill_gradientn(colours = magma(200))+
      gen_theme
  },height = 400, width = 600)
  
  output$krige_plot_var_interp = renderPlot({
    input$update_plots
    ggplot(
      isolate({kriged_data$df}),
      aes(x=x_dim, y=y_dim))+
      ylab(isolate({input$y_dim}))+
      xlab(isolate({input$x_dim}))+
      ggtitle(paste0(isolate({input$krig_param})," variance"))+
      geom_raster(aes(fill=var1.var),interpolate = T)+
      stat_contour(aes(z = var1.var),col = "white")+
      scale_fill_gradientn(colours = magma(200))+
      gen_theme
  },height = 400, width = 600)
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

