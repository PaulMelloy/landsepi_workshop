library(shiny)
library(landsepi)
library(DT)
source("../R/rolling_vec.R")
simul_params <- createSimulParams(outputDir = "../sims/")
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Landsepi workshop UQ"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          sliderInput("Nyears",
                      "Years in simulation",
                      min = 1,
                      max = 50,
                      value = 10,
                      step = 1),
          numericInput("season",
                       "Cropping season",
                       value = 180,
                       min = 60,
                       max = 300,
                       step = 1),
          h2("Pathogen Parameters"),
          numericInput("latent_per",
                       "Latent period",
                       value = 10,
                       min = 3,
                       max = 60,
                       step = 1),
          numericInput("infectious_per",
                       "Infectious period",
                       value = 24,
                       min = 7,
                       max = 200,
                       step = 1),

          actionButton(inputId = "set_params",
                       label = "Set Parameters"),
          h2("Landscape and crop rotation"),
          h3("Rotation proportions"),
          numericInput("p_forest",
                       "Proportion Forest",
                       min = 0,
                       max = 0.8,
                       step = 0.01,
                       0.2),
          numericInput("p_cropping",
                       "Proportion cropping",
                       min = 0,
                       max = 0.8,
                       step = 0.01,
                       value = 0.5),
          numericInput("p_pasture",
                       "Proportion permanent pasture",
                       min = 0,
                       max = 0.8,
                       step = 0.01,
                       value = 0.3),
          h3("Cropping sequence"),
          p("cropID number, seperated by a comma ',' no spaces "),
          textInput("crop_seq",
                    "Cropping sequence",
                    "5,3,3,5,6"),
          actionButton(inputId = "set_rotation",
                       label = "Set Final Parameters"),
          actionButton(inputId = "run_sim",
                       label = "Run simulation")
        ),

        mainPanel(
          h2("Main Panel") ,
          #textOutput(outputId = "dir"),
           dataTableOutput("cultivar_types")
        )
)
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$dir <- renderText(getwd())
    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    #
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white',
    #          xlab = 'Waiting time to next eruption (in mins)',
    #          main = 'Histogram of waiting times')
    # })

  # Prepare parameters

  sim_par <-
    reactive({
cat("here\n")
    simul_params <- setTime(simul_params,
                            Nyears = input$Nyears,
                            nTSpY = input$season)
    # pathogen
    basic_patho_param <- loadPathogen(disease = "rust")
    basic_patho_param$latent_period_mean <- input$latent_per
    basic_patho_param$infectious_period_mean <- input$infectious_per

    simul_params <- setPathogen(simul_params, patho_params = basic_patho_param)
    simul_params <- setInoculum(simul_params, val = 5e-4)

    # Landscape
    landscape <- loadLandscape(id = 2)
    disp_patho_clonal <- loadDispersalPathogen(id = 2)[[1]]
    #length(landscape)^2 == length(disp_patho_clonal) # Check matrix matches
    simul_params <- setLandscape(simul_params, land = landscape)
    simul_params <- setDispersalPathogen(simul_params, disp_patho_clonal)

    # Crops
    cultivar1 <- loadCultivar(name = "VS", type = "growingHost")
    cultivar2 <- loadCultivar(name = "R", type = "growingHost")
    cultivar3 <- loadCultivar(name = "MR_APR", type = "growingHost")
    cultivar4 <- loadCultivar(name = "MSMR_APR", type = "growingHost")
    cultivar5 <- loadCultivar(name = "MS1", type = "growingHost")
    cultivar6 <- loadCultivar(name = "MSMR1", type = "growingHost")
    cultivar7 <- loadCultivar(name = "MR", type = "growingHost")
    cultivar8 <- loadCultivar(name = "MSMR2", type = "growingHost")
    cultivar9 <- loadCultivar(name = "MS2", type = "growingHost")
    cultivar10 <- loadCultivar(name = "S", type = "growingHost")
    Fallow1 <- loadCultivar(name = "UnManaged_Fallow", type = "growingHost")
    Fallow2 <- loadCultivar(name = "Managed_Fallow", type = "nonCrop")
    Forest <- loadCultivar(name = "Other", type = "nonCrop")
    cultivars <- data.frame(rbind(cultivar1, cultivar2, cultivar3, cultivar4, cultivar5,
                                  cultivar6, cultivar7, cultivar8, cultivar9, cultivar10,
                                  Fallow1,Fallow2, Forest)
                            , stringsAsFactors = FALSE)

    cultivars[cultivars$cultivarName == "UnManaged_Fallow", "growth_rate"] <- 0.01
    cultivars[cultivars$cultivarName == "UnManaged_Fallow", "max_density"] <- 0.4
    cultivars[cultivars$cultivarName == "UnManaged_Fallow", "initial_density"] <- 0.01
    cultivars[cultivars$cultivarName == "UnManaged_Fallow", "yield_H"] <- 0
    cultivars[cultivars$cultivarName == "UnManaged_Fallow", "planting_cost"] <- 0
    cultivars[cultivars$cultivarName == "UnManaged_Fallow", "market_value"] <- 0


    # Resistance genes
    gene1.1 <- loadGene(name = "MG 1", type = "majorGene")
    gene1.1$mutation_prob <- 1e-2
    gene1.2 <- loadGene(name = "MG 2", type = "majorGene")
    gene1.2$mutation_prob <- 1e-3
    gene2.1 <- loadGene(name = "APR1", type = "APR")
    gene2.1$time_to_activ_mean <- 60
    gene2.1$time_to_activ_var <- 10
    gene2.1$efficiency <- 0.9
    gene2.2 <- loadGene(name = "APR2", type = "APR")
    gene2.1$time_to_activ_mean <- 90
    gene2.1$time_to_activ_var <- 10
    gene2.1$efficiency <- 0.9
    gene3.1 <- loadGene(name = "QTL 1", type = "QTL")
    gene3.1$time_to_activ_mean <- 10
    gene3.1$time_to_activ_var <- 10
    gene3.1$mutation_prob <- 1e-12
    gene3.1$efficiency <- 0.4
    gene3.2 <- loadGene(name = "QTL 2", type = "QTL")
    gene3.2$time_to_activ_mean <- 1
    gene3.2$time_to_activ_var <- 10
    gene3.2$mutation_prob <- 1e-12
    gene3.2$efficiency <- 0.2
    gene3.3 <- loadGene(name = "QTL 3", type = "QTL")
    gene3.3$time_to_activ_mean <- 10
    gene3.3$time_to_activ_var <- 10
    gene3.3$mutation_prob <- 1e-12
    gene3.3$efficiency <- 0.1
    gene3.4 <- loadGene(name = "QTL 4", type = "QTL")
    gene3.4$time_to_activ_mean <- 20
    gene3.4$time_to_activ_var <- 10
    gene3.4$mutation_prob <- 1e-12
    gene3.4$efficiency <- 0.65

    gene4 <- loadGene(name = "nonhost resistance", type = "immunity")
    genes <- data.frame(rbind(gene1.1,gene1.2,
                              gene2.1, gene2.2,
                              gene3.1,gene3.2, gene3.3,gene3.4,
                              gene4), stringsAsFactors = FALSE)
    simul_params <- setGenes(simul_params, dfGenes = genes)
    simul_params <- setCultivars(simul_params, dfCultivars = cultivars)

    # Allocate genes to cultivars
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "R"
                                          , listGenesNames = c("MG 1"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "MR_APR"
                                          , listGenesNames = c("APR1"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "MSMR_APR"
                                          , listGenesNames = c("APR2"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "MS1"
                                          , listGenesNames = c("QTL 2"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "MSMR1"
                                          , listGenesNames = c("QTL 2","QTL 1"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "MR"
                                          , listGenesNames = c("QTL 4","QTL 1"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "MSMR2"
                                          , listGenesNames = c("QTL 4"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "MS2"
                                          , listGenesNames = c("QTL 2","QTL 3"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "S"
                                          , listGenesNames = c("QTL 3"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "UnManaged_Fallow"
                                          , listGenesNames = c("QTL 3"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "Managed_Fallow"
                                          , listGenesNames = c("nonhost resistance"))
    simul_params <- allocateCultivarGenes(simul_params
                                          , cultivarName = "Other"
                                          , listGenesNames = c("nonhost resistance"))
    # Croptypes
    croptypes <- loadCroptypes(simul_params, names = c("Other",
                                                       "pasture",
                                                       "unmanaged_fallow",
                                                       "managed_fallow",
                                                       "monoculture_R",
                                                       "monoculture_MR_APR",
                                                       "monoculture_MSMR_APR",
                                                       "monoculture_MS1",
                                                       "monoculture_MSMR1",
                                                       "monoculture_MR",
                                                       "monoculture_MSMR2",
                                                       "monoculture_MS2",
                                                       "monoculture_S"))
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_R",
                                           cultivarsInCroptype = "R")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_MR_APR",
                                           cultivarsInCroptype = "MR_APR")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_MSMR_APR",
                                           cultivarsInCroptype = "MSMR_APR")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_MS1",
                                           cultivarsInCroptype = "MS1")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_MSMR1",
                                           cultivarsInCroptype = "MSMR1")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_MR",
                                           cultivarsInCroptype = "MR")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_MSMR2",
                                           cultivarsInCroptype = "MSMR2")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_MS2",
                                           cultivarsInCroptype = "MS2")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "monoculture_S",
                                           cultivarsInCroptype = "S")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "unmanaged_fallow",
                                           cultivarsInCroptype = "UnManaged_Fallow")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "managed_fallow",
                                           cultivarsInCroptype = "Managed_Fallow")
    croptypes <- allocateCroptypeCultivars(croptypes
                                           , croptypeName = "Other"
                                           , cultivarsInCroptype = "Other")
    croptypes <- allocateCroptypeCultivars(croptypes,
                                           croptypeName = "pasture",
                                           cultivarsInCroptype = c("UnManaged_Fallow","Managed_Fallow"),
                                           prop = c(0.3,0.7))



    simul_params <- setCroptypes(simul_params, dfCroptypes = croptypes)
    simul_params


  })


  ct_df <-
    eventReactive(input$set_params,{
      DT::datatable(sim_par()@Croptypes[,1:2])
      })
  output$cultivar_types <- renderDataTable(ct_df())


  # observeEvent(input$cultivar_types_rows_selected,{
  # #   output$rotation_text <- renderText("cultivar_types_rows_selected")
  # output$rotation_text <- renderPrint(cat(sim_par()@Croptypes$croptypeName[input$cultivar_types_rows_selected], sep = ", \n"))
  # })

  # convert crop sequence to a numeric vector
  crop_seq1 <- reactive(as.numeric(unlist(strsplit(input$crop_seq, ","))))

  rotation_ind <- reactive({
    sample(seq_along(crop_seq1()),
           replace = FALSE,
           size = length(crop_seq1())
  )})

  rotation_seq <- reactive({
    lapply(seq_along(crop_seq1()), function(x) {
      # Forest_prop, pasture_prop, cropping_prop
      c1 <- rolling_vec(crop_seq1(),x)[rotation_ind()]
      c(0, 1, c1)
    })
  })

  # proportion in surface
  prop <- reactive({
    lapply(crop_seq1(),function(x2){
      props <- c(input$p_forest,
        input$p_pasture,
        rep(input$p_cropping/length(crop_seq1()),
            length(crop_seq1())))

      if(sum(props)!= 1)stop("proportions don't sum to 1")
      return(props)
    })
  })

  sim_par2 <-
    eventReactive(input$set_rotation ,{
      cat("start_check")
    aggreg <- vector(mode = "list",
                     length = length(crop_seq1))
    aggreg[1:length(crop_seq1())] <- 0.2
    sim_params <-
        allocateLandscapeCroptypes(
          sim_par(),
          rotation_period = 1,
          rotation_sequence = rotation_seq(),
          prop = prop(),
          aggreg = aggreg,
          graphic = TRUE
        )
    cat("end_check")
    sim_params
  })
  observe(sim_par2())

  # Run simulation
  observeEvent(input$run_sim ,{
    # specify outputs
    outputlist <- loadOutputs(epid_outputs = c("audpc_rel","gla_rel", "eco_yield"), evol_outputs = "all")

    #update sim_par
    simulation_pars <- setOutputs(sim_par2(), outputlist)

    checkSimulParams(simulation_pars)
    simulation_pars <- saveDeploymentStrategy(simulation_pars)

    sim <- runSimul(simulation_pars, graphic = TRUE, videoMP4 = FALSE)

    cat("end")
  })





}

# Run the application
shinyApp(ui = ui, server = server)



