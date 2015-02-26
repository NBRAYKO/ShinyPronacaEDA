# ui.R
library(shiny)

shinyUI(fluidPage(
 
 titlePanel("Pronaca exploratory data analysis" ),
 sidebarLayout(
  
  
  #Sidebar panel------
  sidebarPanel(   
   img(src="http://i2.cpcache.com/image/66672167_125x125.png",height="auto", width="auto",inline=TRUE),
   h6("Filter Dataset"),
   selectInput("trial", 
              label = "Study (1a, 2a, 2b)",
              
              choices = c("1a Tet/New", "2a Enro/New", "2b Enro/Old"),
              multiple=TRUE,
              selected = "2a Enro/New"),  
   selectInput("tipo_muestra", 
              label = "Sample type (chicken or litter)",
              choices = c("Chicken", "Litter"),
              multiple=TRUE,
              selected = "Chicken"),
   selectInput("isolate_type", 
              label = "Isolate type (MH, oxa or enro media)",
              choices = c("SA", "EN", "OX"),
              selected = "SA"),

   h6("Choose Variables"),

  #select the first categorical variable (selection pressure, number of chickens, day etc. )
  selectInput("y_cat", 
              label = "Categorical outcome",
              choices = c("Res", "int1", "tetAB", "qnrb"),
              selected = "Res"),
  checkboxInput(
              inputId="class_ag",
              label = "Aggregate by class?",
              value = TRUE),
 
 
  h6("Additional Filtering"),
  
  selectInput("num_pollos", 
              label = "Birds in cage",
              multiple=TRUE,
              selectize=TRUE,
              choices = c("Single", "6", "34"),
              selected = c("Single", "6", "34")),
  selectInput("dia", 
              label = "Day of sampling",
              multiple=TRUE,
              selectize=TRUE,
              choices = c("01", "19", "21", "38"),
              selected = c("01", "19","21", "38")),
  selectInput("treat_full", 
              label = "Treatment",
              multiple=TRUE,
              selectize=TRUE,
              choices = c("Both", "Control",    "Food",   "Water" ),
              selected = c("Both", "Control",    "Food",   "Water")),
  width=2
  ),
 
#Main panel--------
 mainPanel(
  tabsetPanel(
   tabPanel("GENOTYPIC", 
            #GENO SUBTABS-------
            tabsetPanel(
             #subtab1
             tabPanel("Phen/gen dynamics",
                      withTags(div(class='row-fluid',
                                   div(class='span3', 
                                       selectInput("left_stack",
                                                   label="Left stack of stacked char",
                                                   choices = c("genotype",
                                                               "phen_desc_class","gen_str" ),
                                                   selected = "genotype")),
                                   div(class='span5',
                                       selectInput("right_stack",
                                                   label="Right stack",
                                                   choices = c("genotype",
                                                               "phen_desc_class","gen_str" ),
                                                   selected = "phen_desc_class")),
                                   selectInput("x1", label = "Category",choices = 
                                                c( "dia","num_pollos", "litter",
                                                   "treat_full", "treat_water", 
                                                   "treat_feed","trial","int1", 
                                                   "tetAB", "qnrb", "tipo_muestra" ),
                                               selected = "dia")
                      )),
                      plotOutput("PhenGenPlot",inline=TRUE)
             ),
             #subtab 2
             tabPanel("Ordination",
                      h6("Plot ordination fo subsetted data:"),
                      actionButton("ordiButton",
                                   label= img(src="https://img0.etsystatic.com/031/0/7321014/il_fullxfull.617240932_tvru.jpg",
                                              height = 20, width = 20)),
                      selectInput("x1", label = "Category",choices = 
                                   c( "dia","num_pollos", "litter",
                                      "treat_full", "treat_water", 
                                      "treat_feed","trial","int1", 
                                      "tetAB", "qnrb", "tipo_muestra" ),
                                  selected = "dia"),
                      selectInput("ordi_var",
                                  label="Similarity variable",
                                  choices = c("Genotypic (jaccard)", 
                                              "Genotypic (bray)",
                                              "Phenotypic" )),
                      sliderInput("range_x", "Range X:",step=0.5,
                                  min = -10, max = 10, value = c(-5,5)),
                      sliderInput("range_y", "Range Y:",step=0.5,
                                  min = -10, max = 10, value = c(-5,5)),
                      plotOutput("OrdiPlot")
             ),
             tabPanel("Adonis",                       
                      h6("Adonis model specification"),
                      helpText("To make model comprehensive, use filter tools to left to unsubset the data and pool results from trials; note that order of variables matters. "), 
                      selectInput("adonis_vars", 
                                  label = "Select variables:",
                                  multiple=TRUE,
                                  selectize=TRUE,
                                  choices = c("treat_full", "dia" ,"num_pollos", "litter"
                                              ,"tipo_muestra"  ,"int1*dia", "num_pollos*dia",
                                              "int1",  "tetAB", "tetAB*dia", "tetAB*num_pollos" 
                                              , "qnrb",
                                              "isolate_type", "trial"),
                                  selected = c("treat_full", "num_pollos","int1","dia")),
                      checkboxInput("adonis_strat",label="Stratify by trial?",
                                    value=FALSE),
                      h6("Run model"),
                      actionButton("adonisButton",
                                   label= img(src="https://img0.etsystatic.com/031/0/7321014/il_fullxfull.617240932_tvru.jpg",
                                             height = 20, width = 20)),
                      verbatimTextOutput("AdonisMod")
             )
            )
            
   ),
   tabPanel("PHENOTYPIC", 
            #PHENO SUBTABS-------
            tabsetPanel(
             tabPanel("Univariate plots and model",                     
                      withTags(div(class='row-fluid',
                                   selectInput("x1", label = "Category",choices = 
                                                c( "dia","num_pollos", "litter",
                                                   "treat_full", "treat_water", 
                                                   "treat_feed","trial","int1", 
                                                   "tetAB", "qnrb", "tipo_muestra" ),
                                               selected = "dia"),
                                   div(class='span5', 
                                       plotOutput("BarPlotSingle")),
                                   div(class='span5',conditionalPanel(condition = 
                                                                       "input.class_ag ==
                                                                      false",
                                       plotOutput("DensPlot"))
                                       ))), 
                      #select the drug for which you want detailed model results
                      #if results not aggregated by class, choose classes
                      #button to recaulculate model
                      h6("Reset model"),
                      actionButton("modelButton", 
                                   label=img(src="https://img0.etsystatic.com/031/0/7321014/il_fullxfull.617240932_tvru.jpg",
                                             height = 20, width = 20)
                      ),
                      checkboxInput("nest",label="Nested REs specificaiton (1|RE1/RE2)?",value=TRUE),
                      selectInput("Drug", 
                              label = "Drug-specific model results (scroll to bottom)",
                              choices = c("Ampicillin", "Amoxicillin/clavulanate"
                                          ,"Cefotaxime",  "Cephalothin", 
                                          "Chloramphenicol", "Ciprofloxacin"
                                          ,"Trimethoprim","Gentamicin", 
                                          "Streptomycin" , "Enrofloxacin",
                                          "Sulfisoxazole", "Tetracycline"),
                              selected="Enrofloxacin"
                      ),
                      
                      verbatimTextOutput("ModelSum")
                      ),
             tabPanel("Bivariate plots", 
                      
                      selectInput("x1", label = "Category",choices = 
                                   c( "dia","num_pollos", "litter",
                                      "treat_full", "treat_water", 
                                      "treat_feed","trial","int1", 
                                      "tetAB", "qnrb", "tipo_muestra" ),
                                  selected = "dia"),
                      selectInput("x2", label = "Grouping factor",
                                  choices = c("dia", "num_pollos","tipo_muestra", 
                                              "litter","int1", "tetAB", "qnrb"),   
                                  selected = "dia"),
                      selectInput("stacked",
                                    label="Bar chart type",
                                    choices = c("dodge","stack","fill" ),
                                    selected = "dodge"),
                       plotOutput("BarPlotDouble",inline=TRUE),
                       verbatimTextOutput("SummaryTab")
                      )
             ,tabPanel("MV model selection", 
                      selectInput("pheno_main", 
                                  label = "Select main effects:",
                                  multiple=TRUE,
                                  selectize=TRUE,
                                  choices = c("treat_full", "dia" ,"num_pollos", "litter"
                                              ,"tipo_muestra" ,"int1", "tetAB", "qnrb",
                                              "isolate_type", "trial" ),
                                  selected = c("dia","treat_full","num_pollos","int1", "tetAB", "qnrb")),
                      selectInput("fixed_main", 
                                  label = "Always include in models:",
                                  multiple=TRUE,
                                  selectize=TRUE,
                                  choices = c("treat_full", "dia" ,"num_pollos", "litter"
                                              ,"tipo_muestra" ,"int1", "tetAB", "qnrb",
                                              "isolate_type", "trial" ),
                                  selected = c("dia","treat_full")),
                      checkboxInput("time_int",label="Interact w time?",value=TRUE),
                      checkboxInput("tetAB_int",label="Interact w tetAB?",value=FALSE),
                      checkboxInput("int1_int",label="Interact w int1?",value=FALSE),
                      selectInput("REff", 
                                  label = "Select random effects:",
                                  choices = c("Random intercept by sample"= "(1|random.eff)",
                                              "Random intercept by bird"= "(1|random.effBig)", 
                                              "Nested intercept (1|bird/sample)"=
                                               "(1|random.effBig/random.eff)"
                                              )),
                      h6("Note: models will take a while to run (1--20 mins)"),
                      plotOutput("DredgeGraph",height=1200),
                      verbatimTextOutput("DredgeMod")
                      
             )
             ,selected=3)
            ),
   tabPanel("DATA TABLES", 
            tabsetPanel(
             tabPanel("Sample counts",
              h6("Produce sample counts:"),
              selectInput("tab_vars", 
                          label = "Select variables:",
                          multiple=TRUE,
                          selectize=TRUE,
                          choices = c("isolate_type","tipo_muestra", "trial", "dia",
                                      "treat_full" ,"num_pollos", "Res", "int1",
                                      "tetAB", "qnrb"),
                          selected = c("trial", "tipo_muestra","dia")),
              tableOutput("SampleCounts")),
             tabPanel("Sample counts",
              h6("Raw data"),
              dataTableOutput("RawData"))
            )
   ,selected=2),
  width=10
  )
 )
 ))