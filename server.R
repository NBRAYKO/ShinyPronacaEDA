library(shiny)
library(plyr)
source("Pronaca functions.R")
load("PronacaClean.Rdat")
y<-"Zone"
by_var<-"Drug"
random.eff<-"sample_id"
random.effBig<-"pollo_id"
y_cat<-"Res"

shinyServer(function(input, output){
 
 #setting the plotting variables 	
  y<-"Zone"
  by_var<-"Drug"
  random.eff<-"sample_id"
  random.effBig<-"pollo_id"
  y_cat<-reactive({input$y_cat})
  x1<-reactive({input$x1})
  x2<-reactive({input$x2})
  stacked<-reactive({input$stacked})
  Drug<- reactive({input$Drug})
  Class<- reactive({input$Class})
  left_stack<- reactive({input$left_stack})
  right_stack<- reactive({input$right_stack})
  ordi_var<- reactive({input$ordi_var})
  adonis_vars<-reactive({input$adonis_vars})
  class_ag<-reactive({input$class_ag})
  nest<-reactive({input$nest})
  adonis_strat<-reactive({input$adonis_strat})
  range_x<-reactive({input$range_x})
  range_y<-reactive({input$range_y})
  pheno_main<-reactive({input$pheno_main})
  fixed_main<-reactive({input$pheno_main})
  #time_int<-reactive({input$time_int})
  int1_int<-reactive({input$int1_int})
  tetAB_int<-reactive({input$tetAB_int})
  REff<- reactive ({input$REff})
  tab_vars<- reactive ({input$tab_vars})
  
  
  Data <- reactive({  
       #if the outcome variable is not %R, but one of the genes:
       if (input$y_cat != "Res") {
        dfLo=subset(dfLo,dfLo$Drug=="Ampicillin")
        dfLo$Drug="All drugs"
       }
       #if we want to aggregate 
       if (input$class_ag ==TRUE) {
        #collapse at sample level for following varlist
        CollapseVars=c("drug_class","trial","dia", "sample_id", "tipo_muestra", "isolate_type", "num_pollos", "treat_full", "litter","int1", "tetAB", "qnrb", "sample_id", "pollo_id","isolate_code", "id" )
        CollapseForm=as.formula(paste("cbind(Res,Zone)","~",paste(CollapseVars, collapse="+")))
        dfLo=aggregate(CollapseForm,data=dfLo,FUN=max)
        dfLo$Drug=dfLo$drug_class
        if (input$y_cat != "Res") dfLo$Drug="All drugs" #set label to "all drugs" if collapsing
       }
       Data=dfLo[which(dfLo$trial%in%input$trial &
       dfLo$tipo_muestra %in% input$tipo_muestra &
       dfLo$isolate_type==input$isolate_type &
       dfLo$num_pollos%in%input$num_pollos &
       dfLo$dia%in%input$dia &
       dfLo$treat_full%in%input$treat_full),]
       #add somthing here to collapse data by class if a box is ticked; if input$class_ag==TRUE
       
  })
  #wide data 
  DataWide <-reactive({  
   DataWide=dfWi[which(dfWi$trial%in%input$trial &
                    dfWi$tipo_muestra %in% input$tipo_muestra &
                    dfWi$isolate_type==input$isolate_type &
                    dfWi$num_pollos%in%input$num_pollos &
                    dfWi$dia%in%input$dia &
                    dfWi$treat_full%in%input$treat_full),]
  })
  
  #normalized wide data for use in adonis models
  DataNorm <-reactive({
   # trnasform 0-1 to numeric and delete missing vars
   dfWi[MarkVarList]=as.data.frame(lapply(MarkVarList, function(x)
    as.numeric(as.matrix(dfWi[x]))
   ))    
   dfWi[MarkVarList]=na.omit(dfWi[MarkVarList]) #remove missing variables
  
   if (input$ordi_var == "Phenotypic") {
       #if euclidean, scale the zones of inhibitions
       dfWi[TestDrugs]=as.data.frame(lapply(TestDrugs, function(x)
        scale(dfWi[x],center=FALSE,scale=TRUE)))
   }
   #if genotypic ordination, clean and collapse at sample level
   if (input$ordi_var != "Phenotypic") {
    #if not euclidean then:
    dfWi=dfWi[which(dfWi$phat!="0"),] #exclude isolates that were not genotypes
    #collapse at sample level for following varlist
    FactorVars=c("trial","dia", "sample_id", "tipo_muestra", "isolate_type", "num_pollos", "treat_full", "litter","int1", "tetAB", "qnrb")
    dfWi=aggregate(dfWi[MarkVarList], by=dfWi[which(names(dfWi)%in%FactorVars)],FUN=sum) 
   }
    DataNorm=dfWi[which(dfWi$trial%in%input$trial &
                         dfWi$tipo_muestra %in% input$tipo_muestra &
                         dfWi$isolate_type==input$isolate_type &
                         dfWi$num_pollos%in%input$num_pollos &
                         dfWi$dia%in%input$dia &
                         dfWi$treat_full%in%input$treat_full),] 
  })
  
  
  
 # TAB1: Generate a kernel density plot f
 output$DensPlot <- renderPlot({
  my.plotDens(Data(), x1(),y, y_cat(), by_var, random.eff, 
              title=paste(c("Zone kernel density, by"),x1()
                          ))
 }, width = 500, height = 500)
 
 # Generate bar graph for single-factor
 output$BarPlotSingle<- renderPlot({
  #Switch title strings for the y_cat variable
  TitleStr <- switch(input$y_cat, 
                     "Res" = "Percent resistant",
                     "int1" = "int1",
                     "tetAB" = "tetA and tetB",
                     "qnrb" = "qnrB")
  
  my.barChart.single(Data(),x1(),y_cat(),by_var,random.eff, 
                     title=paste(c(TitleStr,", by"),x1()
                                 ))
 }, width = 500, height = 500)
 
 # Generate bar graph for 2-factor analysis
 output$BarPlotDouble<- renderPlot({
  #Switch title strings for the y_cat variable
  TitleStr <- switch(input$y_cat, 
                     "Res" = "Percent resistant",
                     "int1" = "int1",
                     "tetAB" = "tetA and tetB",
                     "qnrb" = "qnrB")
  
  my.barChart.double(Data(),x1(),x2(),y_cat(),by_var,random.eff, 
                     title=paste(c(TitleStr,", by"),x1(),"and", x2()
                                 ),
                     stacked())
 },width = 500, height = 500)
 
 #Tab 2.1 univar modesl-----
 output$ModelSum <- renderPrint({
  # Take a dependency on input$goButton
  input$modelButton
  
  # Use isolate() to avoid dependency on input$obs
  model <- isolate(
   my.mod.list(df=Data(), x=x1(), y, by_var, random.eff, random.effBig, y_cat(), nest_mod=TRUE,always_nest=nest()))

  tab1=ldply(model,function(z){
   if(!is.na(z)) c(z$Mod$Tab, z$ModCat$Tab)
  })
  
  if (input$class_ag==TRUE) { #change the input for Drug() dropdown if aggregating by class
   Drug<- switch(input$Drug, 
                      "Ampicillin"= "Beta-lactam", "Amoxicillin/clavulanate"= "Beta-lactam",
                        "Cefotaxime"= "Beta-lactam", "Cephalothin" = "Beta-lactam",
                      "Ciprofloxacin"="Fluoroquinolone", "Enrofloxacin" = "Fluoroquinolone",
                      "Streptomycin" = "Aminoglycoside", "Gentamicin" = "Aminoglycoside",
                      "Sulfisoxazole" = "Sulfonamide", "Trimethoprim" = "Sulfonamide")
  }
  if(class(model[[Drug]])=="list") {
   modDrugSum=try(model[[Drug]]$ModCat$SumCat)
   modDrugSum$SumCat$call=modDrugSum$SumCat$call[2]
   list(tab1, modDrugSum)  }
  else{
   list(tab1) 
  }
    })

output$ModelDets <- renderPrint({
 # Take a dependency on input$goButton
 input$modelButton
 })

#Tab 2.2 bivar graphs-----
output$SummaryTab <- renderPrint({
 df.sum.cat2(df=Data(),x1(),x2(),y_cat(),by_var,random.eff)
})



#subtab 2.1------
output$PhenGenPlot<- renderPlot({
 TitleStr3 <- switch(input$x1, 
                    "treat_full"="treatment group",
                    "treat_water"="treatment in water", 
                    "treat_feed"="treatment in food", 
                    "dia"="day",
                    "num_pollos"="chickens in a cage", 
                    "tipo_muestra"="sample type" )
 #my.gen.phen(DataWide(),left_stack(),right_stack(),x1(),x1_val="01",random.eff)
 plots=gen.phen.arrange(DataWide(),left_stack(),right_stack(),x1(),random.eff,
                 title=paste(left_stack(),right_stack(), ",by",TitleStr3))
 #grid.layout(plots[[2]],plots[[1]],plots[[3]], nrow=1)
 plots[[2]]
 plots[[1]]
 plots[[3]]
},width = 1000, height = 500)

#subtab 2.2 ordination ------
output$OrdiPlot<- renderPlot({
 input$ordiButton
 isolate({
 ModType <- switch(input$ordi_var, 
                    "Phenotypic" = "euclidean",                
                    "Genotypic (jaccard)" = "jaccard",
                    "Genotypic (bray)" = "bray")
 #NMDS

  nmds.temp=metaMDS(as.matrix(DataNorm()[MarkVarList]), distance=ModType, 
                    autotransform =TRUE, 
                    wascores = TRUE, 
                    expand = TRUE, 
                    noshare=FALSE,
                    zerodist = "ignore", 
                    center=TRUE,
                    trymax=5)
  ordiplot(nmds.temp,display=c("species"), type=c("points"),xlim=range_x(),ylim=range_y(),
           main=paste("NMDS for", ordi_var(),"by",x1() ,"N=",nrow(DataNorm())))
  ordiellipse(nmds.temp, DataNorm()[[x1()]],draw="polygon",
              col=c("grey99"),label=TRUE,cex=0.6)
  })
 },width = 500, height = 500) 

# output$CommunitySim <- renderPrint({
#  #collapse community data by x1 and x2 levels (e.g., day and cage)
#  
#  #calculate richness and evenness
#  
# })
 #subtab 2.3 adonis mod ------
#adonis model: use full dataset! have checkbox of what variables to add
 
 output$AdonisMod <- renderPrint({

 input$adonisButton
 isolate({
  DepVar=paste(adonis_vars(), collapse="+")
  Formula=paste("DataNorm()[MarkVarList]~",DepVar)
  argList=list()
  argList[["data"]]=DataNorm()
  argList[["formula"]]=as.formula(Formula)
  if(adonis_strat()==TRUE&is.element("trail", adonis_vars())==FALSE) argList[["strata"]]="trial"
  #call adonis and pass arguments
  adonisSum=do.call("adonis", argList)
  adonisSum$call[3]="dynamically chosen set"
  adonisSum
 })
 })


# subtab 3.3 mutivar model selection-------

 #create formula from slectized inputs

MVmod<-reactive({
 MainVars=paste0(unique(fixed_main(),pheno_main()), collapse="+")
 
 IntsInt=""
 IntsTet=""
 #IntsTime=""
 if(int1_int()==TRUE) IntsInt=paste0(pheno_main(), "*int1", collapse="+")
 if(tetAB_int()==TRUE) IntsTet=paste0(pheno_main(), "*tetAB", collapse="+")
 #if(time_int()==TRUE) IntsTime=paste0(pheno_main(), "*dia", collapse="+")
 DepVars=paste(MainVars, IntsInt,IntsTet,IntsTime, sep="+")
 Formula=as.formula(paste("Res~",DepVars, REff()))
 #calculate models 
 MVmod=my.MVmod.list(Data(),formula=Formula, fixed_vars=fixed_main(),by_var,random.eff,random.effBig,odds_ratios=TRUE)
 
})


output$DredgeGraph <- renderPlot({
 #compile FE graphs
 FEplots<-llply(MVmod(), function(z){
  if(class(z)=="list") z$FEPlot[[2]][[1]]
 })
 FEplots=FEplots[lapply(FEplots,length)>0]
 args.list <- c(FEplots,list(ncol=2),main="ORs forbest models")
 do.call(grid.arrange, args.list)
},width=800, height=1000)

#render outputs

output$DredgeMod<- renderPrint({
  out=llply(MVmod(), function(z){
   if(class(z)=="list") z$bestModSum
  })
  out
})




#Tab3 data summary --------

output$RawData <- renderDataTable({
 DataWide()
})
output$SampleCounts <- renderTable({
 ddply(Data(),tab_vars(),
       summarise, Chickens=length(unique(pollo_id)),Samples = length(unique(sample_id)),Isolates =length(unique(id)))
})

})
