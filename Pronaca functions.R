#FUNCTIONS for Pronaca data analysis 

#load packages:

#laod packages
libs <- c('perm','lme4','car','ez','MASS', 'scales', 'diptest','mixtools','arm','moments','ggdendro','rgeos','plyr', 'ggparallel','RColorBrewer', 'vegan','grid','gridExtra', 'sjPlot', 'MuMIn')
#lapply(libs, install.packages, character.only = T)
lapply(libs, require, character.only = T)
source("RenseNieFwdSel.R")
#################
#Statistical tests, one RE------
my.t.smallmix<- function(df,x,y,random.eff, y_cat, odds_ratios=TRUE) {
 #assign variables
 df$y <- df[[y]]
 df$x <- df[[x]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$y_cat<-df[[y_cat]]
 #continuos variable ANOVA with 1 random effect
 mix <-lmer(y ~ x + (1|random.eff),data=df,REML=FALSE)
 mix.sum=summary(mix)
 mix.anova = Anova(mix)
 
 #categorical variable: try lmer, if not use glmmPQL
 mix_cat <- try(lmer(y_cat ~ x + (1|random.eff), family=binomial,data=df),
                silent=TRUE)
 if (class(mix_cat) == "try-error" ){
  stop("No bueno")
 }
  else {
   mix_cat.sum=summary(mix_cat)
   if(odds_ratios==TRUE){
    mix_cat.sum$coefficients[,c(1,2)]=exp(mix_cat.sum$coefficients[,c(1,2)])
    attributes(mix_cat.sum$coefficients)$dimnames[[2]][1]="Odds Ratio"
   }
   mix_cat.anova=Anova(mix_cat)
   # assembling the output as list of two table
   out.mix <- prettyNum(c(mix.sum$coefficients[1,1],mix.sum$coefficients[2,1],mix.anova$P,mix.sum$AICtab[1],length(mix.sum$ngrps),mix_cat.anova$P,mix_cat.sum$AICtab[1],length(mix_cat.sum$ngrps)  ))               
   names(out.mix)= c("Intercept","Effect size", "GLMManovaP", "AICanova", "REsanova","GLMMlogitP", "AIClogit","REslogit") 
   return(list(Tab=out.mix, Sum = mix.sum, SumCat=mix_cat.sum, Mod=mix, ModCat=mix_cat))
  }
}


################

my.t.bigmixCross<- function(df,x,y,random.eff, random.effBig, y_cat, odds_ratios=TRUE) {
 #assign variables
 df$y <- df[[y]]
 df$x <- df[[x]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$random.effBig<-df[[random.effBig]]
 df$y_cat<-df[[y_cat]]
 #if(random.eff=="sample_id") df$random.eff=substr(df$random.eff,1,4) #this is tupid, but to make sure effects are truly crossed....
 
 #continuous variables: run anova with 2 crossed REs
 mix <-try(lmer(y ~ x + (1|random.eff) + (1|random.effBig) ,data=df,REML=FALSE),
           silent=TRUE)
 if (class(mix) == "try-error" ) {
  stop("No bueno")
 }
 
 else{
 mix.sum=summary(mix)
 mix.anova = Anova(mix)
 
 #categorical variable: lmer and flmmPQL with crossed REs
  mix_cat <- try(lmer(y_cat ~ x + (1|random.eff)+ (1|random.effBig), family=binomial,data=df),
                 silent=TRUE)
  if (class(mix_cat) == "try-error" ){
   stop("No bueno")
  }
   mix_cat.sum=summary(mix_cat)
   if(odds_ratios==TRUE){
    mix_cat.sum$coefficients[,c(1,2)]=exp(mix_cat.sum$coefficients[,c(1,2)])
    attributes(mix_cat.sum$coefficients)$dimnames[[2]][1]="Odds Ratio"
   }
   mix_cat.anova=Anova(mix_cat)
  # assembling the output as list of two table
   out.mix <- prettyNum(c(mix.sum$coefficients[1,1],mix.sum$coefficients[2,1],mix.anova$P,mix.sum$AICtab[1],length(mix.sum$ngrps),mix_cat.anova$P,mix_cat.sum$AICtab[1],length(mix_cat.sum$ngrps)  ))               
   names(out.mix)= c("Intercept","Effect size", "GLMManovaP", "AICanova", "REsanova","GLMMlogitP", "AIClogit","REslogit") 
  return(list(Tab=out.mix, Sum = mix.sum, SumCat=mix_cat.sum, Mod=mix, ModCat=mix_cat))
 }
 }



################

my.t.bigmixNest<- function(df,x,y,random.eff, random.effBig, y_cat, odds_ratios=TRUE) {
 #assign variables
 df$y <- df[[y]]
 df$x <- df[[x]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$random.effBig<-df[[random.effBig]]
 df$y_cat<-df[[y_cat]]
 
 #continuous variables: run anova with 2 crossed REs
 mix <-try(lmer(y ~ x + (1|random.effBig/random.eff) ,data=df,REML=FALSE),
           silent=TRUE)
 if (class(mix) == "try-error" ) {
  mix <-try(lmer(y ~ x + (1|random.eff) ,data=df,REML=FALSE),
            silent=TRUE)
  if (class(mix) == "try-error" ){
   mix <- lm(y_cat ~ x,data=df)
  }
 }
 mix.sum=summary(mix)
 mix.anova = Anova(mix)
 
 #categorical variable: lmer and flmmPQL with crossed REs
 mix_cat <- try(lmer(y_cat ~ x + (1|random.effBig/random.eff), family=binomial,data=df),
                silent=TRUE)
 if (class(mix_cat) == "try-error" ){
  #stop("No bueno")
  mix_cat=lmer(y_cat ~ x + (1|random.eff), family=binomial,data=df)
 }
 else {
  mix_cat.sum=summary(mix_cat)
  if(odds_ratios==TRUE){
   mix_cat.sum$coefficients[,c(1,2)]=exp(mix_cat.sum$coefficients[,c(1,2)])
   attributes(mix_cat.sum$coefficients)$dimnames[[2]][1]="Odds Ratio"
  }
  mix_cat.anova=Anova(mix_cat)
  # assembling the output as list of two table
  out.mix <- c(mix.sum$coefficients[1,1],mix.sum$coefficients[2,1],mix.anova$P,mix.sum$AICtab[1],length(mix.sum$ngrps),mix_cat.anova$P,mix_cat.sum$AICtab[1],length(mix_cat.sum$ngrps)  )               
  names(out.mix)= c("Intercept","Effect size", "GLMManovaP", "AICanova", "REsanova","GLMMlogitP", "AIClogit","REslogit") 
  
  return(list(Tab=out.mix, Sum = mix.sum, SumCat=mix_cat.sum, Mod=mix, ModCat=mix_cat))
 }
}



#####################
#Stars----
my.star.mix <- function (p.table) {
 p.table$stars <- "" # to indicate no effect, looks cleaner in the plot
 p.table$stars[p.table$X<=  0.1]  <- "-"
 p.table$stars[p.table$X <=  0.05]  <- "*"
 p.table$stars[p.table$X <=  0.01]  <- "**"
 p.table$stars[p.table$X <= 0.001]  <- "***"
 p.text=paste(p.table$stars)
 return(p.text)
}


#####################
#TABLE OF MODELS----

###DATA FOR LABELS


#data for labelling continuous values (P-values from NLM model)
#UNIVARIATE MODEL SELECTION-----------
my.mod.list=function(df, x, y, by_var, random.eff, random.effBig, y_cat, nest_mod=FALSE, always_nest=FALSE){
 df$y <- df[[y]]
 df$x <- df[[x]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$random.effBig<-df[[random.effBig]]
 df$y_cat<-df[[y_cat]]
 
 out=dlply(df, .(by_var),
               function(z) {
                REMod=try(my.t.smallmix(df=z, x, y, random.eff, y_cat) ,silent=TRUE)
                if(class(REMod)=="try-error") return(NA)
                else{
                RECrossMod=try(my.t.bigmixCross(df=z, x, y, random.eff, random.effBig, y_cat)
                    ,silent=TRUE)
                if(class(RECrossMod)=="try-error") return(NA)
                else{
                if(nest_mod==TRUE) {
                 RENestMod=try(my.t.bigmixNest(df=z, x, y, random.eff, random.effBig, y_cat)
                               ,silent=TRUE)
                 if(class(RENestMod)=="try-error") return(NA)
                 temp_modlist=list(NoREMod=REMod,RECrossMod=RECrossMod,RENestMod=RENestMod)
                 #compare the fit of 3 models for categorical and continuous predictors separately
                 anovaModCat=anova(REMod$ModCat, RECrossMod$ModCat, RENestMod$ModCat)
                 anovaMod=anova(REMod$Mod, RECrossMod$Mod, RENestMod$Mod)
                }
                if(nest_mod==FALSE) {
                 temp_modlist=list(NoREMod=REMod,RECrossMod=RECrossMod)
                 #compare the fit of 3 models for categorical and continuous predictors separately
                 anovaModCat=anova(REMod$ModCat, RECrossMod$ModCat)
                 anovaMod=anova(REMod$Mod, RECrossMod$Mod)
                }

                # choose nested model for categorical,if non of chi-quared p-values are <0.05 pick it
                index_bestModCat=which(anovaModCat$Pr<0.05)
                #force the nested to be best if that option is chosen
                if(always_nest==TRUE) index_bestModCat=3
                
                if(length(index_bestModCat)>0) ModCat=temp_modlist[[index_bestModCat]][c("SumCat","Tab")]
                if(length(index_bestModCat)==0) ModCat=temp_modlist[[which(anovaModCat$AIC==min(anovaModCat$AIC))]][c("SumCat","Tab")]
                ModCat$Tab=ModCat$Tab[c("REslogit","GLMMlogitP")]
                
                #best model for non-categorica
                index_bestMod=which(anovaMod$Pr<0.05)
                if(always_nest==TRUE) index_bestMod=3
                
                if(length(index_bestMod)>0) Mod=temp_modlist[[index_bestMod]][c("Sum","Tab")]
                if(length(index_bestMod)==0) Mod=temp_modlist[[which(anovaMod$AIC==min(anovaMod$AIC))]][c("Sum","Tab")]
                Mod$Tab=Mod$Tab[c("Intercept","Effect size", "REsanova","GLMManovaP")]
               
                list(Mod=Mod, ModCat=ModCat)
                }
                }
                })
 #return(my.star.mix(out))
 return(out)
}

#MULTIVARIATE MODEL SELECTION LME4------

my.MVmod.list=function(df,formula, fixed_vars, by_var, random.eff,random.effBig, odds_ratios=TRUE){
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$random.effBig<-df[[random.effBig]]
 
 outlist <- sapply(levels(df$by_var),function(x) NULL)
 #had to use a loop b/c ddply kept giving errors!!!!
 options(na.action = "na.fail")
 for(drug in levels(droplevels(df$by_var))){
  df_sub=df[df$by_var==drug,]
  bestModSum=NA
  bestModSum<-try({
   gloMod <- glmer(as.formula(formula), data=df_sub, family=binomial,na.action = "na.fail")
   gloMod@call$formula <- formula
   gloMod@call$data <- df_sub
   dredgeProc=dredge(gloMod,fixed=fixed_vars)
   bestMod= get.models(dredgeProc, 1)[[1]]
   bestModSum= summary(bestMod)
   })
  if(class(bestModSum)=="try-error") outlist[[drug]]=NA
  else{
   #remove long call string, replace w title
    bestModSum$call$control=paste("fixed vars:", paste0(fixed_vars, 
                                                                       collapse=","), 
                                  "; initial model:",as.character(formula))
    bestModSum$call[3]=paste(toupper(drug))
    #convert to ORs
   if(odds_ratios==TRUE){
    bestModSum$coefficients[,c(1,2)]=exp(bestModSum$coefficients[,c(1,2)])
    attributes(bestModSum$coefficients)$dimnames[[2]][1]="Odds Ratio"
   }
   #plots
    FEPlot=sjp.glmer(bestMod,type="fe",title=paste("Fixed effects,",drug), 
                     fade.ns=TRUE,showPValueLabels=TRUE, printPlot=FALSE)
    REPlot=sjp.glmer(bestMod,type="re",title=paste("Random effects,",drug), 
                     fade.ns=TRUE,showPValueLabels=TRUE, printPlot=FALSE) 
   outlist[[drug]]$bestModSum <-bestModSum
   #outlist[[drug]]$allMods <-dredgeProc
   outlist[[drug]]$FEPlot <-FEPlot
   outlist[[drug]]$REPlot <-REPlot
  }
 }
 options(na.action = "na.omit")
 return(outlist)
}
#  
#  out=dlply(df, .(by_var),
#            function(z) {
#             bestModSum <- try({
#              gloMod<-glmer(formula, data=z, family=binomial,na.action = "na.fail")
#              dredgeProc=dredge(gloMod)
#              bestMod= get.models(dredgeProc, 1)[[1]]
#              summary(bestMod)  
#            # })
#              if(class(bestModSum)=="try-error") return(NA)
#            # else{ 
#              #get title string
#              TitleStr=as.character(unique(z[[by_var]]))
#              #convert to ORs
#              if(odds_ratios==TRUE){
#               bestModSum$coefficients[,c(1,2)]=exp(bestModSum$coefficients[,c(1,2)])
#               attributes(bestModSum$coefficients)$dimnames[[2]][1]="Odds Ratio"
#              }
#              #remove long call string, replace w title
#              bestModSum$call$control=toupper(TitleStr)
#              #plots
#              FEPlot=sjp.glmer(bestMod,type="fe",title=paste("Fixed effects,",TitleStr), 
#                             fade.ns=TRUE,showPValueLabels=TRUE, printPlot=FALSE)
#              REPlot=sjp.glmer(bestMod,type="re",title=paste("Random effects,",TitleStr), 
#                               fade.ns=TRUE,showPValueLabels=TRUE, printPlot=FALSE)
#              #return list
#              list(bestModSum=bestModSum, allMods=dredgeProc, FEPlot=FEPlot, REPlot=REPlot)
#            }
#            }
#            
#           )
#  return(out)
# }




#create table wiht position for labels
# df_lab = data.frame(x=rep(1.5,12), 
#                     y= rep(as.numeric(0.9*max(df$zone,na.rm=TRUE)),12), 
#                     drug_name= levels(df$drug_name),lab=my.star.mix(p.table))
# 
# #data for labelling categorical variables (%R and Chi-sq test)
df.sum.cat<- function(df, x1, y_cat, by_var, random.eff){
 df$x1 <- df[[x1]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$y_cat<-df[[y_cat]]
 
 out=(ddply(df, .(by_var, x1), summarise, 
                   n_sam=length(unique(random.eff)), 
                   n_isol=length(random.eff),
                   n_NS = sum(as.numeric(na.omit(y_cat))),
                   pct_R=100*n_NS/n_isol
                   #pct_R=paste0(round(100*n_NS/n_isol,digits=1),"%")
                   #x.pos=mean(BP_custom)
            ))
 return(out)
}

df.sum.cat2<- function(df, x1, x2, y_cat, by_var, random.eff){
 df$x1 <- df[[x1]]
 df$x2 <- df[[x2]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$y_cat<-df[[y_cat]]
 
 out=(ddply(df, .(by_var, x1, x2), summarise, 
            n_sam=length(unique(random.eff)), 
            n_isol=length(random.eff),
            n_NS = sum(as.numeric(na.omit(y_cat))),
            pct_R=100*n_NS/n_isol
            #pct_R=paste0(round(100*n_NS/n_isol,digits=1),"%")
            #x.pos=mean(BP_custom)
 ))
 return(out)
}
# 
# 
# #####################
# 
# df.pval.cat<- as.data.frame(
#  cbind(drug_name=levels(p.table$drug_name),
#        chi_lab=paste("                    P zone = "
#                      , round(p.table$X,3),'\n', 
#                      "                    P cat = ", 
#                      round(p.table$LogitP,3))))
# 
# 
# #merge the two
# df.lab.cat = unique(merge(x=df.sum.cat,
#                           y=df.pval.cat,
#                           by='drug_name', sort=FALSE, all=FALSE))
# 
# 
# ####GRAPH
# 
# 
# BAR CHARTS-------

my.barChart.single <- function(df, x1, y_cat, by_var, random.eff, title=""){
 df$x1 <- df[[x1]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$y_cat<-df[[y_cat]]
 
 df_coll=df.sum.cat(df,x1 ,y_cat ,by_var, random.eff)
 
 ggplot(data = df_coll, aes(x=x1 , 
            y = pct_R,
            label=paste(round(pct_R,1),"%"))
         )+
  geom_bar(stat = "identity",fill = "dark grey", colour = "black", alpha = 1/3)+
  geom_text( position = position_dodge(width=1), size=3, vjust=-0.5)+ 
  facet_wrap(~ by_var) +
  ylab("Percent of isolates") + 
  scale_y_continuous(limits = c(0,100))+
  theme_bw()+
  theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle(paste(title,
                " (N = ", sum(df_coll[1:length(unique(df_coll$x1)),]$n_isol)," isolates |",sum(df_coll[1:length(unique(df_coll$x1)),]$n_sam)," samples)")) 
}

######

my.barChart.double <- function(df, x1, x2, y_cat, by_var, random.eff, title="", bar_position="dodge"){
 df$x1 <- df[[x1]]
 df$x2 <- df[[x2]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$y_cat<-df[[y_cat]]
 
df_coll=df.sum.cat2(df,x1 ,x2 ,y_cat ,by_var, random.eff)
N_isol=sum(df_coll[1:
                    (length(unique(df_coll$x1))+length(unique(df_coll$x2)))]$n_isol)
N_sam=sum(df_coll[1:
                    (length(unique(df_coll$x1))+length(unique(df_coll$x2)))]$n_sam)

if (bar_position=="dodge"){
Plot=ggplot(aes(x1 , 
           y = 100*n_NS/n_isol, 
           fill=x2,
           label=paste(as.character(round(100*n_NS/n_isol,1)), "%")) ,
       data = df_coll)+
 geom_bar(stat = "identity",position="dodge") + 
 scale_y_continuous(limits = c(0,100))+
 ylab("Percent of isolates")
 
 #geom_text( position = position_dodge(width=1), size=3, angle=45)+
}
if(bar_position!="dodge"){
 Plot=ggplot(aes(x1 , 
                 y = n_NS, 
                 fill=x2,
                 label=paste(as.character(round(100*n_NS/n_isol,1)), "%")) ,
             data = df_coll)+
  geom_bar(stat = "identity", position=bar_position)+
  ylab("N")
}

Plot+
 facet_wrap(~ by_var) +
 theme_bw()+
 theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))+
 ggtitle(paste(title,
               " (N = ", N_isol ," isolates |",N_sam," samples)"))  

}

# 
# ##Density plot w annotations
my.plotDens <- function(df, x1, y, y_cat, by_var, random.eff, title=""){
 df$y <- df[[y]]
 df$x1 <- df[[x1]]
 df$by_var <- df[[by_var]]
 df$random.eff<-df[[random.eff]]
 df$y_cat<-df[[y_cat]]
 
 ggplot(df, aes(y, color = x1)) +
  geom_density(alpha=0.2) +
  #facet_wrap(~ by_var,scales="free_y") +
  facet_wrap(~ by_var,scales="free_x")+
  theme_bw()+
  theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(aes(xintercept = S_bpNew, linetype="dashed",color="grey"),linetype="dashed",color="grey")+
  geom_vline(aes(xintercept = R_bpCust, linetype="dashed",color="blue"),linetype="dashed",color="blue")+
  ggtitle(paste(title," (N = ", length(unique(paste(df$random.eff,df$isolate_code)))," isolates |",length(unique(df$random.eff))," samples)"))
}
#####################
#Gen-phen graph-----
#takes the wide data
#fix the ggparallel function 
source("my.ggparallel.R")

my.gen.phen<- function(df, left_stack, right_stack, x1, x1_val,  random.eff){
 df$right_stack <- df[[right_stack]]
 df$left_stack <- df[[left_stack]]
 df$x1 <- df[[x1]]
 df$random.eff<-df[[random.eff]]
 
 df=subset(df, df$x1==x1_val)
 len1=length(unique(df$left_stack))
 len2=length(unique(df$right_stack))
 label1 = paste("Richness=",len1)
 label2 = paste("Richness=",len2)
 
 IsolateN=nrow(df)
 SampleN=length(unique(df$random.eff))
 titleN= paste(x1_val,"(",IsolateN,"isolates |", SampleN," samples)")
  
 plot=my.ggparallel(list(left_stack, right_stack), data=df,text.angle=0, color="white")+
  scale_color_manual(values = c( clrs.hcl(n=len1+len2), rep("grey80", len1+len2)), guide="none")+
  annotate("text", 
           label = label1, x = 1, y = IsolateN +5 , size = 3, colour = "gray25")+
  annotate("text", 
           label = label2, x = 2, y = IsolateN +5 , size = 3, colour = "gray25") +
  ggtitle(titleN)+ 
  ylab("Isolates")+
  theme_bw()+
  theme(legend.position="none")
 return(plot)
}

gen.phen.arrange<-function(df, left_stack, right_stack, x1,  random.eff, title="" ){

 df$right_stack <- df[[right_stack]]
 df$left_stack <- df[[left_stack]]
 df$x1 <- df[[x1]]
 df$random.eff<-df[[random.eff]]
 
 plots<-llply(unique(as.character(df$x1)), function(x1_val){
  my.gen.phen(df, left_stack, right_stack, x1, x1_val,  random.eff)
 })
 names(plots)=unique(as.character(df$x1))
 #just print ou the 3 plots... 
 args.list <- c(plots,list(nrow=1),main=title)
 do.call(grid.arrange, args.list)
 
}


#####Brewer colors------
clrs.brew.spec <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

###Brewer color palette, default is blues
clrs.brew<- function(n=9,scheme="Blues"){
 colorRampPalette(brewer.pal(9,scheme))(n) 
}
clrs.hcl <- function(n) {                                                                                                                          hcl(h = seq(230, 0, length.out = n),                                                                                                                               c = 60, l = seq(10, 90, length.out = n),                                                                                                                             fixup = TRUE)                                                                                                                              }