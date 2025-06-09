###PLOTS

###PLOTS###
library(ggplot2)
library(tidyverse)
library(igraph)
library(Cairo)

###FIGURE MANUSCRIPT
###FIGURAS DELTAS
###################################
###Delta metrics= initial - final###
####################################
#load data and combine tables

#load data and combine tables
setwd("G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1")
#setwd("C:\\Users\\User\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")

SC1<-read.table("tabla_final_SC1.csv", sep=",", header=T) 
scenario=rep("SC1",nrow(SC1))
SC11<-cbind(scenario,SC1)

SC2<-read.table("tabla_final_SC2.csv", sep=",", header=T) 
scenario=rep("SC2",nrow(SC2))
SC21<-cbind(scenario,SC2)

SC3<-read.table("tabla_final_SC3.csv", sep=",", header=T) 
scenario=rep("SC3",nrow(SC3))
SC31<-cbind(scenario,SC3)

tabla_final<-rbind(SC11,SC21,SC31)

##########
#DELTAS CONTROLES
no_resistance<-read.table("control_no_resistance.csv", sep=",", header=T) 
scenario=rep("no_resistance",nrow(no_resistance))
no_resistance1<-cbind(scenario,no_resistance) 

no_dyn_SC1<-read.table("CONTROL_NO_DYNAMIC_SC1.csv", sep=",", header=T) 
scenario=rep("no_dyn_SC1",nrow(no_dyn_SC1))
no_dyn_SC11<-cbind(scenario,no_dyn_SC1)

no_dyn_SC2<-read.table("CONTROL_NO_DYNAMIC_SC2.csv", sep=",", header=T) 
scenario=rep("no_dyn_SC2",nrow(no_dyn_SC2))
no_dyn_SC21<-cbind(scenario,no_dyn_SC2)

no_dyn_SC3<-read.table("CONTROL_NO_DYNAMIC_SC3.csv", sep=",", header=T) 
scenario=rep("no_dyn_SC3",nrow(no_dyn_SC3))
no_dyn_SC31<-cbind(scenario,no_dyn_SC3)

no_modif_SC1<-read.table("CONTROL_NO_MODIFIERS_SC1.csv", sep=",", header=T) 
scenario=rep("no_modif_SC1",nrow(no_modif_SC1))
no_modif_SC11<-cbind(scenario,no_modif_SC1)

no_modif_SC2<-read.table("CONTROL_NO_MODIFIERS_SC2.csv", sep=",", header=T) 
scenario=rep("no_modif_SC2",nrow(no_modif_SC2))
no_modif_SC21<-cbind(scenario,no_modif_SC2)

no_modif_SC3<-read.table("CONTROL_NO_MODIFIERS_SC3.csv", sep=",", header=T) 
scenario=rep("no_modif_SC3",nrow(no_modif_SC3))
no_modif_SC31<-cbind(scenario,no_modif_SC3)

no_rew_SC1<-read.table("CONTROL_NO_REWIRING_SC1.csv", sep=",", header=T) 
scenario=rep("no_rew_SC1",nrow(no_rew_SC1))
no_rew_SC11<-cbind(scenario,no_rew_SC1)

no_rew_SC2<-read.table("CONTROL_NO_REWIRING_SC2.csv", sep=",", header=T) 
scenario=rep("no_rew_SC2",nrow(no_rew_SC2))
no_rew_SC21<-cbind(scenario,no_rew_SC2)

no_rew_SC3<-read.table("CONTROL_NO_REWIRING_SC3.csv", sep=",", header=T) 
scenario=rep("no_rew_SC3",nrow(no_rew_SC3))
no_rew_SC31<-cbind(scenario,no_rew_SC3)


tabla_final_controles= rbind(SC11,SC21,SC31,no_resistance1, no_dyn_SC11, no_dyn_SC21, no_dyn_SC31,no_modif_SC11,no_modif_SC21,no_modif_SC31,no_rew_SC11,no_rew_SC21,no_rew_SC31 )

tabla_final=tabla_final_controles #la llamo tabla final para no cambiar todo el codigo


###"S_plantas", "S_poli"
#"links", "evenn", "conn","nest","robust", "funct_comp" "N_overHL" "N_overLL" "int_div" "gen" "vul" "H2" 
tabla_final<-na.omit(tabla_final)

##ESTIMATE DELTAS
diff<-list()
for (z in 1:13){
	sc=unique(tabla_final$scenario)
	escenario=tabla_final[tabla_final$scenario==sc[z],]
	delta=list()
	for (k in 1:1000){
		ciclo=escenario[escenario$cycle==k,]
		inicial=ciclo[ciclo$Sumulation==1,]
		final=ciclo[ciclo$Sumulation==16,]
		if (dim(final)[1]==0){diffe = rep(NA,16)}
		if (dim(final)[1]>0){diffe = as.numeric(final[c(2:17)])- as.numeric(inicial[c(2:17)])}
		delta[[k]]= c(sc[z],k,diffe)
	}
	ddd=do.call(rbind,delta)
	diff[[z]]= ddd
}

diferencias=do.call(rbind, diff)
colnames(diferencias)=c("SC","cycle","S_plantas","S_poli","links","evenn","conn","nest", "nodf","robust","robustPl","funct_comp","N_overHL", "N_overLL", "int_div", "gen", "vul", "H2")
write.table(diferencias, "deltas_todos1.csv", sep=",",row.names=F)


###############################################################################################

###FIGURE 2
#setwd("")
#library(emmeans) # v. 1.7.0
#library(magrittr) # v. 2.0.1
#library(multcomp)
setwd("G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1")
aa=read.table("deltas_todos.csv", sep=",", header=T)

colore<-c("lightgoldenrod3","goldenrod3","indianred3")#,"grey50" )
aa1=aa[aa$SC=="1"|aa$SC=="2"|aa$SC=="3",]###
par(mfrow=c(3,2))
AA<-boxplot(aa1$S_plantas~aa1$SC, ylab=expression(Delta~" Plant richness"), xlab=NULL, cex.lab=1.1, col=colore, names = c("HP1", "HP2", "HP3"), las=1)
BB<-boxplot(aa1$S_poli~aa1$SC, ylab=expression(Delta~" Pollinator richness"), xlab=NULL,cex.lab=1.1, las =1, col=colore, names = c("HP1", "HP2", "HP3"))
DD<-boxplot(aa1$evenn~aa1$SC, ylab=expression(Delta~" Interaction evenness"), xlab=NULL,cex.lab=1.1, las=1,col=colore, names = c("HP1", "HP2", "HP3"))
EE<-boxplot(aa1$conn~aa1$SC, ylab=expression(Delta~" Connectance"),xlab=NULL, cex.lab=1.1, las=1,col=colore, names = c("HP1", "HP2", "HP3"))
FF<-boxplot(aa1$nest~aa1$SC, ylab=expression(Delta~" Nestedness (NODF)"),xlab="Herbicide application programmes", cex.lab=1.1, las=1,col=colore, names = c("HP1", "HP2", "HP3"))
GG<-boxplot(aa1$robustPl~aa1$SC, ylab=expression(Delta~" Robustness"),xlab="Herbicide application programmes", cex.lab=1.1,las=1, col=colore, names = c("HP1", "HP2", "HP3"))

#4= no_resistance1,5= no_dyn_SC11, 6=no_dyn_SC21, 7=no_dyn_SC31,8=no_modif_SC11,9=no_modif_SC21,10=no_modif_SC31,11=no_rew_SC11,12=no_rew_SC21,13=no_rew_SC31
####delta suplementarias control
aa1=aa[aa$SC=="1"|aa$SC=="2"|aa$SC=="3"|aa$SC=="11"|aa$SC=="12"|aa$SC=="13",]###
colore<-c("lightgoldenrod3","goldenrod3","indianred3","grey90","grey70","grey50" )
par(mfrow=c(3,2))
AA<-boxplot(aa1$S_plantas~aa1$SC, ylab=expression(Delta~" Plant richness"), xlab=NULL, cex.lab=1.1, col=colore, names = c("HP1", "HP2", "HP3", "NO REW1", "NO REW2","NO REW3"))
BB<-boxplot(aa1$S_poli~aa1$SC, ylab=expression(Delta~" Pollinator richness"), xlab=NULL,cex.lab=1.1, col=colore, names = c("HP1", "HP2", "HP3","NO REW1", "NO REW2","NO REW3"))
DD<-boxplot(aa1$evenn~aa1$SC, ylab=expression(Delta~" Interaction evenness"), xlab=NULL,cex.lab=1.1, col=colore, names = c("HP1", "HP2", "HP3","NO REW1", "NO REW2","NO REW3"))
EE<-boxplot(aa1$conn~aa1$SC, ylab=expression(Delta~" Connectance"),xlab=NULL, cex.lab=1.1, col=colore, names = c("HP1", "HP2", "HP3","NO REW1", "NO REW2","NO REW3"))
FF<-boxplot(aa1$nest~aa1$SC, ylab=expression(Delta~" Nestedness (NODF)"),xlab="HPs and control models", cex.lab=1.1, col=colore, names = c("HP1", "HP2", "HP3","NO REW1", "NO REW2","NO REW3"))
GG<-boxplot(aa1$robustPl~aa1$SC, ylab=expression(Delta~" Robustness"),xlab="HPs and control models", cex.lab=1.1, col=colore, names = c("HP1", "HP2", "HP3", "NO REW1", "NO REW2","NO REW3"))


###################################3
####Figura inicial final
tabla_final$S_poli<-replace(tabla_final[,3], tabla_final[,2]=="small", NA)
tabla_final$S_plantas<-as.numeric(replace(tabla_final[,2], tabla_final[,2]=="small", NA))

#tabla_final 
colore<-c("lightgoldenrod3", "lightgoldenrod3", "goldenrod3", "goldenrod3", "indianred3","indianred3")

par(mfrow=c(3,2))
AA<-boxplot(tabla_finalA$S_plantas ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6), ylab="Plant richness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

BB<-boxplot(tabla_finalA$S_poli ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6), ylab="Pollinator richness", xlab="Herbicide application programmes",xaxs = FALSE,cex.lab=1.1, col=colore, names = c(expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

CC<-boxplot(tabla_finalA$evenn ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6), ylab="Interaction evenness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1,col=colore, names = c(expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

DD<-boxplot(tabla_finalA$conn ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6), ylab="Connectance", xlab="Herbicide application programmes",xaxs = FALSE,cex.lab=1.1, col=colore, names = c(expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

EE<-boxplot(tabla_finalA$nest ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6), ylab="Nestedness (NODF)", xlab="Herbicide application programmes",xaxs = FALSE,cex.lab=1.1, col=colore, names = c(expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

FF<-boxplot(tabla_finalA$robustPl ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6), ylab="Robustness", xlab="Herbicide application programmes",xaxs = FALSE,cex.lab=1.1, col=colore, names = c(expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

####################controles


#4= no_resistance1,5= no_dyn_SC11, 6=no_dyn_SC21, 7=no_dyn_SC31,8=no_modif_SC11,9=no_modif_SC21,10=no_modif_SC31,11=no_rew_SC11,12=no_rew_SC21,13=no_rew_SC31
#SC1"           "SC2"           "SC3"        
#   "no_resistance"
 #[5] "no_dyn_SC1"    "no_dyn_SC2"    "no_dyn_SC3"   
 # "no_modif_SC1" 
 #[9] "no_modif_SC2"  "no_modif_SC3"  
 #"no_rew_SC1"    "no_rew_SC2"   
#[13] "no_rew_SC3"

tabla_final$S_poli<-replace(tabla_final[,3], tabla_final[,2]=="small", NA)
tabla_final$S_plantas<-as.numeric(replace(tabla_final[,2], tabla_final[,2]=="small", NA))
tabla_finalB=tabla_final[tabla_final$scenario=="SC1"|tabla_final$scenario=="SC2"|tabla_final$scenario=="SC3"|tabla_final$scenario=="no_rew_SC1"|tabla_final$scenario=="no_rew_SC2"|tabla_final$scenario=="no_rew_SC3",]#
tabla_finalA = tabla_finalB %>% filter(Sumulation == c(1,16))

colore<-c( "grey50", "grey50","grey70", "grey70","grey90", "grey90","lightgoldenrod3", "lightgoldenrod3", "goldenrod3", "goldenrod3", "indianred3","indianred3")

par(mfrow=c(3,2))
AA<-boxplot(tabla_finalA$S_plantas ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6,7:8,9:10,11:12), ylab="Plant richness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NRew1[i]), expression(NRew1[f]),expression(NRew2[i]), expression(NRew2[f]),expression(NRew3[i]), expression(NRew3[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

BB<-boxplot(tabla_finalA$S_poli ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8,9:10,11:12), ylab="Pollinator richness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NRew1[i]), expression(NRew1[f]),expression(NRew2[i]), expression(NRew2[f]),expression(NRew3[i]), expression(NRew3[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))


CC<-boxplot(tabla_finalA$evenn ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8,9:10,11:12), ylab="Interaction evenness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NRew1[i]), expression(NRew1[f]),expression(NRew2[i]), expression(NRew2[f]),expression(NRew3[i]), expression(NRew3[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))


DD<-boxplot(tabla_finalA$conn ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8,9:10,11:12), ylab="Connectance", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NRew1[i]), expression(NRew1[f]),expression(NRew2[i]), expression(NRew2[f]),expression(NRew3[i]), expression(NRew3[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

EE<-boxplot(tabla_finalA$nest ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8,9:10,11:12), ylab="Nestedness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NRew1[i]), expression(NRew1[f]),expression(NRew2[i]), expression(NRew2[f]),expression(NRew3[i]), expression(NRew3[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

FF<-boxplot(tabla_finalA$robustPl ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8,9:10,11:12), ylab="Robustness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NRew1[i]), expression(NRew1[f]),expression(NRew2[i]), expression(NRew2[f]),expression(NRew3[i]), expression(NRew3[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))


######
tabla_finalB=tabla_final[tabla_final$scenario=="SC1"|tabla_final$scenario=="SC2"|tabla_final$scenario=="SC3"|tabla_final$scenario=="no_resistance",]#
tabla_finalA = tabla_finalB %>% filter(Sumulation == c(1,16))

colore<-c( "grey50", "grey50","lightgoldenrod3", "lightgoldenrod3", "goldenrod3", "goldenrod3", "indianred3","indianred3")

par(mfrow=c(3,2))
AA<-boxplot(tabla_finalA$S_plantas ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6,7:8), ylab="Plant richness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NR[i]), expression(NR[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

BB<-boxplot(tabla_finalA$S_poli ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8), ylab="Pollinator richness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NR[i]), expression(NR[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

CC<-boxplot(tabla_finalA$evenn ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8), ylab="Interaction evenness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NR[i]), expression(NR[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

DD<-boxplot(tabla_finalA$conn ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8), ylab="Connectance", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NR[i]), expression(NR[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

EE<-boxplot(tabla_finalA$nest ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8), ylab="Nestedness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NR[i]), expression(NR[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))

FF<-boxplot(tabla_finalA$robustPl ~ tabla_finalA$Sumulation + tabla_finalA$scenario ,at = c(1:2, 3:4, 5:6, 7:8), ylab="Robustness", xlab="Herbicide application programmes",xaxs = FALSE, cex.lab=1.1, col=colore, names = c(expression(NR[i]), expression(NR[f]),expression(HP1[i]),expression(HP1[f]), expression(HP2[i]), expression(HP2[f]),expression(HP3[i]),expression(HP3[f])))




##FIGURE SUPPLEMENTARY
#First, it is necessary to summarize the data. This can be done in a number of ways, as described on this page. In this case, we’ll use the summarySE() function defined on that page, and also at the bottom of this page. (The code for the summarySE function must be entered before it is called here).
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}






###############}
###REPLACING "SMALL" WITH 1 (PLANTS AND POLLI)
#tabla_final$S_poli<-replace(tabla_final[,3], tabla_final[,2]=="small", 1)
#tabla_final$S_plantas<-as.numeric(replace(tabla_final[,2], tabla_final[,2]=="small", 1))


###"S_plantas", "S_poli"
tabla_final$S_plantas=as.numeric(tabla_final$S_plantas)
na.omit(tabla_final$S_plantas)
#summary for each scenario and simulation
tgc <- summarySE(tabla_final, measurevar=na.omit("S_plantas"),groupvars=c("scenario","Sumulation"))###, groupvars=c("supp","dose"))
#tgc

simulaciones=as.numeric(tgc$Sumulation)-1
tgc2<-cbind(tgc,simulaciones)

#Line graphs
#After the data is summarized, we can make the graph. These are basic line and point graph with error bars representing either the standard error of the mean, or 95% confidence interval.

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.15) # move them .05 to the left and right

##A finished graph with IC of the mean might look like this. The points are drawn last so that the white fill goes on top of the lines and error bars.
#colore<-c("lightgoldenrod3","goldenrod3","indianred3")
AA<-ggplot(tgc2, aes(x=as.numeric(simulaciones), y=S_plantas,colour=scenario)) + 
    geom_errorbar(aes(ymin=S_plantas-ci, ymax=S_plantas+ci), width=.2, position=pd)+#, colour="black") +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=16) + # 21 is filled circle, fill="white")
	scale_color_manual(values = c("SC1" = "lightgoldenrod3", "SC2" = "goldenrod3", "SC3"  = "indianred3"))+#, "SINMODIFSC1"= "olivedrab2", "SINMODIFSC2" = "olivedrab3","SINMODIFSC3" = "olivedrab4", "NODYN" = "black", "NODYN_META" = "grey50"))+#, "NORES"="black"))+
	xlab("time (two years interval)") +
    ylab("Pollinator richness") +
    #scale_colour_hue(name="Scenario of herbicides application",l=50) +    # Legend label, use darker color                     # Use darker colors, lightness=40
    #ggtitle("with family parameter - lambda=0.001") +
    expand_limits(y=0) +                        # Expand y range
    #scale_y_continuous(breaks=0:25*2) +         # Set tick every 4
	scale_x_continuous(breaks=0:15*2) + 
    #theme_bw(panel.grid = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))
	theme(legend.justification=c(1,0),
          legend.position=c(1,0.8))               # Position legend in bottom right



## para estos elimino los na primero
#"links", "evenn", "conn","nest","robustPl", "funct_comp" "N_overHL" "N_overLL" "int_div" "gen" "vul" "H2" 
#tabla_final<-na.omit(tabla_final[tabla_final$links!="NA",])
tabla_final<-na.omit(tabla_final)##pero esto me va a eliminar las small también.Me conviene usar para S_poli y S_plantas la tabla original y luego esta

#summary for each scenario and simulation
tgc <- summarySE(tabla_final, measurevar="robustPl",groupvars=c("scenario","Sumulation"))###, groupvars=c("supp","dose"))
#tgc

simulaciones=as.numeric(tgc$Sumulation)-1
tgc2<-cbind(tgc,simulaciones)

#Line graphs
#After the data is summarized, we can make the graph. These are basic line and point graph with error bars representing either the standard error of the mean, or 95% confidence interval.

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.15) # move them .05 to the left and right

##A finished graph with IC of the mean might look like this. The points are drawn last so that the white fill goes on top of the lines and error bars.
#colore<-c("lightgoldenrod3","goldenrod3","indianred3")
NN<-ggplot(tgc2, aes(x=as.numeric(simulaciones), y=robustPl,colour=scenario)) + 
    geom_errorbar(aes(ymin=robustPl-ci, ymax=robustPl+ci), width=.2, position=pd)+#, colour="black") +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=16) + # 21 is filled circle, fill="white")
	scale_color_manual(values = c("SC1" = "lightgoldenrod3", "SC2" = "goldenrod3", "SC3"  = "indianred3", "SINMODIFSC1"= "olivedrab2", "SINMODIFSC2" = "olivedrab3","SINMODIFSC3" = "olivedrab4"))+###c("SC1" = "lightgoldenrod3", "NORES" = "goldenrod3", "NOREW"  = "indianred3", "NODYN" = "black", "NODYN_META" = "grey50"))+#, "RAND" = "grey50", "NORES"="black"))+
	xlab("time (two years interval)") +
    ylab("robustness") +
    #scale_colour_hue(name="Scenario of herbicides application",l=50) +    # Legend label, use darker color                     # Use darker colors, lightness=40
    #ggtitle("Change on functional complementarity in 30 years of herbicides application") +
    expand_limits(y=0) +                        # Expand y range
    #scale_y_continuous(breaks=0:25*2) +         # Set tick every 4
	scale_x_continuous(breaks=0:15*2) + 
    #theme_bw(panel.grid = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))
	theme(legend.justification=c(1,0),
          legend.position=c(1,0.8))               # Position legend in bottom right





install.packages("gridExtra")
library("gridExtra")
grid.arrange(AA, BB, CC, EE, ncol=2, nrow =2)
grid.arrange(DD, FF, GG, HH, II, ncol=2, nrow =3)
grid.arrange(JJ, KK, LL, MM,  NN, ncol=2, nrow =3)


grid.arrange(II,JJ,KK,ncol=2, nrow =2)
grid.arrange(LL,MM,NN,ncol=2, nrow =2)



