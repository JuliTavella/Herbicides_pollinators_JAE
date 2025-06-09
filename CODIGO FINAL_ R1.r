####Modelling the effect of herbicide application on plant-pollinator networks in agroecosystems 

####Modelos con distintas opciones de dynamics and rewiring for pollinators
###Julia Tavella and Fred Windsor
##Preparing the environment
#install.packages(c("RColorBrewer","tidyverse","sads","bipartite","igraph","gplots","vegan","reshape"))
library(RColorBrewer)
library(tidyverse) 
library(sads) 
library(bipartite)
library(igraph)
library(gplots)
library(vegan)
library(reshape)

###First, some functions...
##list to matrix function
list2matr<-function(sp1,sp2,effects,matr.type="q"){
  nodelist1u<-array(unique(sp1))
  nodelist2u<-array(unique(sp2))
  m<-length(nodelist1u)
  n<-length(nodelist2u)
  matr<-matrix(0,m,n) #Matrix
  colnames(matr)<-nodelist2u   #######
  row.names(matr)<-nodelist1u
  for (i in 1:m){
    iu<-nodelist1u[i]
    ii<-which(sp1==iu)
    nodelist2ii<-sp2[ii]
    link.weightii<-effects[ii]
    for (j in 1:length(ii)){
      jj<-which(nodelist2u==nodelist2ii[j])
      if (matr.type=="b") matr[i,jj]<-1 #Binary matrix
      if (matr.type=="q") matr[i,jj]<-link.weightii[j] #Weighted matrix
    }
  }
  res<-list(matr)
  names(res)<-c("matr")
  return(res)
}


################
####Function for Plant community changes simulation##
plant_simulation<-function(init.plant.comm,duration,trait.losing,trait.expanding){#,resistance){ 
	#'''simulate losses and expansion of plant species
	#	input:a vector of an initial plant comunity, duration=number of cicles, vectors of probabilities of loss, expansion and resistance preasure
	#	output:an array containing the community composition at each time step (duration)
	#	'''
	
	## Storage arrays ####
	plant.comm_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for the data generated from the model
	plant.loss_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for coverage loss data generated from the model
	plant.exp_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for coverage expansion data generated from the model

	for (j in 1:duration){ # loop through each time step
								
		# Select the starting plant community
		if(j<=1){plant.comm=init.plant.comm} # set the plant community as the starting coverage if the time step is 1
		if(j>1){plant.comm=plant.comm_array[j-1,]} # if not, use the plant community from the previous time step
		
		##ESTE NUNCA TIENE QUE SER NAN		
		# calculate the loss of plant coverage due to herbicide application and resistance  #
		plant.loss_array[j,]<-plant.comm-(plant.comm*trait.losing)#lo que queda de la comunidad anterior
		
	    #set a threshold, a sensible proportion at wich we would not expect a plant population that covers <0.01 of the transect will be viable
		plant.loss_array[j,][plant.loss_array[j,]<0.001]<-0 ##fernando recomends 0.1 o 0.5. #!!! We have to set a real threshold here!
		#if 1 is the space we have, 0.01 is the 1% of area occupied by the plant. estoy usando 0.001
		
		space<-1-sum(plant.loss_array[j,]) # calculate the space available due to the loss of plants
			
		# calculate the proportional ability of plant species to expand into the space 
		# before that, I eliminate plant species that we lose at the previous step
		surviving.plants=plant.loss_array[j,] 
		
		surviving.plants[surviving.plants>0]<-1 #binary vector to "delete" lost species from the expanding vector 
		expanding1=trait.expanding * surviving.plants ## parameter of expansion 
		
		if (sum(expanding1)>0){plant.exp.ability<-(expanding1)/sum(expanding1)}
		if (sum(expanding1)==0){plant.exp.ability<- expanding1} #to avoid NaNs
		
		#Calculate the expansion of plants (more resistant plants more expantion chances)	
		plant.exp_array[j,]<-space*plant.exp.ability ##
		
		# calculate the change in plant community coverage after loss and expansion
		if (sum(plant.exp_array[j,])>0){plant.comm_array[j,]<-(plant.loss_array[j,]+plant.exp_array[j,])/sum(plant.loss_array[j,]+plant.exp_array[j,])}
		if (sum(plant.exp_array[j,])==0){plant.comm_array[j,]<-(plant.loss_array[j,]+plant.exp_array[j,])}
	}
	
	colnames(plant.comm_array)<-names(init.plant.comm)
	return(plant.comm_array)  #### return an array with simulations equal to duration 
}


######NETWORK ASSEMBLY function to construct the initial network###
network_assemblyV2<-function(meta_web, nueva_comunidad, N_pol){ 
	#"""Function to generate networks using an initial plant community, and 
	#based on the links with pollinators from the meta-web
	#I used mgen function (in bipartite) to genrate a random network using a probability matrix (based on pollinator frequency and plant species proportions)
	#input: meta_web, nueva_comunidad (initial plant community or plant community at each step of the simulation), N_links (mean number of links  
	#estimated from empirical networks by plot)
	#output: a plant-pollinator network based on the probability matrix and preserving the proportion of links 
	#""" 

	#select, previously, a subset of pollinators of size equal to the mean number of polli in empirical networks
	#give probabilities according to species frequency in meta_Web
	set.seed(NULL)
	polis=sample(1:ncol(meta_web),N_pol,prob=colSums(meta_web),replace=FALSE)
	nomb_polis=colnames(meta_web)[polis]
	
	#create a matrix containing plants from nueva_comunidad interacting with a the subset of pollinators, with values of interactions available in th meta-web
	sub_meta_web=matrix(0,length(nueva_comunidad), N_pol)#storage matrix
	rownames(sub_meta_web)=names(nueva_comunidad) 
	colnames(sub_meta_web)=nomb_polis
	
	for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
		sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta
		if (sp_planta %in% rownames(meta_web)){ #because not all plant species are in the meta_web
			for (j in 1:N_pol){
				sp_poli<- nomb_polis[j]
				sub_meta_web[i,j]=meta_web[sp_planta,sp_poli] #select the value of the interaction in meta_web to copy in sub_meta_web
			}
		}
	}
	
	#define pollinators relative frequency 
	suma_col=colSums(as.matrix(sub_meta_web))###
	pol.frec=suma_col/sum(suma_col)
		
	#define proportion of entomophilous plant species (deleting plants that not interact with pollinators from nueva_comunidad)
	nueva_comunidad2=nueva_comunidad[nueva_comunidad>0] 
	plant.prop=nueva_comunidad2/sum(nueva_comunidad2)

	#interaction probabilities matrix 
	sub_meta_web1=sub_meta_web
	sub_meta_web1[sub_meta_web1>=1]=1##forbidden links matrix.presence-absence links to give probabily 0 to interactions that do not exist in meta_web
	pmat=((plant.prop%*%t(pol.frec))*sub_meta_web1)/sum((plant.prop%*%t(pol.frec))*sub_meta_web1)#construct the matrix of probabilities, combining vectors and forbidden link mat
	rownames(pmat)=names(nueva_comunidad2)
	pmat= empty(pmat)#delete empty rows and columns if necesary
	
	#new random network based on pmat
	#number of links equal to sub_meta_web (proportion of links), all species remain in the network, quantitative network
	#Setting keep.species to FALSE may (but need not) cause a loss of species and networks have a small number of species. With keep.species=TRUE 
	#all selected species remain in the network. Results are quite different.
		set.seed(NULL)
		L=sum(sub_meta_web)*100
		nueva_red=mgen(pmat,n=L,keep.species=TRUE, rep.cell=TRUE, autotransform="sum")### Alternative:create the network usin "sample" function (is quite simmilar)
		rownames(nueva_red)=rownames(pmat)
		colnames(nueva_red)=colnames(pmat)
	
	return(nueva_red)
}


#changes in pollinators
pollinators_dynamV2<-function(nueva_comunidad, comunidad_anterior, red_anterior){
#'''function that simulates changes on abundance of pollinators species due to changes on plant community
	#input: nueva_comunidad=new plant community after simulation, comunidad_anterior=previous plant community (t(n-1)), 
	#red_anterior= previous network
	#output: a vector containing the new pollinator frequencies
	#'''
	#changes on plant proportions
	popul_changes=nueva_comunidad-comunidad_anterior #difference between previous and current proportions of each plant species
	#pollinators total frequency (num of links) in the network t(n-1)
	prev_pol=colSums(red_anterior)
	
	#dynamics on pollinator populations
	red_anterior[red_anterior>0]=1 #left 0 and 1 to multiply with popul changes
	
	new_ab=NULL###matrix(0,dim(red_anterior)[1],dim(red_anterior)[2])
	for (i in 1:length(prev_pol)){ #for each pollinator species in the previous network
		if (dim(red_anterior)[1]==1 & dim(red_anterior)[2]==1 ){
		hosts=rownames(red_anterior)
		}else{
		hosts=names(red_anterior[,i][red_anterior[,i]>0])		#select its host plants
		}
		if (sum(comunidad_anterior[hosts])>0){resources=sum(popul_changes[hosts])/sum(comunidad_anterior[hosts])}#estimates the rate of change in the proportion of pollinator resources
		if (sum(comunidad_anterior[hosts])==0){resources=sum(popul_changes[hosts])}
		
		if (resources>1){abun_change= sum(rbinom(prev_pol[i],1,prob=1))}#if resources increase, pollinators abundance increase (prob=1)
		if (resources< 0){abun_change= -(sum(rbinom(prev_pol[i],1,prob=abs(resources))))} #if we lost resources we substract abun_change (we loose pollin)
		if (resources >=0 && resources <=1){abun_change= sum(rbinom(prev_pol[i],1,prob=abs(resources)))} #if resources increase, the number of pollinators increase
		
		new_ab[i]=prev_pol[i] + abun_change
		}
	names(new_ab)=names(prev_pol)
	return(new_ab)
	}

	

#
#function to reassembly the network ###no hay rewiring, la nueva abundancia de polinators establece interacciones en la misma proporcion que en la red anterior

nueva_red<-function(nueva_comunidad, nuevas_abundancias, red_anterior){
	red_anterior_norm=sweep(red_anterior, 2, colSums(red_anterior), "/")#normalizo por columna (prop de c poli en sus plantas hospederas en la red anterior tn-1)
	
	new_mat=matrix(0,length(nueva_comunidad), length(nuevas_abundancias)) #empty matrix with the size of the new community 27/1/25
	rownames(new_mat)=names(nueva_comunidad) 
	colnames(new_mat)=names(nuevas_abundancias)
	
			for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
				sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta
				if (sp_planta %in% rownames(red_anterior)){ #because not all plant species are in the meta_web
				for (j in 1:length(nuevas_abundancias)){
					sp_poli<- as.character(names(nuevas_abundancias[j]))
					if (red_anterior[sp_planta,sp_poli]>0){
						new_mat[sp_planta,sp_poli]= round(red_anterior_norm[sp_planta,sp_poli]*nuevas_abundancias[j],0)
						}else{
						new_mat[sp_planta,sp_poli]= 0
						}
		        }
	            }
			}
		red_final=new_mat
		return(red_final)
}


new_reconectar_without_rewiring2<-function(nueva_comunidad, comunidad_anterior, red_anterior){
	if (is.array(red_anterior) && (dim(red_anterior)[1]) > 2 && (dim(red_anterior)[2])>2 && sum(nueva_comunidad>0)>2 && sum(comunidad_anterior>0)>2){
	#estimates new pollinators frequency 
		nuevas_abundancias=pollinators_dynamV2(nueva_comunidad, comunidad_anterior, red_anterior)
		if (length(nuevas_abundancias)>2){
		red_new=empty(nueva_red(nueva_comunidad, nuevas_abundancias, red_anterior))
			if (is.array(red_new) && (dim(red_new)[1]) > 2 && (dim(red_new)[2])>2){
			red_final=red_new
			}else{red_final="small"
			}
		
				}else{
		red_final="small" 
		}
	}else{
 red_final="small"}
 return(red_final)
}



setwd("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")
#setwd("C:\\Users\\User\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")

#Load plant species traits...
plant.traits<-read.table("planilla_parametros_FINAL.csv",sep=";", header=T)
plant.traits=plant.traits[1:142,]
#load herbicides traits
herbicides<- read.table("herbicides.csv",sep=",", header=T)
#load plant-pollinator interactions dataset
total<-read.table("inter_plant_pol.csv",sep=";", header=T)

####CREATE THE PLANT POLLINATOR META-WEB
#grouping spp
agrupar<-group_by(total,Plant,Pollinator)#agrupar filas iguales
#agrupar
#countingr
suma_filas<-summarise(agrupar,sum(frequency))
#suma_filas
metared=data.frame(suma_filas)

matriz_meta<-list2matr(metared$Plant,metared$Pollinator,metared$sum.frequency.,matr.type="q")
meta_web1<- as.matrix(as.data.frame(matriz_meta,col.names = names(metared$Pollinator))) 
meta_web<-meta_web1/sum(meta_web1)#INTERACTION MATRIX


escenario="SCE1" #set scenario "SCE1", SCE2, SCE3
simulations<-1000
#lists to save plant community simmulations and estimated network matrics
plant_pol_net_values<-list()
spp_plantas_simul=list()
todas_redes=list()
for (i in 1:simulations){ 

	#1) CHANGES ON PLANT COMMUNITIES 
	#1.1)choose a plant subset
		Nspp= 20 #number of species to be picked to create de initial plant community ##
		prop.p=plant.traits[,"proportion3"]#plant proportions###
		set.seed(NULL)
		rand.plant <- sample(1:142, Nspp, replace=FALSE, prob=prop.p) #sample a subset of plant species #
		#alternative option: all plant species with the same probability
		#rand.plant <- sample(1:142, Nspp, replace=FALSE) 
		
		#create the initial data subset to work with
		set.plants=plant.traits[rand.plant,]
		init.plant.comm<-set.plants$proportion/sum(set.plants$proportion)# Plant species proportions for the initial plant community
		names(init.plant.comm)<-set.plants$plant_species
		init.plant.comm  ##intial (t0) plant community
		
	#1.2 PLANT PARAMETERS FOR LOCAL EXTINCTION OR EXPANSION
		#ALTERNATIVA DE PARÁMETROS
		#weed_risk
		
		life.cyc<-(set.plants$cycle_value5-0.01)/max(set.plants$cycle_value5)
		banco<-(set.plants$seed_value2-0.01)/max(set.plants$seed_value2)
		dispers<-(set.plants$dispersal_value2-0.01)/max(set.plants$dispersal_value2)
		familia<-(set.plants$family_value2-0.01)/max(set.plants$family_value2)
				
		life.cyc_loss<-(set.plants$inv_cycle_value5-0.01)/max(set.plants$inv_cycle_value5)# 
		banco_loss<-(set.plants$inv_seed_value2-0.01)/max(set.plants$inv_seed_value2)
		familia_loss<-(set.plants$inv_family_value2-0.01)/max(set.plants$inv_family_value2) 
			
		weed_risk_perdida<- (life.cyc_loss*banco_loss*familia_loss)/max(life.cyc_loss*banco_loss*familia_loss)
		weed_risk_expansion<-(life.cyc*banco*dispers*familia)/max(life.cyc*banco*dispers*familia)
		
		#herbicide_risk
		biotypes1<-as.matrix(as.data.frame(set.plants[31:38],row.names = set.plants[,1]))
		biotypes=biotypes1/max(biotypes1)
		#herbicides characteristics and scenarios
		her<- as.matrix(as.data.frame(herbicides[,2:9],row.names = herbicides[,1]))
		presion<-(her["select_preasure",])/max(her["select_preasure",])
		dosis<-her[escenario,]/max(her[escenario,])#number of applications of each herbicide (different dosses) in 2 years 
		
		herb_risk= (rowSums(biotypes*presion*dosis))/max(rowSums(biotypes*presion*dosis))
		
		#Management
		###mmodifiers de Moss et al 2019
		if(escenario=="SCE3"){modifiers=1}#*1
		if(escenario=="SCE2"){modifiers=0.67}#0.67
		if(escenario=="SCE1"){modifiers=0.33}#0.33
		
		trait.losing<- ((weed_risk_perdida*(1-herb_risk))/max(weed_risk_perdida*(1-herb_risk)))*modifiers #modifiers va al final para dar valor a #todo el vector, sino se estandariza
		trait.expanding<-((weed_risk_expansion*herb_risk)/max(weed_risk_expansion*herb_risk))*modifiers
		
		## PLANT COMMUNITY MODEL (non-random starting community) ####
		## Time parameters
		years<-15 # how many years you want the model to run for  ###each scenario was designed for two years, so we have 15*2=30 years
		duration<-years##

		#plant simulation
		simulando<-plant_simulation(init.plant.comm,duration,trait.losing,trait.expanding)#,resistance) #here we obtain an array with plant proportions at each time step
		simulando1=rbind(as.vector(init.plant.comm),simulando) ##adding initial plant community to the array 
		spp_plantas_simul[[i]]=simulando1 #to save in the final list
		
		
		#plot curves
		#simu <- data.frame(x = seq_along(simulando1[, 1]),simulando1)
		#Formato long
		#simu <- melt(simu, id.vars = "x")
		#ggplot(simu, aes(x = x, y = value, color = variable)) +
		#geom_line(lwd=1.5)

	#2)NETWORK SIMULATION
	#2.1)INITIAL NETWORK
		#initial plant community
		comunidad_inicial1=simulando1[1,]
				
		#mean number of links by plots in empirical networks to create a random network.. mean 25, median 16
		N_links=25# to be used un network_assemblyV2 function	
		N_pol=15# to be used un network_assemblyV2 function	
		
		#create the initial random network using the initial plant community
		red_inicial=network_assemblyV2(meta_web, nueva_comunidad=comunidad_inicial1, N_pol)#selecting first a subset of pollinators and n_links=sum(web)*100
		#plotweb(red_inicial, arrow="down")
		#red_inicial
		
		#2.2) REASEMBLING NETWORKS AFTER SIMULATIONS
		##reassembling the networks at each cicle of "duration"
		red=list() #strore networks in a list
		for (k in 1:nrow(simulando)){
			comunidad_simulada=as.vector(simulando[k,])#select the corresponding plant community (tn)
			names(comunidad_simulada)=names(simulando[k,])
			
		#Reassembly the network	#using pollinators from previuos temporal network t(n-1).
			if (k==1){red[[k]]=new_reconectar_without_rewiring2(nueva_comunidad=comunidad_simulada, comunidad_anterior=comunidad_inicial1, red_anterior=red_inicial)} 
			if (k>1){
				if(any(red[[k-1]]=="small")){ 
				red[[k]]="small"
				}else{
					red[[k]]=new_reconectar_without_rewiring2(nueva_comunidad=comunidad_simulada, comunidad_anterior=simulando[k-1,],red_anterior=red[[k-1]])}
			}
			}
		#plotweb(red[[15]], arrow="down")
		#red

		
	#2.3)Estimate network metrics
		red_inicial1=list(red_inicial)
		redes=append(red_inicial1,red)###add initial network to the list of networks
		todas_redes[[i]]=redes
		
		plant_pol_net_values_cycle<-array(dim=c((years+1),18)) #storage array ###acá, para 30 años son 16 filas, para 50 años  es 26
		for (n in 1:length(redes)){ 
			
			#if (sum(redes[[n]]!=0)<=1){resumen_final=c(rep(NA,14), paste0("simul",k)) }#esto no va...
			
			if (any(redes[[n]]=="small")){resumen_final=cbind("small", NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, n, i)
			}else{

			#if (dim(redes[[n]])[2]>=2 && dim(redes[[n]])[1]>=2){#exc
				##metrics estimation
				#robustness
				extincion1<-second.extinct(redes[[n]], participant="both", method="abundance", nrep=100) #con both no funciona abundance
				robustez1=robustness(extincion1)
				robustez1
				
				#robustness - only plants
				extincion2<-second.extinct(redes[[n]], participant="lower", method="abundance", nrep=100) #abundance es de menos a mas #abundante
				robustez2=robustness(extincion2)
				robustez2
				
				dimen1<-dim(redes[[n]])
				S_plantas1=dimen1[1]
				S_poli1=dimen1[2]
				links1=sum(redes[[n]])
				evenn1 = networklevel(redes[[n]],index="interaction evenness",intereven="prod")##intereven="prod" or "sum"
				#evenn2 = networklevel(redes[[n]],index="interaction evenness", intereven="prod", effective = TRUE)##with hill numbers Jost 2010
				conn1 = networklevel(redes[[n]],index="connectance")
				nest1 = networklevel(redes[[n]], index="weighted NODF")
				nodf = networklevel(redes[[n]], index="NODF")#unweighted
				fun_com1 = fc(redes[[n]],dist="euclidean", method="average", weighted=TRUE)##### canberra?euclidea
				n_overHL = networklevel(redes[[n]],index="niche overlap")[1]
 				n_overLL = networklevel(redes[[n]],index="niche overlap")[2]
				int_div = networklevel(redes[[n]],index="Shannon diversity")
				gen = networklevel(redes[[n]],index="generality")[1]
				vul = networklevel(redes[[n]],index="generality")[2]
				H2 = H2fun(redes[[n]])[1]

				resumen_final=cbind(S_plantas1, S_poli1,links1, evenn1, conn1,nest1, nodf,robustez1,robustez2,fun_com1,n_overHL, n_overLL, int_div, gen, vul, H2, n, i)
							
							
			}
			plant_pol_net_values_cycle[n,]<-resumen_final
			}
		
		plant_pol_net_values[[i]]<-plant_pol_net_values_cycle
		
	}

#create a final table containing all the results
tabla_final=do.call(rbind, plant_pol_net_values) #deshago la lista
colnames(tabla_final)<- c("S_plantas", "S_poli","links", "evenn", "conn","nest","nodf","robust1","robust", "funct_comp", "N_overHL", "N_overLL", "int_div", "gen", "vul", "H2","Sumulation","cycle")

tabla_final_modifiers=data.frame(tabla_final)
tabla_final_modifiers
#save table
write.table(tabla_final_modifiers, file = "G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\tabla_final_SC1_norew.csv", sep = ",", row.names=F)

#saving the frequency of each plant species in final simulations
spp_plantas_simul	
spp_finales=list()
for (z in 1:simulations){spp_finales[[z]]<-names(spp_plantas_simul[[z]][16,][spp_plantas_simul[[z]][16,]>0])}
frecuencias<-unlist(spp_finales)
suma_frecuencias=data.frame(table(frecuencias))	
suma_frec_ordenado <- suma_frecuencias[order(suma_frecuencias$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_plantas_SC1_norew.csv",sep=";")

#counting networks that became very small
small_net<-which(tabla_final_modifiers[,1]=="small")
length(small_net)

#saving frequency of each plant and pollinator species in final and initial simulations 
#todas_redes
poli_iniciales=list()
for (w in 1:simulations){poli_iniciales[[w]]<-colnames(todas_redes[[w]][[1]])}
frec_poli_inic<-unlist(poli_iniciales)
suma_frec_pol_in=data.frame(table(frec_poli_inic))
suma_frec_ordenado_po_i <- suma_frec_pol_in[order(suma_frec_pol_in$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_po_i,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_poli_iniciales_SC1_norew.csv",sep=";")

poli_finales=list()
for (w in 1:simulations){poli_finales[[w]]<-colnames(todas_redes[[w]][[16]])}
frec_poli_fin<-unlist(poli_finales)
suma_frec_pol_fin=data.frame(table(frec_poli_fin))
suma_frec_ordenado_po_f <- suma_frec_pol_fin[order(suma_frec_pol_fin$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_po_f,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_poli_finales_SC1_norew.csv",sep=";")

plantas_iniciales=list()
for (w in 1:simulations){plantas_iniciales[[w]]<-rownames(todas_redes[[w]][[1]])}
frec_plantas_inic<-unlist(plantas_iniciales)
suma_frec_pl_in=data.frame(table(frec_plantas_inic))
suma_frec_ordenado_pl_i <- suma_frec_pl_in[order(suma_frec_pl_in$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_pl_i,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_pl_iniciales_SC1_norew.csv",sep=";")

plantas_finales=list()
for (w in 1:simulations){plantas_finales[[w]]<-rownames(todas_redes[[w]][[16]])}
frec_plantas_fin<-unlist(plantas_finales)
suma_frec_pl_fin=data.frame(table(frec_plantas_fin))
suma_frec_ordenado_pl_f <- suma_frec_pl_fin[order(suma_frec_pl_fin$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_pl_f,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_pl_finales_SC1_norew.csv",sep=";")

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
# SIN DYNAMIC DE POLLINATORS NI REWIRING
####NO SIRVE- NO SIRVE - NO SIRVE

####3/2/2025  trabajando en este
#function to reassembly the network ###no hay rewiring, la nueva abundancia de polinators establece interacciones en la misma proporcion que en la red anterior

nueva_red<-function(nueva_comunidad, nuevas_abundancias, red_anterior){
	red_anterior_norm=sweep(red_anterior, 2, colSums(red_anterior), "/")
	
	new_mat=matrix(0,length(nueva_comunidad), length(nuevas_abundancias)) #empty matrix with the size of the new community 27/1/25
	rownames(new_mat)=names(nueva_comunidad) 
	colnames(new_mat)=names(nuevas_abundancias)
	
			for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
				sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta
				for (j in 1:length(nuevas_abundancias)){
				sp_poli<- as.character(names(nuevas_abundancias[j]))
					if (sp_planta %in% rownames(red_anterior)){ #because not all plant species are in the meta_web#
					
					if (red_anterior[sp_planta,sp_poli]>0){
						new_mat[sp_planta,sp_poli]= round(red_anterior_norm[sp_planta,sp_poli]*nuevas_abundancias[j],0)
						}else{
						new_mat[sp_planta,sp_poli]= 0
						}
		        }else{
				new_mat[sp_planta,sp_poli]<-0
				}
			}
			}
		red_final=new_mat
		return(red_final)
}



new_reconectar_without_dynam_nochange<-function(nueva_comunidad, comunidad_anterior, red_anterior){
	if (is.array(red_anterior) && (dim(red_anterior)[1]) > 2 && (dim(red_anterior)[2])>2 && sum(nueva_comunidad>0)>2 && sum(comunidad_anterior>0)>2){
	#estimates new pollinators frequency 
		nuevas_abundancias = colSums(red_anterior)
		if (length(nuevas_abundancias)>2){
		red_new=empty(nueva_red(nueva_comunidad, nuevas_abundancias, red_anterior))
			if (is.array(red_new) && (dim(red_new)[1]) > 2 && (dim(red_new)[2])>2){
			red_final=red_new
			}else{red_final="small"
			}
		
				}else{
		red_final="small"
		}
	}else{
 red_final="small"}
 return(red_final)
}



setwd("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")
#setwd("C:\\Users\\User\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")

#Load plant species traits...
plant.traits<-read.table("planilla_parametros_FINAL.csv",sep=";", header=T)
plant.traits=plant.traits[1:142,]
#load herbicides traits
herbicides<- read.table("herbicides.csv",sep=",", header=T)
#load plant-pollinator interactions dataset
total<-read.table("inter_plant_pol.csv",sep=";", header=T)

####CREATE THE PLANT POLLINATOR META-WEB
#grouping spp
agrupar<-group_by(total,Plant,Pollinator)#agrupar filas iguales
#agrupar
#countingr
suma_filas<-summarise(agrupar,sum(frequency))
#suma_filas
metared=data.frame(suma_filas)

matriz_meta<-list2matr(metared$Plant,metared$Pollinator,metared$sum.frequency.,matr.type="q")
meta_web1<- as.matrix(as.data.frame(matriz_meta,col.names = names(metared$Pollinator))) 
meta_web<-meta_web1/sum(meta_web1)#INTERACTION MATRIX


escenario="SCE1" #set scenario "SCE1", SCE2, SCE3
simulations<-1000
#lists to save plant community simmulations and estimated network matrics
plant_pol_net_values<-list()
spp_plantas_simul=list()
todas_redes=list()
for (i in 1:simulations){ 

	#1) CHANGES ON PLANT COMMUNITIES 
	#1.1)choose a plant subset
		Nspp= 20 #number of species to be picked to create de initial plant community ##
		prop.p=plant.traits[,"proportion3"]#plant proportions###
		set.seed(NULL)
		rand.plant <- sample(1:142, Nspp, replace=FALSE, prob=prop.p) #sample a subset of plant species #
		#alternative option: all plant species with the same probability
		#rand.plant <- sample(1:142, Nspp, replace=FALSE) 
		
		#create the initial data subset to work with
		set.plants=plant.traits[rand.plant,]
		init.plant.comm<-set.plants$proportion/sum(set.plants$proportion)# Plant species proportions for the initial plant community
		names(init.plant.comm)<-set.plants$plant_species
		init.plant.comm  ##intial (t0) plant community
		
	#1.2 PLANT PARAMETERS FOR LOCAL EXTINCTION OR EXPANSION
		#ALTERNATIVA DE PARÁMETROS
		#weed_risk
		
		life.cyc<-(set.plants$cycle_value5-0.01)/max(set.plants$cycle_value5)
		banco<-(set.plants$seed_value2-0.01)/max(set.plants$seed_value2)
		dispers<-(set.plants$dispersal_value2-0.01)/max(set.plants$dispersal_value2)
		familia<-(set.plants$family_value2-0.01)/max(set.plants$family_value2)
				
		life.cyc_loss<-(set.plants$inv_cycle_value5-0.01)/max(set.plants$inv_cycle_value5)# 
		banco_loss<-(set.plants$inv_seed_value2-0.01)/max(set.plants$inv_seed_value2)
		familia_loss<-(set.plants$inv_family_value2-0.01)/max(set.plants$inv_family_value2) 
			
		weed_risk_perdida<- (life.cyc_loss*banco_loss*familia_loss)/max(life.cyc_loss*banco_loss*familia_loss)
		weed_risk_expansion<-(life.cyc*banco*dispers*familia)/max(life.cyc*banco*dispers*familia)
		
		#herbicide_risk
		biotypes1<-as.matrix(as.data.frame(set.plants[31:38],row.names = set.plants[,1]))
		biotypes=biotypes1/max(biotypes1)
		#herbicides characteristics and scenarios
		her<- as.matrix(as.data.frame(herbicides[,2:9],row.names = herbicides[,1]))
		presion<-(her["select_preasure",])/max(her["select_preasure",])
		dosis<-her[escenario,]/max(her[escenario,])#number of applications of each herbicide (different dosses) in 2 years 
		
		herb_risk= (rowSums(biotypes*presion*dosis))/max(rowSums(biotypes*presion*dosis))
		
		#Management
		###mmodifiers de Moss et al 2019
		if(escenario=="SCE3"){modifiers=1}#*1
		if(escenario=="SCE2"){modifiers=0.67}#0.67
		if(escenario=="SCE1"){modifiers=0.33}#0.33
		
		trait.losing<- ((weed_risk_perdida*(1-herb_risk))/max(weed_risk_perdida*(1-herb_risk)))*modifiers #modifiers va al final para dar valor a #todo el vector, sino se estandariza
		trait.expanding<-((weed_risk_expansion*herb_risk)/max(weed_risk_expansion*herb_risk))*modifiers
		
		## PLANT COMMUNITY MODEL (non-random starting community) ####
		## Time parameters
		years<-15 # how many years you want the model to run for  ###each scenario was designed for two years, so we have 15*2=30 years
		duration<-years##

		#plant simulation
		simulando<-plant_simulation(init.plant.comm,duration,trait.losing,trait.expanding)#,resistance) #here we obtain an array with plant proportions at each time step
		simulando1=rbind(as.vector(init.plant.comm),simulando) ##adding initial plant community to the array 
		spp_plantas_simul[[i]]=simulando1 #to save in the final list
		
		
		#plot curves
		#simu <- data.frame(x = seq_along(simulando1[, 1]),simulando1)
		#Formato long
		#simu <- melt(simu, id.vars = "x")
		#ggplot(simu, aes(x = x, y = value, color = variable)) +
		#geom_line(lwd=1.5)

	#2)NETWORK SIMULATION
	#2.1)INITIAL NETWORK
		#initial plant community
		comunidad_inicial1=simulando1[1,]
				
		#mean number of links by plots in empirical networks to create a random network.. mean 25, median 16
		N_links=25# to be used un network_assemblyV2 function	
		N_pol=15# to be used un network_assemblyV2 function	
		
		#create the initial random network using the initial plant community
		red_inicial=network_assemblyV2(meta_web, nueva_comunidad=comunidad_inicial1, N_pol)#selecting first a subset of pollinators and n_links=sum(web)*100
		#plotweb(red_inicial, arrow="down")
		#red_inicial
		
		#2.2) REASEMBLING NETWORKS AFTER SIMULATIONS
		##reassembling the networks at each cicle of "duration"
		red=list() #strore networks in a list
		for (k in 1:nrow(simulando)){
			comunidad_simulada=as.vector(simulando[k,])#select the corresponding plant community (tn)
			names(comunidad_simulada)=names(simulando[k,])
			
		#Reassembly the network	#using pollinators from previuos temporal network t(n-1).
			if (k==1){red[[k]]=new_reconectar_without_dynam_nochange(nueva_comunidad=comunidad_simulada, comunidad_anterior=comunidad_inicial1, red_anterior=red_inicial)} 
			if (k>1){
				if(any(red[[k-1]]=="small")){ 
				red[[k]]="small"
				}else{
					red[[k]]=new_reconectar_without_dynam_nochange(nueva_comunidad=comunidad_simulada, comunidad_anterior=simulando[k-1,],red_anterior=red[[k-1]])}
			}
			}
		#plotweb(red[[15]], arrow="down")
		#red

		
	#2.3)Estimate network metrics
		red_inicial1=list(red_inicial)
		redes=append(red_inicial1,red)###add initial network to the list of networks
		todas_redes[[i]]=redes
		
		plant_pol_net_values_cycle<-array(dim=c((years+1),18)) #storage array ###acá, para 30 años son 16 filas, para 50 años  es 26
		for (n in 1:length(redes)){ 
			
			#if (sum(redes[[n]]!=0)<=1){resumen_final=c(rep(NA,14), paste0("simul",k)) }#esto no va...
			
			if (any(redes[[n]]=="small")){resumen_final=cbind("small", NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, n, i)
			}else{

			#if (dim(redes[[n]])[2]>=2 && dim(redes[[n]])[1]>=2){#exc
				##metrics estimation
				#robustness
				extincion1<-second.extinct(redes[[n]], participant="both", method="abundance", nrep=100) #con both no funciona abundance
				robustez1=robustness(extincion1)
				robustez1
				
				#robustness - only plants
				extincion2<-second.extinct(redes[[n]], participant="lower", method="abundance", nrep=100) #abundance es de menos a mas #abundante
				robustez2=robustness(extincion2)
				robustez2
				
				dimen1<-dim(redes[[n]])
				S_plantas1=dimen1[1]
				S_poli1=dimen1[2]
				links1=sum(redes[[n]])
				evenn1 = networklevel(redes[[n]],index="interaction evenness",intereven="prod")##intereven="prod" or "sum"
				#evenn2 = networklevel(redes[[n]],index="interaction evenness", intereven="prod", effective = TRUE)##with hill numbers Jost 2010
				conn1 = networklevel(redes[[n]],index="connectance")
				nest1 = networklevel(redes[[n]], index="weighted NODF")
				nodf = networklevel(redes[[n]], index="NODF")#unweighted
				fun_com1 = fc(redes[[n]],dist="euclidean", method="average", weighted=TRUE)##### canberra?euclidea
				n_overHL = networklevel(redes[[n]],index="niche overlap")[1]
 				n_overLL = networklevel(redes[[n]],index="niche overlap")[2]
				int_div = networklevel(redes[[n]],index="Shannon diversity")
				gen = networklevel(redes[[n]],index="generality")[1]
				vul = networklevel(redes[[n]],index="generality")[2]
				H2 = H2fun(redes[[n]])[1]

				resumen_final=cbind(S_plantas1, S_poli1,links1, evenn1, conn1,nest1, nodf,robustez1,robustez2,fun_com1,n_overHL, n_overLL, int_div, gen, vul, H2, n, i)
							
							
			}
			plant_pol_net_values_cycle[n,]<-resumen_final
			}
		
		plant_pol_net_values[[i]]<-plant_pol_net_values_cycle
		
	}

#create a final table containing all the results
tabla_final=do.call(rbind, plant_pol_net_values) #deshago la lista
colnames(tabla_final)<- c("S_plantas", "S_poli","links", "evenn", "conn","nest", "nodf","robust1","robust2", "funct_comp", "N_overHL", "N_overLL", "int_div", "gen", "vul", "H2","Sumulation","cycle")

tabla_final_modifiers=data.frame(tabla_final)
tabla_final_modifiers
#save table
write.table(tabla_final_modifiers, file = "G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\tabla_final_SC1_nodyn.csv", sep = ",", row.names=F)

#saving the frequency of each plant species in final simulations
spp_plantas_simul	
spp_finales=list()
for (z in 1:simulations){spp_finales[[z]]<-names(spp_plantas_simul[[z]][16,][spp_plantas_simul[[z]][16,]>0])}
frecuencias<-unlist(spp_finales)
suma_frecuencias=data.frame(table(frecuencias))	
suma_frec_ordenado <- suma_frecuencias[order(suma_frecuencias$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_plantas_SC1_nodyn.csv",sep=";")

#counting networks that became very small
small_net<-which(tabla_final_modifiers[,1]=="small")
length(small_net)

#saving frequency of each plant and pollinator species in final and initial simulations 
#todas_redes
poli_iniciales=list()
for (w in 1:simulations){poli_iniciales[[w]]<-colnames(todas_redes[[w]][[1]])}
frec_poli_inic<-unlist(poli_iniciales)
suma_frec_pol_in=data.frame(table(frec_poli_inic))
suma_frec_ordenado_po_i <- suma_frec_pol_in[order(suma_frec_pol_in$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_po_i,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_poli_iniciales_SC1_nodyn.csv",sep=";")

poli_finales=list()
for (w in 1:simulations){poli_finales[[w]]<-colnames(todas_redes[[w]][[16]])}
frec_poli_fin<-unlist(poli_finales)
suma_frec_pol_fin=data.frame(table(frec_poli_fin))
suma_frec_ordenado_po_f <- suma_frec_pol_fin[order(suma_frec_pol_fin$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_po_f,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_poli_finales_SC1_nodyn.csv",sep=";")

plantas_iniciales=list()
for (w in 1:simulations){plantas_iniciales[[w]]<-rownames(todas_redes[[w]][[1]])}
frec_plantas_inic<-unlist(plantas_iniciales)
suma_frec_pl_in=data.frame(table(frec_plantas_inic))
suma_frec_ordenado_pl_i <- suma_frec_pl_in[order(suma_frec_pl_in$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_pl_i,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_pl_iniciales_SC1_nodyn.csv",sep=";")

plantas_finales=list()
for (w in 1:simulations){plantas_finales[[w]]<-rownames(todas_redes[[w]][[16]])}
frec_plantas_fin<-unlist(plantas_finales)
suma_frec_pl_fin=data.frame(table(frec_plantas_fin))
suma_frec_ordenado_pl_f <- suma_frec_pl_fin[order(suma_frec_pl_fin$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_pl_f,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_pl_finales_SC1_nodyn.csv",sep=";")



##############################################################################################################################################################################################################################
#ahora la abundancia de poli no cambia pero pueden ocupar cualquier planta con la que interactue en la meta_web

#function to reassembly the network ###no hay rewiring, la nueva abundancia de polinators establece interacciones en la misma proporcion que en la red anterior

nueva_red<-function(nueva_comunidad, nuevas_abundancias, red_anterior){
	
	#esta debería estar armada con los valores de la meta_web
	
	sub_meta_web=matrix(0,length(nueva_comunidad), length(nuevas_abundancias))#storage matrix
	rownames(sub_meta_web)=names(nueva_comunidad) 
	colnames(sub_meta_web)=names (nuevas_abundancias)
	
	for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
		sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta
		if (sp_planta %in% rownames(meta_web)){ #because not all plant species are in the meta_web
			for (j in 1:length(nuevas_abundancias)){
				sp_poli<- names(nuevas_abundancias)[j]
				
				###ESTO TIENE Q SER A 0-1 Y LUEGO MULTIP * NUEVA COMUNIAD DE PLANTAS Y LUEGO NORMALIZO.
						
				sub_meta_web[i,j]=meta_web[sp_planta,sp_poli] #select the value of the interaction in meta_web to copy in sub_meta_web
				}}}
				
	 sub_meta_web_bin<- ifelse(sub_meta_web > 0, 1, 0)
	 prob_por_planta= sub_meta_web_bin * nueva_comunidad #con esto doy prob a cada planta posible de la red de ser visitada según su proporcion en el plot
	 sub_meta_norm=sweep(prob_por_planta, 2, colSums(prob_por_planta), "/")	
	
	new_mat=matrix(0,length(nueva_comunidad), length(nuevas_abundancias)) #empty matrix with the size of the new community 27/1/25
	rownames(new_mat)=names(nueva_comunidad) 
	colnames(new_mat)=names(nuevas_abundancias)
	
			for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
				sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta
				for (j in 1:length(nuevas_abundancias)){
				sp_poli<- as.character(names(nuevas_abundancias[j]))
					if (sp_planta %in% rownames(red_anterior)){ #because not all plant species are in the meta_web#
					
					if (red_anterior[sp_planta,sp_poli]>0){
						new_mat[sp_planta,sp_poli]= round(sub_meta_norm[sp_planta,sp_poli]*nuevas_abundancias[j],0)
						}else{
						new_mat[sp_planta,sp_poli]= 0
						}
		        }else{
				new_mat[sp_planta,sp_poli]<-0
				}
			}
			}
		red_final=new_mat
		return(red_final)
}



new_reconectar_without_dynam_meta<-function(nueva_comunidad, comunidad_anterior, red_anterior){
	if (is.array(red_anterior) && (dim(red_anterior)[1]) > 2 && (dim(red_anterior)[2])>2 && sum(nueva_comunidad>0)>2 && sum(comunidad_anterior>0)>2){
	#estimates new pollinators frequency 
		nuevas_abundancias = colSums(red_anterior)
		if (length(nuevas_abundancias)>2){
		red_new=empty(nueva_red(nueva_comunidad, nuevas_abundancias, red_anterior))
			if (is.array(red_new) && (dim(red_new)[1]) > 2 && (dim(red_new)[2])>2){
			red_final=red_new
			}else{red_final="small"
			}
		
				}else{
		red_final="small"
		}
	}else{
 red_final="small"}
 return(red_final)
}



setwd("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")
#setwd("C:\\Users\\User\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")

#Load plant species traits...
plant.traits<-read.table("planilla_parametros_FINAL.csv",sep=";", header=T)
plant.traits=plant.traits[1:142,]
#load herbicides traits
herbicides<- read.table("herbicides.csv",sep=",", header=T)
#load plant-pollinator interactions dataset
total<-read.table("inter_plant_pol.csv",sep=";", header=T)

####CREATE THE PLANT POLLINATOR META-WEB
#grouping spp
agrupar<-group_by(total,Plant,Pollinator)#agrupar filas iguales
#agrupar
#countingr
suma_filas<-summarise(agrupar,sum(frequency))
#suma_filas
metared=data.frame(suma_filas)

matriz_meta<-list2matr(metared$Plant,metared$Pollinator,metared$sum.frequency.,matr.type="q")
meta_web1<- as.matrix(as.data.frame(matriz_meta,col.names = names(metared$Pollinator))) 
meta_web<-meta_web1/sum(meta_web1)#INTERACTION MATRIX


escenario="SCE3" #set scenario "SCE1", SCE2, SCE3
simulations<-1000
#lists to save plant community simmulations and estimated network matrics
plant_pol_net_values<-list()
spp_plantas_simul=list()
todas_redes=list()
for (i in 1:simulations){ 

	#1) CHANGES ON PLANT COMMUNITIES 
	#1.1)choose a plant subset
		Nspp= 20 #number of species to be picked to create de initial plant community ##
		prop.p=plant.traits[,"proportion3"]#plant proportions###
		set.seed(NULL)
		rand.plant <- sample(1:142, Nspp, replace=FALSE, prob=prop.p) #sample a subset of plant species #
		#alternative option: all plant species with the same probability
		#rand.plant <- sample(1:142, Nspp, replace=FALSE) 
		
		#create the initial data subset to work with
		set.plants=plant.traits[rand.plant,]
		init.plant.comm<-set.plants$proportion/sum(set.plants$proportion)# Plant species proportions for the initial plant community
		names(init.plant.comm)<-set.plants$plant_species
		init.plant.comm  ##intial (t0) plant community
		
	#1.2 PLANT PARAMETERS FOR LOCAL EXTINCTION OR EXPANSION
		#ALTERNATIVA DE PARÁMETROS
		#weed_risk
		
		life.cyc<-(set.plants$cycle_value5-0.01)/max(set.plants$cycle_value5)
		banco<-(set.plants$seed_value2-0.01)/max(set.plants$seed_value2)
		dispers<-(set.plants$dispersal_value2-0.01)/max(set.plants$dispersal_value2)
		familia<-(set.plants$family_value2-0.01)/max(set.plants$family_value2)
				
		life.cyc_loss<-(set.plants$inv_cycle_value5-0.01)/max(set.plants$inv_cycle_value5)# 
		banco_loss<-(set.plants$inv_seed_value2-0.01)/max(set.plants$inv_seed_value2)
		familia_loss<-(set.plants$inv_family_value2-0.01)/max(set.plants$inv_family_value2) 
			
		weed_risk_perdida<- (life.cyc_loss*banco_loss*familia_loss)/max(life.cyc_loss*banco_loss*familia_loss)
		weed_risk_expansion<-(life.cyc*banco*dispers*familia)/max(life.cyc*banco*dispers*familia)
		
		#herbicide_risk
		biotypes1<-as.matrix(as.data.frame(set.plants[31:38],row.names = set.plants[,1]))
		biotypes=biotypes1/max(biotypes1)
		#herbicides characteristics and scenarios
		her<- as.matrix(as.data.frame(herbicides[,2:9],row.names = herbicides[,1]))
		presion<-(her["select_preasure",])/max(her["select_preasure",])
		dosis<-her[escenario,]/max(her[escenario,])#number of applications of each herbicide (different dosses) in 2 years 
		
		herb_risk= (rowSums(biotypes*presion*dosis))/max(rowSums(biotypes*presion*dosis))
		
		#Management
		###mmodifiers de Moss et al 2019
		if(escenario=="SCE3"){modifiers=1}#*1
		if(escenario=="SCE2"){modifiers=0.67}#0.67
		if(escenario=="SCE1"){modifiers=0.33}#0.33
		
		trait.losing<- ((weed_risk_perdida*(1-herb_risk))/max(weed_risk_perdida*(1-herb_risk)))*modifiers #modifiers va al final para dar valor a #todo el vector, sino se estandariza
		trait.expanding<-((weed_risk_expansion*herb_risk)/max(weed_risk_expansion*herb_risk))*modifiers
		
		## PLANT COMMUNITY MODEL (non-random starting community) ####
		## Time parameters
		years<-15 # how many years you want the model to run for  ###each scenario was designed for two years, so we have 15*2=30 years
		duration<-years##

		#plant simulation
		simulando<-plant_simulation(init.plant.comm,duration,trait.losing,trait.expanding)#,resistance) #here we obtain an array with plant proportions at each time step
		simulando1=rbind(as.vector(init.plant.comm),simulando) ##adding initial plant community to the array 
		spp_plantas_simul[[i]]=simulando1 #to save in the final list
		
		
		#plot curves
		#simu <- data.frame(x = seq_along(simulando1[, 1]),simulando1)
		#Formato long
		#simu <- melt(simu, id.vars = "x")
		#ggplot(simu, aes(x = x, y = value, color = variable)) +
		#geom_line(lwd=1.5)

	#2)NETWORK SIMULATION
	#2.1)INITIAL NETWORK
		#initial plant community
		comunidad_inicial1=simulando1[1,]
				
		#mean number of links by plots in empirical networks to create a random network.. mean 25, median 16
		N_links=25# to be used un network_assemblyV2 function	
		N_pol=15# to be used un network_assemblyV2 function	
		
		#create the initial random network using the initial plant community
		red_inicial=network_assemblyV2(meta_web, nueva_comunidad=comunidad_inicial1, N_pol)#selecting first a subset of pollinators and n_links=sum(web)*100
		#plotweb(red_inicial, arrow="down")
		#red_inicial
		
		#2.2) REASEMBLING NETWORKS AFTER SIMULATIONS
		##reassembling the networks at each cicle of "duration"
		red=list() #strore networks in a list
		for (k in 1:nrow(simulando)){
			comunidad_simulada=as.vector(simulando[k,])#select the corresponding plant community (tn)
			names(comunidad_simulada)=names(simulando[k,])
			
		#Reassembly the network	#using pollinators from previuos temporal network t(n-1).
			if (k==1){red[[k]]=new_reconectar_without_dynam_meta(nueva_comunidad=comunidad_simulada, comunidad_anterior=comunidad_inicial1, red_anterior=red_inicial)} 
			if (k>1){
				if(any(red[[k-1]]=="small")){ 
				red[[k]]="small"
				}else{
					red[[k]]=new_reconectar_without_dynam_meta(nueva_comunidad=comunidad_simulada, comunidad_anterior=simulando[k-1,],red_anterior=red[[k-1]])}
			}
			}
		#plotweb(red[[15]], arrow="down")
		#red

		
	#2.3)Estimate network metrics
		red_inicial1=list(red_inicial)
		redes=append(red_inicial1,red)###add initial network to the list of networks
		todas_redes[[i]]=redes
		
			plant_pol_net_values_cycle<-array(dim=c((years+1),18)) #storage array ###acá, para 30 años son 16 filas, para 50 años  es 26
		for (n in 1:length(redes)){ 
			
			#if (sum(redes[[n]]!=0)<=1){resumen_final=c(rep(NA,14), paste0("simul",k)) }#esto no va...
			
			if (any(redes[[n]]=="small")){resumen_final=cbind("small", NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, n, i)
			}else{

			#if (dim(redes[[n]])[2]>=2 && dim(redes[[n]])[1]>=2){#exc
				##metrics estimation
				#robustness
				extincion1<-second.extinct(redes[[n]], participant="both", method="abundance", nrep=100) #con both no funciona abundance
				robustez1=robustness(extincion1)
				robustez1
				
				#robustness - only plants
				extincion2<-second.extinct(redes[[n]], participant="lower", method="abundance", nrep=100) #abundance es de menos a mas #abundante
				robustez2=robustness(extincion2)
				robustez2
				
				dimen1<-dim(redes[[n]])
				S_plantas1=dimen1[1]
				S_poli1=dimen1[2]
				links1=sum(redes[[n]])
				evenn1 = networklevel(redes[[n]],index="interaction evenness",intereven="prod")##intereven="prod" or "sum"
				#evenn2 = networklevel(redes[[n]],index="interaction evenness", intereven="prod", effective = TRUE)##with hill numbers Jost 2010
				conn1 = networklevel(redes[[n]],index="connectance")
				nest1 = networklevel(redes[[n]], index="weighted NODF")
				nodf = networklevel(redes[[n]], index="NODF")#unweighted
				fun_com1 = fc(redes[[n]],dist="euclidean", method="average", weighted=TRUE)##### canberra?euclidea
				n_overHL = networklevel(redes[[n]],index="niche overlap")[1]
 				n_overLL = networklevel(redes[[n]],index="niche overlap")[2]
				int_div = networklevel(redes[[n]],index="Shannon diversity")
				gen = networklevel(redes[[n]],index="generality")[1]
				vul = networklevel(redes[[n]],index="generality")[2]
				H2 = H2fun(redes[[n]])[1]

				resumen_final=cbind(S_plantas1, S_poli1,links1, evenn1, conn1,nest1, nodf,robustez1,robustez2,fun_com1,n_overHL, n_overLL, int_div, gen, vul, H2, n, i)
							
			}
			plant_pol_net_values_cycle[n,]<-resumen_final
			}
		plant_pol_net_values_cycle
		plant_pol_net_values[[i]]<-plant_pol_net_values_cycle
		
	}

#create a final table containing all the results
tabla_final=do.call(rbind, plant_pol_net_values) #deshago la lista
colnames(tabla_final)<- c("S_plantas", "S_poli","links", "evenn", "conn","nest", "nodf","robust1","robust2", "funct_comp", "N_overHL", "N_overLL", "int_div", "gen", "vul", "H2","Sumulation","cycle")

tabla_final_modifiers=data.frame(tabla_final)
tabla_final_modifiers
#save table
write.table(tabla_final_modifiers, file = "G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\tabla_final_SC3_nodyn_meta.csv", sep = ",", row.names=F)

#saving the frequency of each plant species in final simulations
spp_plantas_simul	
spp_finales=list()
for (z in 1:simulations){spp_finales[[z]]<-names(spp_plantas_simul[[z]][16,][spp_plantas_simul[[z]][16,]>0])}
frecuencias<-unlist(spp_finales)
suma_frecuencias=data.frame(table(frecuencias))	
suma_frec_ordenado <- suma_frecuencias[order(suma_frecuencias$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_plantas_SC3_nodyn_meta.csv",sep=";")

#counting networks that became very small
small_net<-which(tabla_final_modifiers[,1]=="small")
length(small_net)

#saving frequency of each plant and pollinator species in final and initial simulations 
#todas_redes
poli_iniciales=list()
for (w in 1:simulations){poli_iniciales[[w]]<-colnames(todas_redes[[w]][[1]])}
frec_poli_inic<-unlist(poli_iniciales)
suma_frec_pol_in=data.frame(table(frec_poli_inic))
suma_frec_ordenado_po_i <- suma_frec_pol_in[order(suma_frec_pol_in$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_po_i,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_poli_iniciales_SC3_nodyn_meta.csv",sep=";")

poli_finales=list()
for (w in 1:simulations){poli_finales[[w]]<-colnames(todas_redes[[w]][[16]])}
frec_poli_fin<-unlist(poli_finales)
suma_frec_pol_fin=data.frame(table(frec_poli_fin))
suma_frec_ordenado_po_f <- suma_frec_pol_fin[order(suma_frec_pol_fin$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_po_f,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_poli_finales_SC3_nodyn_meta.csv",sep=";")

plantas_iniciales=list()
for (w in 1:simulations){plantas_iniciales[[w]]<-rownames(todas_redes[[w]][[1]])}
frec_plantas_inic<-unlist(plantas_iniciales)
suma_frec_pl_in=data.frame(table(frec_plantas_inic))
suma_frec_ordenado_pl_i <- suma_frec_pl_in[order(suma_frec_pl_in$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_pl_i,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_pl_iniciales_SC3_nodyn_meta.csv",sep=";")

plantas_finales=list()
for (w in 1:simulations){plantas_finales[[w]]<-rownames(todas_redes[[w]][[16]])}
frec_plantas_fin<-unlist(plantas_finales)
suma_frec_pl_fin=data.frame(table(frec_plantas_fin))
suma_frec_ordenado_pl_f <- suma_frec_pl_fin[order(suma_frec_pl_fin$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_pl_f,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_pl_finales_SC3_nodyn_meta.csv",sep=";")

############################################################################################################################################################################################################################################################################################################################################################################################################################################################
###12/2/2025
###SIN DINÁMICA (NO CAMBIO ABUND DE POLIS PERO EL MISMO REWIRING QUE EL ORIGINAL

library(RColorBrewer)
library(tidyverse) 
library(sads) 
library(bipartite)
library(igraph)
library(gplots)
library(vegan)
library(reshape)

###First, some functions...
##list to matrix function
list2matr<-function(sp1,sp2,effects,matr.type="q"){
  nodelist1u<-array(unique(sp1))
  nodelist2u<-array(unique(sp2))
  m<-length(nodelist1u)
  n<-length(nodelist2u)
  matr<-matrix(0,m,n) #Matrix
  colnames(matr)<-nodelist2u   #######
  row.names(matr)<-nodelist1u
  for (i in 1:m){
    iu<-nodelist1u[i]
    ii<-which(sp1==iu)
    nodelist2ii<-sp2[ii]
    link.weightii<-effects[ii]
    for (j in 1:length(ii)){
      jj<-which(nodelist2u==nodelist2ii[j])
      if (matr.type=="b") matr[i,jj]<-1 #Binary matrix
      if (matr.type=="q") matr[i,jj]<-link.weightii[j] #Weighted matrix
    }
  }
  res<-list(matr)
  names(res)<-c("matr")
  return(res)
}


################NUEVA PLANT SIMULATION	11/10/2023
####Function for Plant community changes simulation##
plant_simulation<-function(init.plant.comm,duration,trait.losing,trait.expanding){#,resistance){ 
	#'''simulate losses and expansion of plant species
	#	input:a vector of an initial plant comunity, duration=number of cicles, vectors of probabilities of loss, expansion and resistance preasure
	#	output:an array containing the community composition at each time step (duration)
	#	'''
	
	## Storage arrays ####
	plant.comm_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for the data generated from the model
	plant.loss_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for coverage loss data generated from the model
	plant.exp_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for coverage expansion data generated from the model

	for (j in 1:duration){ # loop through each time step
								
		# Select the starting plant community
		if(j<=1){plant.comm=init.plant.comm} # set the plant community as the starting coverage if the time step is 1
		if(j>1){plant.comm=plant.comm_array[j-1,]} # if not, use the plant community from the previous time step
		
		##ESTE NUNCA TIENE QUE SER NAN		
		# calculate the loss of plant coverage due to herbicide application and resistance  #
		plant.loss_array[j,]<-plant.comm-(plant.comm*trait.losing)#lo que queda de la comunidad anterior
		
	    #set a threshold, a sensible proportion at wich we would not expect a plant population that covers <0.01 of the transect will be viable
		plant.loss_array[j,][plant.loss_array[j,]<0.001]<-0 ##fernando recomends 0.1 o 0.5. #!!! We have to set a real threshold here!
		#if 1 is the space we have, 0.01 is the 1% of area occupied by the plant. estoy usando 0.001
		
		space<-1-sum(plant.loss_array[j,]) # calculate the space available due to the loss of plants
			
		# calculate the proportional ability of plant species to expand into the space 
		# before that, I eliminate plant species that we lose at the previous step
		surviving.plants=plant.loss_array[j,] 
		
		surviving.plants[surviving.plants>0]<-1 #binary vector to "delete" lost species from the expanding vector 
		expanding1=trait.expanding * surviving.plants ## parameter of expansion 
		
		if (sum(expanding1)>0){plant.exp.ability<-(expanding1)/sum(expanding1)}
		if (sum(expanding1)==0){plant.exp.ability<- expanding1} #to avoid NaNs
		
		#Calculate the expansion of plants (more resistant plants more expantion chances)	
		plant.exp_array[j,]<-space*plant.exp.ability ##
		
		# calculate the change in plant community coverage after loss and expansion
		if (sum(plant.exp_array[j,])>0){plant.comm_array[j,]<-(plant.loss_array[j,]+plant.exp_array[j,])/sum(plant.loss_array[j,]+plant.exp_array[j,])}
		if (sum(plant.exp_array[j,])==0){plant.comm_array[j,]<-(plant.loss_array[j,]+plant.exp_array[j,])}
	}
	
	colnames(plant.comm_array)<-names(init.plant.comm)
	return(plant.comm_array)  #### return an array with simulations equal to duration 
}


######NETWORK ASSEMBLY function to construct the initial network###
network_assemblyV2<-function(meta_web, nueva_comunidad, N_pol){ 
	#"""Function to generate networks using an initial plant community, and 
	#based on the links with pollinators from the meta-web
	#I used mgen function (in bipartite) to genrate a random network using a probability matrix (based on pollinator frequency and plant species proportions)
	#input: meta_web, nueva_comunidad (initial plant community or plant community at each step of the simulation), N_links (mean number of links  
	#estimated from empirical networks by plot)
	#output: a plant-pollinator network based on the probability matrix and preserving the proportion of links 
	#""" 

	#select, previously, a subset of pollinators of size equal to the mean number of polli in empirical networks
	#give probabilities according to species frequency in meta_Web
	set.seed(NULL)
	polis=sample(1:ncol(meta_web),N_pol,prob=colSums(meta_web),replace=FALSE)
	nomb_polis=colnames(meta_web)[polis]
	
	#create a matrix containing plants from nueva_comunidad interacting with a the subset of pollinators, with values of interactions available in th meta-web
	sub_meta_web=matrix(0,length(nueva_comunidad), N_pol)#storage matrix
	rownames(sub_meta_web)=names(nueva_comunidad) 
	colnames(sub_meta_web)=nomb_polis
	
	for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
		sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta
		if (sp_planta %in% rownames(meta_web)){ #because not all plant species are in the meta_web
			for (j in 1:N_pol){
				sp_poli<- nomb_polis[j]
				sub_meta_web[i,j]=meta_web[sp_planta,sp_poli] #select the value of the interaction in meta_web to copy in sub_meta_web
			}
		}
	}
	
	#define pollinators relative frequency 
	suma_col=colSums(as.matrix(sub_meta_web))###
	pol.frec=suma_col/sum(suma_col)
		
	#define proportion of entomophilous plant species (deleting plants that not interact with pollinators from nueva_comunidad)
	nueva_comunidad2=nueva_comunidad[nueva_comunidad>0] 
	plant.prop=nueva_comunidad2/sum(nueva_comunidad2)

	#interaction probabilities matrix 
	sub_meta_web1=sub_meta_web
	sub_meta_web1[sub_meta_web1>=1]=1##forbidden links matrix.presence-absence links to give probabily 0 to interactions that do not exist in meta_web
	pmat=((plant.prop%*%t(pol.frec))*sub_meta_web1)/sum((plant.prop%*%t(pol.frec))*sub_meta_web1)#construct the matrix of probabilities, combining vectors and forbidden link mat
	rownames(pmat)=names(nueva_comunidad2)
	pmat= empty(pmat)#delete empty rows and columns if necesary
	
	#new random network based on pmat
	#number of links equal to sub_meta_web (proportion of links), all species remain in the network, quantitative network
	#Setting keep.species to FALSE may (but need not) cause a loss of species and networks have a small number of species. With keep.species=TRUE 
	#all selected species remain in the network. Results are quite different.
		set.seed(NULL)
		L=sum(sub_meta_web)*100
		nueva_red=mgen(pmat,n=L,keep.species=TRUE, rep.cell=TRUE, autotransform="sum")### Alternative:create the network usin "sample" function (is quite simmilar)
		rownames(nueva_red)=rownames(pmat)
		colnames(nueva_red)=colnames(pmat)
	
	return(nueva_red)
}

##function for interaction rewiring probabilities
rew.rule <- function(M){
	#function to build a forbidden links matrix based on similarity of interacctions between species (from Maia et al 2021)
	#Input: MAtrix of interactions
	#output: Matrix built based on similarity between species interactions
  rede <- t(M) #aca transpuse , tambien deberia transponer la salida
  rede <- as.matrix(rede)
  
  PP <- 1-as.matrix(vegdist(rede, method="jaccard", binary=FALSE, diag=FALSE, upper=TRUE, na.rm = FALSE)) # similarity
  
  pairs <- expand.grid(rownames(rede),rownames(rede)) # list of insect pairs
  pairs <- pairs[pairs$Var1!=pairs$Var2,] 
  
  for(j in 1:nrow(pairs)){ # for each pair
    prob <- PP[rownames(PP)==pairs$Var1[j],colnames(PP)==pairs$Var2[j]] # similarity
    exc.var1 <- which(M[rownames(M)==pairs$Var1[j],]==1&M[rownames(M)==pairs$Var2[j],]==0) # interactions exclusive to insect 1
    exc.var2 <- which(M[rownames(M)==pairs$Var2[j],]==1&M[rownames(M)==pairs$Var1[j],]==0) # interactions exclusive to insect 2
    rede[rownames(rede)==pairs$Var1[j], exc.var2] <- sample(c(1,0), length(exc.var2), prob=c(prob,1-prob), replace=TRUE) # int. exc. to ins. 2 copied by ins. 1
    rede[rownames(rede)==pairs$Var2[j], exc.var1] <- sample(c(1,0), length(exc.var1), prob=c(prob,1-prob), replace=TRUE) # int. exc. to ins. 1 copied by ins. 2
  }
  return(rede)
}	

#function to reassembly the network 
new_reconectar_rewiring_sin_dynamics<-function(nueva_comunidad, comunidad_anterior, red_anterior, fbmat){
	#'''function to reassembly the new plants and pollinators communities at each cicle allowing the rewiring of interactions
	#  based on a matrix of probabilities
	#input: nueva_comunidad=new plant community after simulation, comunidad_anterior=previous plant community (t(n-1)), 
	#red_anterior= previous network, fbmat= matrix of forbidden links    
	#output: a new plant-pollinator interactions network'''
	
	if (is.array(red_anterior) && (dim(red_anterior)[1]) >= 2 && (dim(red_anterior)[2])>=2){
	
		#estimates new pollinators frequency 
		nuevas_abundancias=colSums(red_anterior)
		
		if (length(nuevas_abundancias)<=2){
			red_final="small"
			#print("small1")
			
		}else{
		#estimates rewiring probabilities
		#for each pollinator, the prob of rewiring is proportionally inverse to the availability to their resouses (previous links) P=1-n 
		Prew=NULL 
		for (i in 1:ncol(red_anterior)){
			hosts=names(red_anterior[,i][red_anterior[,i]>0]) #which plant resources used the pollinator in the previuos time step??
			resources=sum(nueva_comunidad[hosts])   ### sum the resources availability in the new community
			Prew[i]=1-resources #estimates probability P=1-n 
		}
		names(Prew)=colnames(red_anterior)
		
		#combine with a vector of plants proportions in a prbability matrix
		Abmat=nueva_comunidad%*%t(Prew)
		rownames(Abmat)=names(nueva_comunidad)
		
		if (is.array(Abmat) && (dim(Abmat)[1]) >= 2 && (dim(Abmat)[2])>=2) {
			#reduce the fbmat (forbidden links matrix) to the subset of species we have in the new network 
			sub_fbmat=matrix(0,length(nueva_comunidad), length(nuevas_abundancias)) #empty matrix
			rownames(sub_fbmat)=names(nueva_comunidad) 
			colnames(sub_fbmat)=names(nuevas_abundancias)
			
			for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
				sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta 
				if (sp_planta %in% rownames(fbmat)){ #because not all plant species in nueva_comunidad are in fbmat
					for (j in 1:length(nuevas_abundancias)){
						sp_poli<- names(nuevas_abundancias)[j]#pollinators name
						sub_fbmat[i,j]=fbmat[sp_planta,sp_poli] #coopy cell value from fbmat to sub_fbmat
					}
				}
			}
			#Final probability matrix to re-assembly the network combining the previous matrices
			Pmat=(Abmat*sub_fbmat)/sum(Abmat*sub_fbmat)
			
			if (is.array(Pmat) && (dim(Pmat)[1]) >= 2 && (dim(Pmat)[2])>=2) {
			
				#re-assembly the network
				#option with mgen
				set.seed(NULL)
				L=sum(nuevas_abundancias) #number of links
				Pmat2=empty(Pmat)/sum(empty(Pmat))#normalize te matrix
				
				if (is.array(Pmat2) && (dim(Pmat2)[1]) >= 2 && (dim(Pmat2)[2])>=2) {
				#generate the new network
				red_final3=mgen(Pmat2,n=L,keep.species=TRUE, rep.cell=TRUE, autotransform="sum")##
				rownames(red_final3)=rownames(Pmat2)
				colnames(red_final3)=colnames(Pmat2)
				red_final3= empty(red_final3)#delete empty rows and columns
				
				if (is.array(red_final3) && (dim(red_final3)[1]) >= 2 && (dim(red_final3)[2])>=2) {
					red_final=red_final3
					
					}else{
					#print("small2")
					red_final="small"
					}
			}else{
			#print("small3")
			red_final="small"
			}
	}else{
			#print("small4")
			red_final= "small"}#cierra if Pmat
					
	}else{
	red_final= "small"
	#print("small5")
	
	}	#cierra if dimen de Abmat
				
	}
	}else{
	red_final="small"
	#print("small8")
	}
	
return(red_final)
}#

##################
####SIMULATION####
##################

setwd("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")
#setwd("C:\\Users\\User\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")

#Load plant species traits...
plant.traits<-read.table("planilla_parametros_FINAL.csv",sep=";", header=T)
plant.traits=plant.traits[1:142,]
#load herbicides traits
herbicides<- read.table("herbicides.csv",sep=",", header=T)
#load plant-pollinator interactions dataset
total<-read.table("inter_plant_pol.csv",sep=";", header=T)

####CREATE THE PLANT POLLINATOR META-WEB
#grouping spp
agrupar<-group_by(total,Plant,Pollinator)#agrupar filas iguales
#agrupar
#countingr
suma_filas<-summarise(agrupar,sum(frequency))
#suma_filas
metared=data.frame(suma_filas)

matriz_meta<-list2matr(metared$Plant,metared$Pollinator,metared$sum.frequency.,matr.type="q")
meta_web1<- as.matrix(as.data.frame(matriz_meta,col.names = names(metared$Pollinator))) 
meta_web<-meta_web1/sum(meta_web1)#INTERACTION MATRIX

####MATRIX OF FORBIDEN LINKS 
#forbidden links matrix based on pollinators similarity in the use of resources (to be used in reassembling networks)
fbmat=t(rew.rule(meta_web))#Based on MAia et al 2021 function. Transposed because the matrix is inversed 
#fbmat

escenario="SCE1" #set scenario "SCE1", SCE2, SCE3
simulations<-1000
#lists to save plant community simmulations and estimated network matrics
plant_pol_net_values<-list()
spp_plantas_simul=list()
todas_redes=list()
for (i in 1:simulations){ 

	#1) CHANGES ON PLANT COMMUNITIES 
	#1.1)choose a plant subset
		Nspp= 20 #number of species to be picked to create de initial plant community ##
		prop.p=plant.traits[,"proportion3"]#plant proportions###
		set.seed(NULL)
		rand.plant <- sample(1:142, Nspp, replace=FALSE, prob=prop.p) #sample a subset of plant species #
		#alternative option: all plant species with the same probability
		#rand.plant <- sample(1:142, Nspp, replace=FALSE) 
		
		#create the initial data subset to work with
		set.plants=plant.traits[rand.plant,]
		init.plant.comm<-set.plants$proportion/sum(set.plants$proportion)# Plant species proportions for the initial plant community
		names(init.plant.comm)<-set.plants$plant_species
		init.plant.comm  ##intial (t0) plant community
		
	#1.2 PLANT PARAMETERS FOR LOCAL EXTINCTION OR EXPANSION
		#ALTERNATIVA DE PARÁMETROS
		#weed_risk
		
		life.cyc<-(set.plants$cycle_value5-0.01)/max(set.plants$cycle_value5)
		banco<-(set.plants$seed_value2-0.01)/max(set.plants$seed_value2)
		dispers<-(set.plants$dispersal_value2-0.01)/max(set.plants$dispersal_value2)
		familia<-(set.plants$family_value2-0.01)/max(set.plants$family_value2)
				
		life.cyc_loss<-(set.plants$inv_cycle_value5-0.01)/max(set.plants$inv_cycle_value5)# 
		banco_loss<-(set.plants$inv_seed_value2-0.01)/max(set.plants$inv_seed_value2)
		familia_loss<-(set.plants$inv_family_value2-0.01)/max(set.plants$inv_family_value2) 
			
		weed_risk_perdida<- (life.cyc_loss*banco_loss*familia_loss)/max(life.cyc_loss*banco_loss*familia_loss)
		weed_risk_expansion<-(life.cyc*banco*dispers*familia)/max(life.cyc*banco*dispers*familia)
		
		#herbicide_risk
		biotypes1<-as.matrix(as.data.frame(set.plants[31:38],row.names = set.plants[,1]))
		biotypes=biotypes1/max(biotypes1)
		#herbicides characteristics and scenarios
		her<- as.matrix(as.data.frame(herbicides[,2:9],row.names = herbicides[,1]))
		presion<-(her["select_preasure",])/max(her["select_preasure",])
		dosis<-her[escenario,]/max(her[escenario,])#number of applications of each herbicide (different dosses) in 2 years 
		
		herb_risk= (rowSums(biotypes*presion*dosis))/max(rowSums(biotypes*presion*dosis))
		
		#Management
		###mmodifiers de Moss et al 2019
		if(escenario=="SCE3"){modifiers=1}#*1
		if(escenario=="SCE2"){modifiers=0.67}#0.67
		if(escenario=="SCE1"){modifiers=0.33}#0.33
		
		trait.losing<- ((weed_risk_perdida*(1-herb_risk))/max(weed_risk_perdida*(1-herb_risk)))*modifiers #modifiers va al final para dar valor a #todo el vector, sino se estandariza
		trait.expanding<-((weed_risk_expansion*herb_risk)/max(weed_risk_expansion*herb_risk))*modifiers
		
		## PLANT COMMUNITY MODEL (non-random starting community) ####
		## Time parameters
		years<-15 # how many years you want the model to run for  ###each scenario was designed for two years, so we have 15*2=30 years
		duration<-years##

		#plant simulation
		simulando<-plant_simulation(init.plant.comm,duration,trait.losing,trait.expanding)#,resistance) #here we obtain an array with plant proportions at each time step
		simulando1=rbind(as.vector(init.plant.comm),simulando) ##adding initial plant community to the array 
		spp_plantas_simul[[i]]=simulando1 #to save in the final list
		
		
		#plot curves
		#simu <- data.frame(x = seq_along(simulando1[, 1]),simulando1)
		#Formato long
		#simu <- melt(simu, id.vars = "x")
		#ggplot(simu, aes(x = x, y = value, color = variable)) +
		#geom_line(lwd=1.5)

	#2)NETWORK SIMULATION
	#2.1)INITIAL NETWORK
		#initial plant community
		comunidad_inicial1=as.vector(init.plant.comm)
		names(comunidad_inicial1)=names(init.plant.comm)
				
		#mean number of links by plots in empirical networks to create a random network.. mean 25, median 16
		N_links=25# to be used un network_assemblyV2 function	
		N_pol=15# to be used un network_assemblyV2 function	
		
		#create the initial random network using the initial plant community
		red_inicial=network_assemblyV2(meta_web, nueva_comunidad=comunidad_inicial1, N_pol)#selecting first a subset of pollinators and n_links=sum(web)*100
		#plotweb(red_inicial, arrow="down")
		#red_inicial
		
		#2.2) REASEMBLING NETWORKS AFTER SIMULATIONS
		##reassembling the networks at each cicle of "duration"
		red=list() #strore networks in a list
		for (k in 1:nrow(simulando)){
			comunidad_simulada=as.vector(simulando[k,])#select the corresponding plant community (tn)
			names(comunidad_simulada)=names(simulando[k,])
			
		#Reassembly the network	#using pollinators from previuos temporal network t(n-1).
			if (k==1){red[[k]]=new_reconectar_rewiring_sin_dynamics(nueva_comunidad=comunidad_simulada, comunidad_anterior=comunidad_inicial1, red_anterior=red_inicial, fbmat)} 
			if (k>1){
				if(any(red[[k-1]]=="small")){ 
				red[[k]]="small"
				}else{
					red[[k]]=new_reconectar_rewiring_sin_dynamics(nueva_comunidad=comunidad_simulada, comunidad_anterior=simulando[k-1,],red_anterior=red[[k-1]], fbmat)}
			}
			}
		#plotweb(red[[15]], arrow="down")
		#red

		
	#2.3)Estimate network metrics
		red_inicial1=list(red_inicial)
		redes=append(red_inicial1,red)###add initial network to the list of networks
		todas_redes[[i]]=redes
		
		plant_pol_net_values_cycle<-array(dim=c((years+1),18)) #storage array ###acá, para 30 años son 16 filas, para 50 años  es 26
		for (n in 1:length(redes)){ 
			
			#if (sum(redes[[n]]!=0)<=1){resumen_final=c(rep(NA,14), paste0("simul",k)) }#esto no va...
			
			if (any(redes[[n]]=="small")){resumen_final=cbind("small", NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, n, i)
			}else{

			#if (dim(redes[[n]])[2]>=2 && dim(redes[[n]])[1]>=2){#exc
				##metrics estimation
				#robustness
				extincion1<-second.extinct(redes[[n]], participant="both", method="abundance", nrep=100) #con both no funciona abundance
				robustez1=robustness(extincion1)
				robustez1
				
				#robustness - only plants
				extincion2<-second.extinct(redes[[n]], participant="lower", method="abundance", nrep=100) #abundance es de menos a mas #abundante
				robustez2=robustness(extincion2)
				robustez2
				
				dimen1<-dim(redes[[n]])
				S_plantas1=dimen1[1]
				S_poli1=dimen1[2]
				links1=sum(redes[[n]])
				evenn1 = networklevel(redes[[n]],index="interaction evenness",intereven="prod")##intereven="prod" or "sum"
				#evenn2 = networklevel(redes[[n]],index="interaction evenness", intereven="prod", effective = TRUE)##with hill numbers Jost 2010
				conn1 = networklevel(redes[[n]],index="connectance")
				nest1 = networklevel(redes[[n]], index="weighted NODF")
				nodf = networklevel(redes[[n]], index="NODF")#unweighted
				fun_com1 = fc(redes[[n]],dist="euclidean", method="average", weighted=TRUE)##### canberra?euclidea
				n_overHL = networklevel(redes[[n]],index="niche overlap")[1]
 				n_overLL = networklevel(redes[[n]],index="niche overlap")[2]
				int_div = networklevel(redes[[n]],index="Shannon diversity")
				gen = networklevel(redes[[n]],index="generality")[1]
				vul = networklevel(redes[[n]],index="generality")[2]
				H2 = H2fun(redes[[n]])[1]

				resumen_final=cbind(S_plantas1, S_poli1,links1, evenn1, conn1,nest1, nodf,robustez1,robustez2,fun_com1,n_overHL, n_overLL, int_div, gen, vul, H2, n, i)
							
			}
			plant_pol_net_values_cycle[n,]<-resumen_final
			}
		
		plant_pol_net_values[[i]]<-plant_pol_net_values_cycle
		
	}

#create a final table containing all the results
tabla_final=do.call(rbind, plant_pol_net_values) #deshago la lista
colnames(tabla_final)<- c("S_plantas", "S_poli","links", "evenn","conn","nest", "nodf","robust","robustPl", "funct_comp", "N_overHL", "N_overLL", "int_div", "gen", "vul", "H2","Sumulation","cycle")


tabla_final_modifiers=data.frame(tabla_final)
tabla_final_modifiers
#save table
write.table(tabla_final_modifiers, file = "G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\tabla_final_SC1_sindyn.csv", sep = ",", row.names=F)

#saving the frequency of each plant species in final simulations
spp_plantas_simul	
spp_finales=list()
for (z in 1:simulations){spp_finales[[z]]<-names(spp_plantas_simul[[z]][16,][spp_plantas_simul[[z]][16,]>0])}
frecuencias<-unlist(spp_finales)
suma_frecuencias=data.frame(table(frecuencias))	
suma_frec_ordenado <- suma_frecuencias[order(suma_frecuencias$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_plantas_SC1_sindyn.csv",sep=";")

#counting networks that became very small
small_net<-which(tabla_final_modifiers[,1]=="small")
length(small_net)

#saving frequency of each plant and pollinator species in final and initial simulations 
#todas_redes
poli_iniciales=list()
for (w in 1:simulations){poli_iniciales[[w]]<-colnames(todas_redes[[w]][[1]])}
frec_poli_inic<-unlist(poli_iniciales)
suma_frec_pol_in=data.frame(table(frec_poli_inic))
suma_frec_ordenado_po_i <- suma_frec_pol_in[order(suma_frec_pol_in$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_po_i,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_poli_iniciales_SC1_sindyn.csv",sep=";")

poli_finales=list()
for (w in 1:simulations){poli_finales[[w]]<-colnames(todas_redes[[w]][[16]])}
frec_poli_fin<-unlist(poli_finales)
suma_frec_pol_fin=data.frame(table(frec_poli_fin))
suma_frec_ordenado_po_f <- suma_frec_pol_fin[order(suma_frec_pol_fin$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_po_f,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_poli_finales_SC1_sindyn.csv",sep=";")

plantas_iniciales=list()
for (w in 1:simulations){plantas_iniciales[[w]]<-rownames(todas_redes[[w]][[1]])}
frec_plantas_inic<-unlist(plantas_iniciales)
suma_frec_pl_in=data.frame(table(frec_plantas_inic))
suma_frec_ordenado_pl_i <- suma_frec_pl_in[order(suma_frec_pl_in$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_pl_i,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_pl_iniciales_SC1_sindyn.csv",sep=";")

plantas_finales=list()
for (w in 1:simulations){plantas_finales[[w]]<-rownames(todas_redes[[w]][[16]])}
frec_plantas_fin<-unlist(plantas_finales)
suma_frec_pl_fin=data.frame(table(frec_plantas_fin))
suma_frec_ordenado_pl_f <- suma_frec_pl_fin[order(suma_frec_pl_fin$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado_pl_f,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_pl_finales_SC1_sindyn.csv",sep=";")


####################################################################################################################################################################
########################################
########################
##########################
#SIN modifiers
####Modelling the effect of herbicide application on plant-pollinator networks in agroecosystems 
###Julia Tavella and Fred Windsor
##Preparing the environment
#install.packages(c("RColorBrewer","tidyverse","sads","bipartite","igraph","gplots","vegan","reshape"))
library(RColorBrewer)
library(tidyverse) 
library(sads) 
library(bipartite)
library(igraph)
library(gplots)
library(vegan)
library(reshape)

###First, some functions...
##list to matrix function
list2matr<-function(sp1,sp2,effects,matr.type="q"){
  nodelist1u<-array(unique(sp1))
  nodelist2u<-array(unique(sp2))
  m<-length(nodelist1u)
  n<-length(nodelist2u)
  matr<-matrix(0,m,n) #Matrix
  colnames(matr)<-nodelist2u   #######
  row.names(matr)<-nodelist1u
  for (i in 1:m){
    iu<-nodelist1u[i]
    ii<-which(sp1==iu)
    nodelist2ii<-sp2[ii]
    link.weightii<-effects[ii]
    for (j in 1:length(ii)){
      jj<-which(nodelist2u==nodelist2ii[j])
      if (matr.type=="b") matr[i,jj]<-1 #Binary matrix
      if (matr.type=="q") matr[i,jj]<-link.weightii[j] #Weighted matrix
    }
  }
  res<-list(matr)
  names(res)<-c("matr")
  return(res)
}


################NUEVA PLANT SIMULATION	11/10/2023
####Function for Plant community changes simulation##
plant_simulation<-function(init.plant.comm,duration,trait.losing,trait.expanding){#,resistance){ 
	#'''simulate losses and expansion of plant species
	#	input:a vector of an initial plant comunity, duration=number of cicles, vectors of probabilities of loss, expansion and resistance preasure
	#	output:an array containing the community composition at each time step (duration)
	#	'''
	
	## Storage arrays ####
	plant.comm_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for the data generated from the model
	plant.loss_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for coverage loss data generated from the model
	plant.exp_array<-array(dim=c(duration,length(init.plant.comm))) # A storage array for coverage expansion data generated from the model

	for (j in 1:duration){ # loop through each time step
								
		# Select the starting plant community
		if(j<=1){plant.comm=init.plant.comm} # set the plant community as the starting coverage if the time step is 1
		if(j>1){plant.comm=plant.comm_array[j-1,]} # if not, use the plant community from the previous time step
		
		##ESTE NUNCA TIENE QUE SER NAN		
		# calculate the loss of plant coverage due to herbicide application and resistance  #
		plant.loss_array[j,]<-plant.comm-(plant.comm*trait.losing)#lo que queda de la comunidad anterior
		
	    #set a threshold, a sensible proportion at wich we would not expect a plant population that covers <0.01 of the transect will be viable
		plant.loss_array[j,][plant.loss_array[j,]<0.001]<-0 ##fernando recomends 0.1 o 0.5. #!!! We have to set a real threshold here!
		#if 1 is the space we have, 0.01 is the 1% of area occupied by the plant. estoy usando 0.001
		
		space<-1-sum(plant.loss_array[j,]) # calculate the space available due to the loss of plants
			
		# calculate the proportional ability of plant species to expand into the space 
		# before that, I eliminate plant species that we lose at the previous step
		surviving.plants=plant.loss_array[j,] 
		
		surviving.plants[surviving.plants>0]<-1 #binary vector to "delete" lost species from the expanding vector 
		expanding1=trait.expanding * surviving.plants ## parameter of expansion 
		
		if (sum(expanding1)>0){plant.exp.ability<-(expanding1)/sum(expanding1)}
		if (sum(expanding1)==0){plant.exp.ability<- expanding1} #to avoid NaNs
		
		#Calculate the expansion of plants (more resistant plants more expantion chances)	
		plant.exp_array[j,]<-space*plant.exp.ability ##
		
		# calculate the change in plant community coverage after loss and expansion
		if (sum(plant.exp_array[j,])>0){plant.comm_array[j,]<-(plant.loss_array[j,]+plant.exp_array[j,])/sum(plant.loss_array[j,]+plant.exp_array[j,])}
		if (sum(plant.exp_array[j,])==0){plant.comm_array[j,]<-(plant.loss_array[j,]+plant.exp_array[j,])}
	}
	
	colnames(plant.comm_array)<-names(init.plant.comm)
	return(plant.comm_array)  #### return an array with simulations equal to duration 
}


######NETWORK ASSEMBLY function to construct the initial network###
network_assemblyV2<-function(meta_web, nueva_comunidad, N_pol){ 
	#"""Function to generate networks using an initial plant community, and 
	#based on the links with pollinators from the meta-web
	#I used mgen function (in bipartite) to genrate a random network using a probability matrix (based on pollinator frequency and plant species proportions)
	#input: meta_web, nueva_comunidad (initial plant community or plant community at each step of the simulation), N_links (mean number of links  
	#estimated from empirical networks by plot)
	#output: a plant-pollinator network based on the probability matrix and preserving the proportion of links 
	#""" 

	#select, previously, a subset of pollinators of size equal to the mean number of polli in empirical networks
	#give probabilities according to species frequency in meta_Web
	set.seed(NULL)
	polis=sample(1:ncol(meta_web),N_pol,prob=colSums(meta_web),replace=FALSE)
	nomb_polis=colnames(meta_web)[polis]
	
	#create a matrix containing plants from nueva_comunidad interacting with a the subset of pollinators, with values of interactions available in th meta-web
	sub_meta_web=matrix(0,length(nueva_comunidad), N_pol)#storage matrix
	rownames(sub_meta_web)=names(nueva_comunidad) 
	colnames(sub_meta_web)=nomb_polis
	
	for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
		sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta
		if (sp_planta %in% rownames(meta_web)){ #because not all plant species are in the meta_web
			for (j in 1:N_pol){
				sp_poli<- nomb_polis[j]
				sub_meta_web[i,j]=meta_web[sp_planta,sp_poli] #select the value of the interaction in meta_web to copy in sub_meta_web
			}
		}
	}
	
	#define pollinators relative frequency 
	suma_col=colSums(as.matrix(sub_meta_web))###
	pol.frec=suma_col/sum(suma_col)
		
	#define proportion of entomophilous plant species (deleting plants that not interact with pollinators from nueva_comunidad)
	nueva_comunidad2=nueva_comunidad[nueva_comunidad>0] 
	plant.prop=nueva_comunidad2/sum(nueva_comunidad2)

	#interaction probabilities matrix 
	sub_meta_web1=sub_meta_web
	sub_meta_web1[sub_meta_web1>=1]=1##forbidden links matrix.presence-absence links to give probabily 0 to interactions that do not exist in meta_web
	pmat=((plant.prop%*%t(pol.frec))*sub_meta_web1)/sum((plant.prop%*%t(pol.frec))*sub_meta_web1)#construct the matrix of probabilities, combining vectors and forbidden link mat
	rownames(pmat)=names(nueva_comunidad2)
	pmat= empty(pmat)#delete empty rows and columns if necesary
	
	#new random network based on pmat
	#number of links equal to sub_meta_web (proportion of links), all species remain in the network, quantitative network
	#Setting keep.species to FALSE may (but need not) cause a loss of species and networks have a small number of species. With keep.species=TRUE 
	#all selected species remain in the network. Results are quite different.
		set.seed(NULL)
		L=sum(sub_meta_web)*100
		nueva_red=mgen(pmat,n=L,keep.species=TRUE, rep.cell=TRUE, autotransform="sum")### Alternative:create the network usin "sample" function (is quite simmilar)
		rownames(nueva_red)=rownames(pmat)
		colnames(nueva_red)=colnames(pmat)
	
	return(nueva_red)
}

##function for interaction rewiring probabilities
rew.rule <- function(M){
	#function to build a forbidden links matrix based on similarity of interacctions between species (from Maia et al 2021)
	#Input: MAtrix of interactions
	#output: Matrix built based on similarity between species interactions
  rede <- t(M) #aca transpuse , tambien deberia transponer la salida
  rede <- as.matrix(rede)
  
  PP <- 1-as.matrix(vegdist(rede, method="jaccard", binary=FALSE, diag=FALSE, upper=TRUE, na.rm = FALSE)) # similarity
  
  pairs <- expand.grid(rownames(rede),rownames(rede)) # list of insect pairs
  pairs <- pairs[pairs$Var1!=pairs$Var2,] 
  
  for(j in 1:nrow(pairs)){ # for each pair
    prob <- PP[rownames(PP)==pairs$Var1[j],colnames(PP)==pairs$Var2[j]] # similarity
    exc.var1 <- which(M[rownames(M)==pairs$Var1[j],]==1&M[rownames(M)==pairs$Var2[j],]==0) # interactions exclusive to insect 1
    exc.var2 <- which(M[rownames(M)==pairs$Var2[j],]==1&M[rownames(M)==pairs$Var1[j],]==0) # interactions exclusive to insect 2
    rede[rownames(rede)==pairs$Var1[j], exc.var2] <- sample(c(1,0), length(exc.var2), prob=c(prob,1-prob), replace=TRUE) # int. exc. to ins. 2 copied by ins. 1
    rede[rownames(rede)==pairs$Var2[j], exc.var1] <- sample(c(1,0), length(exc.var1), prob=c(prob,1-prob), replace=TRUE) # int. exc. to ins. 1 copied by ins. 2
  }
  return(rede)
}	

#changes in pollinators
pollinators_dynamV2<-function(nueva_comunidad, comunidad_anterior, red_anterior){
#'''function that simulates changes on abundance of pollinators species due to changes on plant community
	#input: nueva_comunidad=new plant community after simulation, comunidad_anterior=previous plant community (t(n-1)), 
	#red_anterior= previous network
	#output: a vector containing the new pollinator frequencies
	#'''
	#changes on plant proportions
	popul_changes=nueva_comunidad-comunidad_anterior #difference between previous and current proportions of each plant species
	#pollinators total frequency (num of links) in the network t(n-1)
	prev_pol=colSums(red_anterior)
	
	#dynamics on pollinator populations
	red_anterior[red_anterior>0]=1 #left 0 and 1 to multiply with popul changes
	
	new_ab=NULL###matrix(0,dim(red_anterior)[1],dim(red_anterior)[2])
	for (i in 1:length(prev_pol)){ #for each pollinator species in the previous network
		if (dim(red_anterior)[1]==1 & dim(red_anterior)[2]==1 ){
		hosts=rownames(red_anterior)
		}else{
		hosts=names(red_anterior[,i][red_anterior[,i]>0])		#select its host plants
		}
		if (sum(comunidad_anterior[hosts])>0){resources=sum(popul_changes[hosts])/sum(comunidad_anterior[hosts])}#estimates the rate of change in the proportion of pollinator resources
		if (sum(comunidad_anterior[hosts])==0){resources=sum(popul_changes[hosts])}
		
		if (resources>1){abun_change= sum(rbinom(prev_pol[i],1,prob=1))}#if resources increase, pollinators abundance increase (prob=1)
		if (resources< 0){abun_change= -(sum(rbinom(prev_pol[i],1,prob=abs(resources))))} #if we lost resources we substract abun_change (we loose pollin)
		if (resources >=0 && resources <=1){abun_change= sum(rbinom(prev_pol[i],1,prob=abs(resources)))} #if resources increase, the number of pollinators increase
		
		new_ab[i]=prev_pol[i] + abun_change
		}
	names(new_ab)=names(prev_pol)
	return(new_ab)
	}
	
#function to reassembly the network 
new_reconectar_rewiring<-function(nueva_comunidad, comunidad_anterior, red_anterior, fbmat){
	#'''function to reassembly the new plants and pollinators communities at each cicle allowing the rewiring of interactions
	#  based on a matrix of probabilities
	#input: nueva_comunidad=new plant community after simulation, comunidad_anterior=previous plant community (t(n-1)), 
	#red_anterior= previous network, fbmat= matrix of forbidden links    
	#output: a new plant-pollinator interactions network'''
	
	if (is.array(red_anterior) && (dim(red_anterior)[1]) >= 2 && (dim(red_anterior)[2])>=2){
	
		#estimates new pollinators frequency 
		nuevas_abundancias=pollinators_dynamV2(nueva_comunidad, comunidad_anterior, red_anterior)
		
		if (length(nuevas_abundancias)<=2){
			red_final="small"
			#print("small1")
			
		}else{
		#estimates rewiring probabilities
		#for each pollinator, the prob of rewiring is proportionally inverse to the availability to their resouses (previous links) P=1-n 
		Prew=NULL 
		for (i in 1:ncol(red_anterior)){
			hosts=names(red_anterior[,i][red_anterior[,i]>0]) #which plant resources used the pollinator in the previuos time step??
			resources=sum(nueva_comunidad[hosts])   ### sum the resources availability in the new community
			Prew[i]=1-resources #estimates probability P=1-n 
		}
		names(Prew)=colnames(red_anterior)
		
		#combine with a vector of plants proportions in a prbability matrix
		Abmat=nueva_comunidad%*%t(Prew)
		rownames(Abmat)=names(nueva_comunidad)
		
		if (is.array(Abmat) && (dim(Abmat)[1]) >= 2 && (dim(Abmat)[2])>=2) {
			#reduce the fbmat (forbidden links matrix) to the subset of species we have in the new network 
			sub_fbmat=matrix(0,length(nueva_comunidad), length(nuevas_abundancias)) #empty matrix
			rownames(sub_fbmat)=names(nueva_comunidad) 
			colnames(sub_fbmat)=names(nuevas_abundancias)
			
			for (i in 1:length(nueva_comunidad)){ #para cada planta de la nueva camunidad
				sp_planta<-as.character(names(nueva_comunidad[i]))#nombre de la planta 
				if (sp_planta %in% rownames(fbmat)){ #because not all plant species in nueva_comunidad are in fbmat
					for (j in 1:length(nuevas_abundancias)){
						sp_poli<- names(nuevas_abundancias)[j]#pollinators name
						sub_fbmat[i,j]=fbmat[sp_planta,sp_poli] #coopy cell value from fbmat to sub_fbmat
					}
				}
			}
			#Final probability matrix to re-assembly the network combining the previous matrices
			Pmat=(Abmat*sub_fbmat)/sum(Abmat*sub_fbmat)
			
			if (is.array(Pmat) && (dim(Pmat)[1]) >= 2 && (dim(Pmat)[2])>=2) {
			
				#re-assembly the network
				#option with mgen
				set.seed(NULL)
				L=sum(nuevas_abundancias) #number of links
				Pmat2=empty(Pmat)/sum(empty(Pmat))#normalize te matrix
				
				if (is.array(Pmat2) && (dim(Pmat2)[1]) >= 2 && (dim(Pmat2)[2])>=2) {
				#generate the new network
				red_final3=mgen(Pmat2,n=L,keep.species=TRUE, rep.cell=TRUE, autotransform="sum")##
				rownames(red_final3)=rownames(Pmat2)
				colnames(red_final3)=colnames(Pmat2)
				red_final3= empty(red_final3)#delete empty rows and columns
				
				if (is.array(red_final3) && (dim(red_final3)[1]) >= 2 && (dim(red_final3)[2])>=2) {
					red_final=red_final3
					
					}else{
					#print("small2")
					red_final="small"
					}
			}else{
			#print("small3")
			red_final="small"
			}
	}else{
			#print("small4")
			red_final= "small"}#cierra if Pmat
					
	}else{
	red_final= "small"
	#print("small5")
	
	}	#cierra if dimen de Abmat
				
	}
	}else{
	red_final="small"
	#print("small8")
	}
	
return(red_final)
}#

##################
####SIMULATION####
##################

setwd("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")
#setwd("C:\\Users\\User\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")

#Load plant species traits...
plant.traits<-read.table("planilla_parametros_FINAL.csv",sep=";", header=T)
plant.traits=plant.traits[1:142,]
#load herbicides traits
herbicides<- read.table("herbicides.csv",sep=",", header=T)
#load plant-pollinator interactions dataset
total<-read.table("inter_plant_pol.csv",sep=";", header=T)

####CREATE THE PLANT POLLINATOR META-WEB
#grouping spp
agrupar<-group_by(total,Plant,Pollinator)#agrupar filas iguales
#agrupar
#countingr
suma_filas<-summarise(agrupar,sum(frequency))
#suma_filas
metared=data.frame(suma_filas)

matriz_meta<-list2matr(metared$Plant,metared$Pollinator,metared$sum.frequency.,matr.type="q")
meta_web1<- as.matrix(as.data.frame(matriz_meta,col.names = names(metared$Pollinator))) 
meta_web<-meta_web1/sum(meta_web1)#INTERACTION MATRIX

####MATRIX OF FORBIDEN LINKS 
#forbidden links matrix based on pollinators similarity in the use of resources (to be used in reassembling networks)
fbmat=t(rew.rule(meta_web))#Based on MAia et al 2021 function. Transposed because the matrix is inversed 
#fbmat

escenario="SCE1" #set scenario "SCE1", SCE2, SCE3
simulations<-1000
#lists to save plant community simmulations and estimated network matrics
plant_pol_net_values<-list()
spp_plantas_simul=list()
todas_redes=list()
for (i in 1:simulations){ 

	#1) CHANGES ON PLANT COMMUNITIES 
	#1.1)choose a plant subset
		Nspp= 20 #number of species to be picked to create de initial plant community ##
		prop.p=plant.traits[,"proportion3"]#plant proportions###
		set.seed(NULL)
		rand.plant <- sample(1:142, Nspp, replace=FALSE, prob=prop.p) #sample a subset of plant species #
		#alternative option: all plant species with the same probability
		#rand.plant <- sample(1:142, Nspp, replace=FALSE) 
		
		#create the initial data subset to work with
		set.plants=plant.traits[rand.plant,]
		init.plant.comm<-set.plants$proportion/sum(set.plants$proportion)# Plant species proportions for the initial plant community
		names(init.plant.comm)<-set.plants$plant_species
		init.plant.comm  ##intial (t0) plant community
		
	#1.2 PLANT PARAMETERS FOR LOCAL EXTINCTION OR EXPANSION
		#ALTERNATIVA DE PARÁMETROS
		#weed_risk
		
		life.cyc<-(set.plants$cycle_value5-0.01)/max(set.plants$cycle_value5)
		banco<-(set.plants$seed_value2-0.01)/max(set.plants$seed_value2)
		dispers<-(set.plants$dispersal_value2-0.01)/max(set.plants$dispersal_value2)
		familia<-(set.plants$family_value2-0.01)/max(set.plants$family_value2)
				
		life.cyc_loss<-(set.plants$inv_cycle_value5-0.01)/max(set.plants$inv_cycle_value5)# 
		banco_loss<-(set.plants$inv_seed_value2-0.01)/max(set.plants$inv_seed_value2)
		familia_loss<-(set.plants$inv_family_value2-0.01)/max(set.plants$inv_family_value2) 
			
		weed_risk_perdida<- (life.cyc_loss*banco_loss*familia_loss)/max(life.cyc_loss*banco_loss*familia_loss)
		weed_risk_expansion<-(life.cyc*banco*dispers*familia)/max(life.cyc*banco*dispers*familia)
		
		#herbicide_risk
		biotypes1<-as.matrix(as.data.frame(set.plants[31:38],row.names = set.plants[,1]))
		biotypes=biotypes1/max(biotypes1)
		#herbicides characteristics and scenarios
		her<- as.matrix(as.data.frame(herbicides[,2:9],row.names = herbicides[,1]))
		presion<-(her["select_preasure",])/max(her["select_preasure",])
		dosis<-her[escenario,]/max(her[escenario,])#number of applications of each herbicide (different dosses) in 2 years 
		
		herb_risk= (rowSums(biotypes*presion*dosis))/max(rowSums(biotypes*presion*dosis))
		
		#Management
		###mmodifiers de Moss et al 2019
		#if(escenario=="SCE3"){modifiers=1}#*1
		#if(escenario=="SCE2"){modifiers=0.67}#0.67
		#if(escenario=="SCE1"){modifiers=0.33}#0.33
		
		trait.losing<- ((weed_risk_perdida*(1-herb_risk))/max(weed_risk_perdida*(1-herb_risk)))#*modifiers #modifiers va al final para dar valor a #todo el vector, sino se estandariza
		trait.expanding<-((weed_risk_expansion*herb_risk)/max(weed_risk_expansion*herb_risk))#*modifiers
		
		## PLANT COMMUNITY MODEL (non-random starting community) ####
		## Time parameters
		years<-15 # how many years you want the model to run for  ###each scenario was designed for two years, so we have 15*2=30 years
		duration<-years##

		#plant simulation
		simulando<-plant_simulation(init.plant.comm,duration,trait.losing,trait.expanding)#,resistance) #here we obtain an array with plant proportions at each time step
		simulando1=rbind(as.vector(init.plant.comm),simulando) ##adding initial plant community to the array 
		spp_plantas_simul[[i]]=simulando1 #to save in the final list
		
		
		#plot curves
		#simu <- data.frame(x = seq_along(simulando1[, 1]),simulando1)
		#Formato long
		#simu <- melt(simu, id.vars = "x")
		#ggplot(simu, aes(x = x, y = value, color = variable)) +
		#geom_line(lwd=1.5)

	#2)NETWORK SIMULATION
	#2.1)INITIAL NETWORK
		#initial plant community
		comunidad_inicial1=as.vector(init.plant.comm)
		names(comunidad_inicial1)=names(init.plant.comm)
				
		#mean number of links by plots in empirical networks to create a random network.. mean 25, median 16
		N_links=25# to be used un network_assemblyV2 function	
		N_pol=15# to be used un network_assemblyV2 function	
		
		#create the initial random network using the initial plant community
		red_inicial=network_assemblyV2(meta_web, nueva_comunidad=comunidad_inicial1, N_pol)#selecting first a subset of pollinators and n_links=sum(web)*100
		#plotweb(red_inicial, arrow="down")
		#red_inicial
		
		#2.2) REASEMBLING NETWORKS AFTER SIMULATIONS
		##reassembling the networks at each cicle of "duration"
		red=list() #strore networks in a list
		for (k in 1:nrow(simulando)){
			comunidad_simulada=as.vector(simulando[k,])#select the corresponding plant community (tn)
			names(comunidad_simulada)=names(simulando[k,])
			
		#Reassembly the network	#using pollinators from previuos temporal network t(n-1).
			if (k==1){red[[k]]=new_reconectar_rewiring(nueva_comunidad=comunidad_simulada, comunidad_anterior=comunidad_inicial1, red_anterior=red_inicial, fbmat)} 
			if (k>1){
				if(any(red[[k-1]]=="small")){ 
				red[[k]]="small"
				}else{
					red[[k]]=new_reconectar_rewiring(nueva_comunidad=comunidad_simulada, comunidad_anterior=simulando[k-1,],red_anterior=red[[k-1]], fbmat)}
			}
			}
		#plotweb(red[[15]], arrow="down")
		#red

		
	#2.3)Estimate network metrics
		red_inicial1=list(red_inicial)
		redes=append(red_inicial1,red)###add initial network to the list of networks
		todas_redes[[i]]=redes
		
		plant_pol_net_values_cycle<-array(dim=c((years+1),18)) #storage array ###acá, para 30 años son 16 filas, para 50 años  es 26
		for (n in 1:length(redes)){ 
			
			#if (sum(redes[[n]]!=0)<=1){resumen_final=c(rep(NA,14), paste0("simul",k)) }#esto no va...
			
			if (any(redes[[n]]=="small")){resumen_final=cbind("small", NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, n, i)
			}else{

			#if (dim(redes[[n]])[2]>=2 && dim(redes[[n]])[1]>=2){#exc
				##metrics estimation
				#robustness
				extincion1<-second.extinct(redes[[n]], participant="both", method="abundance", nrep=100) #con both no funciona abundance
				robustez1=robustness(extincion1)
				robustez1
				
				#robustness - only plants
				extincion2<-second.extinct(redes[[n]], participant="lower", method="abundance", nrep=100) #abundance es de menos a mas #abundante
				robustez2=robustness(extincion2)
				robustez2
				
				dimen1<-dim(redes[[n]])
				S_plantas1=dimen1[1]
				S_poli1=dimen1[2]
				links1=sum(redes[[n]])
				evenn1 = networklevel(redes[[n]],index="interaction evenness",intereven="prod")##intereven="prod" or "sum"
				#evenn2 = networklevel(redes[[n]],index="interaction evenness", intereven="prod", effective = TRUE)##with hill numbers Jost 2010
				conn1 = networklevel(redes[[n]],index="connectance")
				nest1 = networklevel(redes[[n]], index="weighted NODF")
				nodf = networklevel(redes[[n]], index="NODF")#unweighted
				fun_com1 = fc(redes[[n]],dist="euclidean", method="average", weighted=TRUE)##### canberra?euclidea
				n_overHL = networklevel(redes[[n]],index="niche overlap")[1]
 				n_overLL = networklevel(redes[[n]],index="niche overlap")[2]
				int_div = networklevel(redes[[n]],index="Shannon diversity")
				gen = networklevel(redes[[n]],index="generality")[1]
				vul = networklevel(redes[[n]],index="generality")[2]
				H2 = H2fun(redes[[n]])[1]

				resumen_final=cbind(S_plantas1, S_poli1,links1, evenn1, conn1,nest1, nodf,robustez1,robustez2,fun_com1,n_overHL, n_overLL, int_div, gen, vul, H2, n, i)
							
			}
			plant_pol_net_values_cycle[n,]<-resumen_final
			}
		
		plant_pol_net_values[[i]]<-plant_pol_net_values_cycle
		
	}

#create a final table containing all the results
tabla_final=do.call(rbind, plant_pol_net_values) #deshago la lista
colnames(tabla_final)<- c("S_plantas", "S_poli","links", "evenn","conn","nest", "nodf","robust","robustPl", "funct_comp", "N_overHL", "N_overLL", "int_div", "gen", "vul", "H2","Sumulation","cycle")

tabla_final_modifiers=data.frame(tabla_final)
tabla_final_modifiers
#save table
write.table(tabla_final_modifiers, file = "G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\tabla_final_SC1_SINMODIF.csv", sep = ",", row.names=F)

#saving the frequency of each plant species in final simulations
spp_plantas_simul	
spp_finales=list()
for (z in 1:simulations){spp_finales[[z]]<-names(spp_plantas_simul[[z]][16,][spp_plantas_simul[[z]][16,]>0])}
frecuencias<-unlist(spp_finales)
suma_frecuencias=data.frame(table(frecuencias))	
suma_frec_ordenado <- suma_frecuencias[order(suma_frecuencias$Freq,decreasing = TRUE), ]
write.table(suma_frec_ordenado,"G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\lista_final_frecuencia_plantas_SC1_SINMODIF.csv",sep=";")

#counting networks that became very small
small_net<-which(tabla_final_modifiers[,1]=="small")
length(small_net)
















############################################################################################################################################################################################################################################################################################################################################################################################################################################################


###PLOTS###
library(ggplot2)
library(tidyverse)
library(igraph)
library(Cairo)



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


#load data and combine tables
setwd("G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1")
#setwd("C:\\Users\\User\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET")

SC1<-read.table("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET\\MODIFIERS\\junio2024\\tabla_final_SC1.csv", sep=",", header=T) 
scenario=rep("SC1",nrow(SC1))
SC11<-cbind(scenario,SC1)

SC2<-read.table("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET\\MODIFIERS\\junio2024\\tabla_final_SC2.csv", sep=",", header=T) 
scenario=rep("SC2",nrow(SC2))
SC21<-cbind(scenario,SC2)

SC3<-read.table("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET\\MODIFIERS\\junio2024\\tabla_final_SC3.csv", sep=",", header=T) 
scenario=rep("SC3",nrow(SC3))
SC31<-cbind(scenario,SC3)

#tabla_final<-rbind(SC11,SC21,SC31)

NORES<-read.table("C:\\Users\\Usuario\\Dropbox\\data redes herbicidas planta herbivoro\\PAPER EN MARCHA\\DATA SET\\MODIFIERS\\junio2024\\tabla_final_random.csv", sep=",", header=T) 
scenario=rep("NORES",nrow(NORES))
NORES1<-cbind(scenario,NORES)

NOREWSC1<-read.table("tabla_final_SC1_norew.csv", sep=",", header=T) #veo solo sc1 para comparar con los otros modelos
scenario=rep("NOREWSC1",nrow(NOREWSC1))
NOREWSC11<-cbind(scenario,NOREWSC1)
colnames(NOREWSC11)=colnames(SC11)

NOREWSC2<-read.table("tabla_final_SC2_norew.csv", sep=",", header=T) #veo solo sc1 para comparar con los otros modelos
scenario=rep("NOREWSC2",nrow(NOREWSC2))
NOREWSC21<-cbind(scenario,NOREWSC2)
colnames(NOREWSC21)=colnames(SC11)

NOREWSC3<-read.table("tabla_final_SC3_norew.csv", sep=",", header=T) #veo solo sc1 para comparar con los otros modelos
scenario=rep("NOREWSC3",nrow(NOREWSC3))
NOREWSC31<-cbind(scenario,NOREWSC3)
colnames(NOREWSC31)=colnames(SC11)

NODYNSC1<-read.table("tabla_final_SC1_nodyn.csv", sep=",", header=T) 
scenario=rep("NODYNSC1",nrow(NODYNSC1))
NODYNSC11<-cbind(scenario,NODYNSC1)
colnames(NODYNSC11)=colnames(SC11)

NODYNSC2<-read.table("tabla_final_SC2_nodyn.csv", sep=",", header=T) 
scenario=rep("NODYNSC2",nrow(NODYNSC2))
NODYNSC21<-cbind(scenario,NODYNSC2)
colnames(NODYNSC21)=colnames(SC11)

NODYNSC3<-read.table("tabla_final_SC3_nodyn.csv", sep=",", header=T) 
scenario=rep("NODYNSC3",nrow(NODYNSC3))
NODYNSC31<-cbind(scenario,NODYNSC3)
colnames(NODYNSC31)=colnames(SC11)

NODYN_METASC1<-read.table("tabla_final_SC1_nodyn_meta.csv", sep=",", header=T) 
scenario=rep("NODYN_METASC1",nrow(NODYN_METASC1))
NODYN_METASC11<-cbind(scenario,NODYN_METASC1)
colnames(NODYN_METASC11)=colnames(SC11)

NODYN_METASC2<-read.table("tabla_final_SC2_nodyn_meta.csv", sep=",", header=T) 
scenario=rep("NODYN_METASC2",nrow(NODYN_METASC2))
NODYN_METASC21<-cbind(scenario,NODYN_METASC2)
colnames(NODYN_METASC21)=colnames(SC21)

NODYN_METASC3<-read.table("tabla_final_SC3_nodyn_meta.csv", sep=",", header=T) 
scenario=rep("NODYN_METASC3",nrow(NODYN_METASC3))
NODYN_METASC31<-cbind(scenario,NODYN_METASC3)
colnames(NODYN_METASC31)=colnames(SC11)


SINDYNSC1<-read.table("tabla_final_SC1_sindyn.csv", sep=",", header=T) 
scenario=rep("SINDYNSC1",nrow(SINDYNSC1))
SINDYNSC11<-cbind(scenario,SINDYNSC1)
colnames(SINDYNSC11)=colnames(SC11)

SINDYNSC2<-read.table("tabla_final_SC2_sindyn.csv", sep=",", header=T) 
scenario=rep("SINDYNSC2",nrow(SINDYNSC2))
SINDYNSC21<-cbind(scenario,SINDYNSC2)
colnames(SINDYNSC21)=colnames(SC11)

SINDYNSC3<-read.table("tabla_final_SC3_sindyn.csv", sep=",", header=T) 
scenario=rep("SINDYNSC3",nrow(SINDYNSC3))
SINDYNSC31<-cbind(scenario,SINDYNSC3)
colnames(SINDYNSC31)=colnames(SC11)

SINMODIFSC1<-read.table("tabla_final_SC1_SINMODIF.csv", sep=",", header=T) 
scenario=rep("SINMODIFSC1",nrow(SINMODIFSC1))
SINMODIFSC11<-cbind(scenario,SINMODIFSC1)
colnames(SINMODIFSC11)=colnames(SC11)

SINMODIFSC2<-read.table("tabla_final_SC2_SINMODIF.csv", sep=",", header=T) 
scenario=rep("SINMODIFSC2",nrow(SINMODIFSC2))
SINMODIFSC21<-cbind(scenario,SINMODIFSC2)
colnames(SINMODIFSC21)=colnames(SC11)

SINMODIFSC3<-read.table("tabla_final_SC3_SINMODIF.csv", sep=",", header=T) 
scenario=rep("SINMODIFSC3",nrow(SINMODIFSC3))
SINMODIFSC31<-cbind(scenario,SINMODIFSC3)
colnames(SINMODIFSC31)=colnames(SC11)



tabla_final<-rbind(SC11,SC21, SC31, NORES1, NOREWSC11,NOREWSC21,NOREWSC31, NODYNSC11, SINDYNSC11, SINDYNSC21, SINDYNSC31,NODYNSC21,NODYNSC31,NODYN_METASC11,NODYN_METASC21,NODYN_METASC31, SINMODIFSC11, SINMODIFSC21, SINMODIFSC31)

tabla_final<-rbind(SC11,SC21, SC31, SINMODIFSC11, SINMODIFSC21, SINMODIFSC31)





###############}
###REPLACING "SMALL" WITH 1 (PLANTS AND POLLI)
tabla_final$S_poli<-replace(tabla_final[,3], tabla_final[,2]=="small", 1)
tabla_final$S_plantas<-as.numeric(replace(tabla_final[,2], tabla_final[,2]=="small", 1))


###"S_plantas", "S_poli"
#summary for each scenario and simulation
tgc <- summarySE(tabla_final, measurevar="S_poli",groupvars=c("scenario","Sumulation"))###, groupvars=c("supp","dose"))
#tgc

simulaciones=as.numeric(tgc$Sumulation)-1
tgc2<-cbind(tgc,simulaciones)

#Line graphs
#After the data is summarized, we can make the graph. These are basic line and point graph with error bars representing either the standard error of the mean, or 95% confidence interval.

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.15) # move them .05 to the left and right

##A finished graph with IC of the mean might look like this. The points are drawn last so that the white fill goes on top of the lines and error bars.
#colore<-c("lightgoldenrod3","goldenrod3","indianred3")
AA<-ggplot(tgc2, aes(x=as.numeric(simulaciones), y=S_poli,colour=scenario)) + 
    geom_errorbar(aes(ymin=S_poli-ci, ymax=S_poli+ci), width=.2, position=pd)+#, colour="black") +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=16) + # 21 is filled circle, fill="white")
	scale_color_manual(values = c("SC1" = "lightgoldenrod3", "SC2" = "goldenrod3", "SC3"  = "indianred3", "SINMODIFSC1"= "olivedrab2", "SINMODIFSC2" = "olivedrab3","SINMODIFSC3" = "olivedrab4"))+#, "NODYN" = "black", "NODYN_META" = "grey50"))+#, "NORES"="black"))+
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
