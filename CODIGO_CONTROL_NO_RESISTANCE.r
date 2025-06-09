##CONTROL - NO RESISTANCE PARAMETERS

###################################################################################################################
###################################################################################################################
############MODELO CONTROL#########################################################################################
###################################################################################################################
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

#use original functions
#list2matr, 
#plant_simulation, 
#network_assemblyV2, 
#rew.rule, 
#pollinators_dynamV2, 
#new_reconectar_rewiring]
#

#escenario="SCE1" #set scenario "SCE1", SCE2, SCE3
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
		#biotypes1<-as.matrix(as.data.frame(set.plants[31:38],row.names = set.plants[,1]))
		#biotypes=biotypes1/max(biotypes1)
		#herbicides characteristics and scenarios
		#her<- as.matrix(as.data.frame(herbicides[,2:9],row.names = herbicides[,1]))
		#presion<-(her["select_preasure",])/max(her["select_preasure",])
		#dosis<-her[escenario,]/max(her[escenario,])#number of applications of each herbicide (different dosses) in 2 years 
		
		#herb_risk= (rowSums(biotypes*presion*dosis))/max(rowSums(biotypes*presion*dosis))
		herb_risk = 0.5		
		#Management
		###mmodifiers de Moss et al 2019
		#if(escenario=="SCE3"){modifiers=1}#*1
		#if(escenario=="SCE2"){modifiers=0.67}#0.67
		#if(escenario=="SCE1"){modifiers=0.33}#0.33
		modifiers = 0.5
		
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
			
			if (any(redes[[n]]=="small")){resumen_final=cbind(redes[[n]][1], redes[[n]][2],redes[[n]][3],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, n, i)
			}else{

			#if (dim(redes[[n]])[2]>=2 && dim(redes[[n]])[1]>=2){#exc
				##metrics estimation
				#robustness
				extincion1<-second.extinct(redes[[n]], participant="both", method="abundance", nrep=100) #participant=low not working. try method="abundance"\"degree"?
				robustez1=robustness(extincion1)
				robustez1

				extincion2<-second.extinct(redes[[n]], participant="lower", method="abundance", nrep=100) #participant=low not working. try method="abundance"\"degree"?
				robustez2=robustness(extincion2)
				robustez2


				dimen1<-dim(redes[[n]])
				S_plantas1=dimen1[1]
				S_poli1=dimen1[2]
				links1=sum(redes[[n]])
				evenn1 = networklevel(redes[[n]],index="interaction evenness",intereven="prod")##intereven="prod" or "sum"
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

				resumen_final=cbind(S_plantas1, S_poli1,links1, evenn1, conn1,nest1,nodf,robustez1,robustez2,fun_com1,n_overHL, n_overLL, int_div, gen, vul, H2, n, i)
							
			}
			plant_pol_net_values_cycle[n,]<-resumen_final
			}
		
		plant_pol_net_values[[i]]<-plant_pol_net_values_cycle
		
	}

#create a final table containing all the results
tabla_final=do.call(rbind, plant_pol_net_values) #deshago la lista
colnames(tabla_final)<- c("S_plantas", "S_poli","links", "evenn", "conn","nest","nodf","robust", "robustPl","funct_comp", "N_overHL", "N_overLL", "int_div", "gen", "vul", "H2","Sumulation","cycle")

tabla_final_modifiers=data.frame(tabla_final)
tabla_final_modifiers
#save table
write.table(tabla_final_modifiers, file = "G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\control_no_resistance.csv", sep = ",", row.names=F)

#counting networks that became very small
small_net<-which(tabla_final_modifiers[,3]=="small")
length(small_net)
which(tabla_final_modifiers[,3]=="small")
small_net<-which(tabla_final_modifiers[,1]=="small")
length(small_net)
which(tabla_final_modifiers[,1]=="small")

