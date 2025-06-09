#CONTROL NO DYNAMIC - sin cambios en la dinámica de los polinizadores. Modifico #pollinators_dynamV2

#use original functions
#list2matr, 
#plant_simulation, 
#network_assemblyV2, 
#rew.rule, 

#nueva funcion en reemplazo de new_reconectar_rewiring
 
#new_reconectar_rewiring
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
			red_final= c(1, length(nuevas_abundancias), "small")
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
					#red_final="small"
					if(is.array(red_final3)==TRUE){red_final=c(dim(red_final3),"small")}else{red_final=c(length(red_final3),1,"small")}
					}
			}else{
			#print("small3")
			#red_final="small"
			if(is.array(Pmat2)==TRUE){red_final=c(dim(Pmat2),"small")}else{red_final=c(length(Pmat2),1, "small")}			
			}
	}else{
			#print("small4")
			#red_final= "small"}#cierra if Pmat
			if(is.array(Pmat)==TRUE){red_final=c(dim(Pmat),"small")}else{red_final=c(length(Pmat),1,"small")}	
			}		
	}else{
	#red_final= "small"
	#print("small5")
	if(is.array(Abmat)==TRUE){red_final=c(dim(Abmat),"small")}else{red_final=c(length(Abmat),1, "small")}
	
	}	#cierra if dimen de Abmat
				
	}
	}else{
	#red_final="small"
	#print("small8")
	if(is.array(red_anterior)==TRUE){red_final=c(dim(red_anterior),"small")}else{red_final=c(length(red_anterior),1, "small")}
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
			
			if (any(redes[[n]]=="small")){resumen_final=cbind(redes[[n]][1], redes[[n]][2],redes[[n]][3],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA, n, i)
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
write.table(tabla_final_modifiers, file = "G:\\Mi unidad\\Herbicides_pollination_Julia\\JAE\\JAE R1\\CONTROL_NO_DYNAMIC_SC3.csv", sep = ",", row.names=F)

##counting networks that became very small
small_net<-which(tabla_final_modifiers[,3]=="small")
length(small_net)
which(tabla_final_modifiers[,3]=="small")
small_net<-which(tabla_final_modifiers[,1]=="small")
length(small_net)
which(tabla_final_modifiers[,1]=="small")


