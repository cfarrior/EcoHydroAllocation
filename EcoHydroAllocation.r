#Last updated 22 NOV 2014

#Rcode used to find ESS allocation strategies under different rainfall regimes. 
#and to make plots 

#to use follow these steps: 
#(1) Create two folders: "~/Documents/Farrior_etal_EcoHydroAllocation/" and "~/Documents/Farrior_etal_EcoHydroAllocation/Figures"
#(2) Save this text as a text document, named "~/Documents/Farrior_etal_EcoHydroAllocation/EcoHydroAllocation.r"
#(3) Copy and paste the following into R (without the first hashmarks)
#   setwd("~/Documents/Farrior_etal_EcoHydroAllocation/")
#	  source("EcoHydroAllocation.r")
#(4) Change any parameter values you want to, by writing their desired values into R (i.e. "t = 180" changes the growing season length to 180 days). 
#   See below and the paper for defaults.
#(5) Pick a name for the set of runs; here "runName". Running the following functions will find the ESSs and then produce the plots from the paper. 
#   FindESSs("runName")  #this will take a couple of hours; use a name to designate your run
# 	### for some reason, sometimes you have to exit the RGui and open again before the next step
#	  PaperPlots("runName") #this will produce the figures of the paper. 

#all combinations of below rainfall parameters will be run
rainV = seq(600,1800,by=100)  #total annual rainfall
lambdaV = seq(0.1,1,by=0.1)  #lambda

#Use this finer grid to get the same level of detail as the paper.   This will take much longer to run FindESSs()
#	rainV = seq(600,1600,by=25); lambdaV = seq(0.025,1,by=0.025)

#Parameter Values 

#for translation from the variables in the paper to these parameters: 
#subscript is indicated by "_"
#"," becomes "."
#zero subscripts are written "not"
#greek letters are written out and capitalized for capitals, all lowercase for lowercase.

#Table of soil moisture parameters for different soil textures
soil_text = data.frame(Soil = c("Sand","Loamy Sand", "Sandy Loam", "Loam","Clay"),b=c(4.05,4.38,4.9,5.39,11.4), Ks = c(2500,1000,800,200,50),n=c(.35,.42,.43,.45,.5),Beta=c(12.1,12.7,13.8,14.8,26.8), s_h = c(.08,.08,.14,.19,.47),s_w=c(.11,.11,.18,.24,.52),s_fc=c(.35,.52,.56,.65,1)) #data from Laio et al. 2001 Table 1 

#Select a soil texture (the code will then use the data from the table)
soil = "Loam"
Zr = 700 #mm;	rooting depth
delta = 2 #mm (trees, IRIbook 299);  interception of water by tree crowns
Ew = .1 #mm/day #Laio et al. 2001; evaporation directly from soil

#Individual tree properties - parameter sources including forest service and other sources from the literature are detailed in Appendix A Table A.1
mu_c = 0.016 #canopy tree mortality rate
mu_u = 0.038 #understory tree mortality rate
F_c = 0.0071  #seedlings/m2/year
F_u = 0  #understory trees do not reproduce (changing this parameter will only change the growth rate of understory trees.  the code has not been written to incorporate differences in population dynamics with understory reproduction)
H = 3.6 #meters/cm^(0.5) height = HD^(gama-1)
alpha_w = 0.2 # m^2/cm^(1.5) crown area = alpha_w*D^gama 
alpha_s = 48.3 # gC/cm^(2.5) structural biomass = alpha_s*D^(gamma+1)
gama = 1.5 #allometric exponent (spelled incorrectly because gamma is a soil moisture model parameter)

t = 240	#days, growing season length
L_not = 1200 #Light at the top of the canopy, MJ PAR/m2/day

V <<- Vdefault <<- 1.45 #maximum leaf level photosynthetic rate gC/m2/day 
alpha_fdefault<<-alpha_f <<- 0.00329  #gC/MJ PhAR
k = .5 #light extinction coefficient 
l_tilda = 1/k*log(alpha_f*L_not/V) #leaf layer at which leaves become light limited


#leaf costs all here
rho_l.m = 0.0915
d = 0.85
LMA_not = 9/240 #gC/m2/day 
LMA_max = 32.5/240 #gC/m2/day
rho_sw.m = 0.23 #gC/m2/day

#leaf cost calculations
cl_not = t*((LMA_max-LMA_not)*(l_tilda+1/k) + 1/k*rho_l.m*d + (rho_sw.m)*(l_tilda+1/k))
cl_lin = t*(LMA_not + rho_l.m*(1-d))
cl_exp = t*(-1/k*((LMA_max-LMA_not)*exp(k*l_tilda) + rho_l.m*d) - 1/k*(rho_sw.m)*exp(k*l_tilda))	

#given light level and the number of leaf layers, leafmass_tot computes leaf mass per unit crown area of the plant
leafmass_tot = function(L,l){
	t*((LMA_max-LMA_not)*(l_tilda+1/k) + LMA_not*l -1/k*(LMA_max-LMA_not)*exp(-k*(l-l_tilda)))

	} 
#given light level and the number of leaf layers, c_l.tot computes cost of holding l leaves per unit crown area 
c_l.tot = function(L,l) cl_not + cl_lin*l +cl_exp*exp(-k*l) #given light level and the number of leaf layers, c_l.tot computes the total yearly cost of having that LAI (l). 

RMA = 22.46 #gC/m2  fine root mass per unit area
r_l = 2 #the average lifetime of a fine-root  #work on this here!?
c_r = RMA*(1/r_l+1.2) #gC/m2/year  cost of fine root biomass per unit crown area per year             

c_b.g = 0.33 #gC/gC  growth respiration of structural biomass
	
c_f = 4870 #cost to the reproducing plant, gC/seedling 

p = 0.75 #the proportion of leaves of a canopy tree that shade the understory 

WUE <<- WUEdefault <<- 2.75  #gc/mmH20 water use efficiency 
K_p = 2.357  #mm/m2/day     conductance from soil to leaf

#the CO2 effects - enhancing alpha_f, V, and WUE.  Increase the enhancement by changing these numbers
alphaboost = 1.12
Vboost = 1.44 #this is the Vboost, if respiration doesn't change
WUEboost = 1.57  #approximate, the increase in WUE should be similar to the increase in CO2 (assuming, 350-> 550)


#other things that make calculations convenient
nums = 15
	
maxLAI = 1/k*log(t*alpha_f*alphaboost*L_not/(c_l.tot(L_not,9)-c_l.tot(L_not,8)))
maxRAI = t*V/k*(1+log(alpha_f*L_not/V) - alpha_f*L_not/V*exp(-k*maxLAI))/c_r 
bit.r = maxRAI/1000
bit.l = maxLAI/1000

rho_l = 1 #the nitrogen in the plant per LAI - not used in this paper. 

r_cfine = l_cfine = 0.01

rain.lab = expression(paste(plain("Total annual rainfall (mm "),plain(yr^-1),plain(")"),sep=""))
aV = c("(a)","(b)","(c)","(d)","(e)","(f)")

PaperPlots = function(filestemAdd="",p.cex=1,ylims=FALSE,xlims=FALSE){
	figurefilestem = paste("Figures/",filestemAdd,sep="")

	ylim=xlim=NULL

	run_ControlFert(filestemAdd)
	
	WUEonly = read.table(paste(filestemAdd,"_WUEonly.txt",sep=""),sep="\t",header=TRUE)
	WUEonly = WUEonly[order(WUEonly$rain),]
	peonly = read.table(paste(filestemAdd,"_peonly.txt",sep=""),sep="\t",header=TRUE)
	peonly = peonly[order(peonly$rain),]
	peWUE = read.table(paste(filestemAdd,"_peWUE.txt",sep=""),sep="\t",header=TRUE)
	peWUE = peWUE[order(peWUE$rain),]
	
	
	
	
	run_ControlFert(paste(filestemAdd,"_Optim",sep=""))
	
	O.WUEonly = read.table(paste(filestemAdd,"_Optim_WUEonly.txt",sep=""),sep="\t",header=TRUE)
	O.peonly = read.table(paste(filestemAdd,"_Optim_peonly.txt",sep=""),sep="\t",header=TRUE)
	O.peWUE = read.table(paste(filestemAdd,"_Optim_peWUE.txt",sep=""),sep="\t",header=TRUE)
	
	##################################################	
	{#Figure 1: ESS allocation across rainfall regimes 
	p.cex = .5
	data = WUEonly
	data = data[data$closedcanopy==1,]
	data = data[data$rain<=1600,]
	wd = 87; ht = 3
	X11(width=wd/26,height=ht)
	par(mfrow=c(2,2),mar = c(1.5,2,2,1)+0.2, oma = c(2,1,0,0),tcl=0.5,ps=10)
	
	for(i in seq(1,3)){
		if(i==1){
			y = leafmass_tot(L=L_not,l=data$l_c)
			ylab = expression(paste(plain("Leaf NPP (gC  "),plain(m^-2),plain(" "),plain(yr^-1),plain(")"),sep=""))
			ylab="Leaf NPP"
		}
		if(i==3){
			y = data$r_c*RMA/r_l
			ylab = expression(paste(plain("Fine-root NPP (gC  "),plain(m^-2),plain(" "),plain(yr^-1),plain(")"),sep=""))
			ylab="Fine-root NPP"
		}
		if(i==2){
			y = data$Anc/(1+c_b.g)
			ylab = expression(paste(plain("Wood NPP (gC  "),plain(m^-2),plain(" "),plain(yr^-1),plain(")"),sep=""))		
			ylab="Wood NPP"
		}

		xbuffer = (max(data$rain) - min(data$rain))*0.1
		ybuffer = (max(y,na.rm=TRUE) - min(y,na.rm=TRUE))*0.1

		if(i!=3) par(mfg=c(1,i)) 
		if(i==3) par(mfg=c(2,1))
		plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain)-xbuffer,max(data$rain)+xbuffer),ylim=c(min(y,na.rm=TRUE)-ybuffer,max(y,na.rm=TRUE)+ybuffer))
		mtext(ylab,line=1.6,side=2,cex=.8)
		axis(1,padj=-1)
		axis(2,padj=1)
		mtext(aV[i],adj=0.1,line=-1)	
		plot3d(y,data,p.cex=p.cex)
	}
	mtext(rain.lab,side=1,line=2,outer=FALSE,adj=0,cex=.9)
	
	par(mfg = c(2,2))
	plot3d(data = data.frame(lambda=unique(data$lambda)),legend=TRUE,legx=1.06,p.cex=p.cex)
	mtext(expression(paste(lambda,plain(" ( "),plain(day^-1),plain(")"),sep="")),line=1.5,side=2,cex=0.85)
#	mtext(aV[4],adj=-.4,line=1)
	
	dev.copy(pdf,paste(figurefilestem,"alloc_fig.pdf",sep=""),width=wd/26,height=ht)
	dev.off()
	}#end Figure 1

	##################################################	
	{#Figure 2: Carbon storage of ES communities 
	p.cex = .5
	data = WUEonly
	data = data[data$closedcanopy==1,]
	data = data[data$rain<=1600,]
	X11(width=100/26,height=2.25)
	par(mfrow=c(1,2),mar = c(0,3.5,2,.2)+0.2, oma = c(4,0,0,0),tcl=0.5,ps=10)
	
	y = data$cStore
	ylab = expression(paste(plain("Carbon storage (kgC  "),plain(m^-2),plain(" "),plain(")"),sep=""))

	xbuffer = (max(data$rain) - min(data$rain))*0.1

	plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain)-xbuffer,max(data$rain)+xbuffer),ylim=c(0,55))
	mtext(ylab,line=1.6,side=2,cex=0.8)
	axis(1,padj=-1)
	axis(2,padj=1)
	plot3d(y,data,p.cex=p.cex)
	mtext(rain.lab,side=1,line=2,outer=TRUE,adj=0.25,cex=.9)
	
	plot3d(data = data.frame(lambda=unique(data$lambda)),legend=TRUE,p.cex=p.cex)
	mtext(expression(paste(lambda,plain(" ( "),plain(day^-1),plain(")"),sep="")),line=1.5,side=2,cex=0.8)
	
	dev.copy(pdf,paste(figurefilestem,"cStore.pdf",sep=""),width=100/26,height=2.25)
	dev.off()
	}#end Figure 2
	
	##################################################	
	{#Figure 3: Carbon sinks
	p.cex = .5
	wd = 87; ht = 3
	X11(width=wd/26,height=ht)
	par(mfrow=c(2,2),mar = c(1.5,2,2,1)+0.2, oma = c(2,1,0,0),tcl=0.5,ps=10)
	for(i in seq(1,3)){
		if(i==1){
			data = WUEonly
			toplab = expression(paste(plain("enhanced "),italic(omega),sep=""))
		}
		if(i==2){
			data = peonly
			toplab = expression(paste(plain("enhanced "),italic(alpha),scriptstyle(f),plain(" and "),italic("V"),sep=""))
		}
		if(i==3){
			data = peWUE
			toplab = expression(paste(plain("enhanced "),italic(omega),", ",italic(alpha),scriptstyle(f),plain(", and "),italic("V"),sep=""))
		}
		data = data[data$closedcanopy==1,]
		data = data[data$r_cUP>0,]
		data = data[data$rain<=1600,]
		xbuffer = (max(data$rain,na.rm=TRUE) - min(data$rain,na.rm=TRUE))*0.1
		
		if(i!=3) par(mfg=c(1,i))
		if(i==3) par(mfg=c(2,1))
		data = data[data$l_cUP>1,]
		y = data$cStoreUP - data$cStore
		plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain,na.rm=TRUE)-xbuffer,max(data$rain,na.rm=TRUE)+xbuffer),ylim=c(-20,50))
		axis(1,padj=-.5)
		axis(2,padj=.5)
		if(i==1) mtext("Change in carbon storage",line=0,outer=TRUE,side=2,cex=0.9)
#		if(i==1) mtext(expression(paste(plain("(kgC "),plain(m^-2),plain(")"),sep="")),line=2,side=2,cex=0.85)
		mtext(toplab,side=3,line=.8)
		mtext(aV[i],adj=0.1,line=-1)	
		plot3d(y,data,p.cex=p.cex)
		abline(0,0)
	}
	mtext(rain.lab,side=1,line=2,outer=FALSE,adj=0,cex=.9)
	
	par(mfg = c(2,2))
	plot3d(data = data.frame(lambda=unique(data$lambda)),legend=TRUE,legx=1.06,p.cex=p.cex)
	mtext(expression(paste(lambda,plain(" ( "),plain(day^-1),plain(")"),sep="")),line=1.5,side=2,cex=0.85)
#	mtext(aV[4],adj=-.4,line=1)	
	dev.copy(pdf,paste(figurefilestem,"sinks.pdf",sep=""),width=wd/26,height=ht)
	dev.off()
	}#end Figure 3

	##################################################	
	{#Figure 4:  Carbon sink breakdown
	p.cex = .5
	data = peWUE
	data = data[data$closedcanopy==1,]
	data = data[data$rain<=1600,]
	p.cex = .5
	wd = 87; ht = 3
	X11(width=wd/26,height=ht)
	par(mfrow=c(2,2),mar = c(1.5,2,2,1)+0.2, oma = c(2,2,0,0),tcl=0.5,ps=10)
	for(i in seq(1,3)){
		if(i==1){
			NPP = leafmass_tot(L_not,data$l_c)+RMA*data$r_c/r_l+data$Anc
			NPPup = leafmass_tot(L_not,data$l_cUP)+RMA*data$r_cUP/r_l + data$Ancup
			y = NPPup/NPP
			ylab="NPP"
			xtick=1400
		}
		if(i==2){
			average.cRes = (leafmass_tot(L_not,data$l_c) + RMA*data$r_c*r_l/r_l + data$Anc/mu_c)/NPP
			average.cResUP = (leafmass_tot(L_not,data$l_cUP) + RMA*data$r_cUP*r_l/r_l + data$Ancup/mu_c)/NPPup
			y = average.cResUP/average.cRes			
			ylab="Residence time"
			xtick=1000
		}
		if(i==3){
			y = data$cStoreUP/data$cStore
			ylab="Storage"
			xtick=1400
		}
		xbuffer = (max(data$rain) - min(data$rain))*0.1
		
		if(i!=3) par(mfg=c(1,i))
		if(i==3) par(mfg=c(2,1))
		plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain)-xbuffer,max(data$rain)+xbuffer),ylim=c(0,9))
		axis(1,padj=-1)
		axis(2,padj=1)
		if(i==1) mtext("property at elevated CO2/\nproperty at ambient CO2",line=-0.5,side=2,cex=0.9,outer=TRUE)
#		text(xtick,8,ylab,adj=0,cex=.8)
		mtext(ylab,side=3,line=.8,adj=0)
		mtext(aV[i],adj=0.1,line=-1)	
		plot3d(y,data,p.cex=p.cex)
		abline(1,0)
	}
	mtext(rain.lab,side=1,line=2,outer=FALSE,adj=0,cex=.9)
	
	par(mfg = c(2,2))
	plot3d(data = data.frame(lambda=unique(data$lambda)),legend=TRUE,legx=1.06,p.cex=p.cex)
	mtext(expression(paste(lambda,plain(" ( "),plain(day^-1),plain(")"),sep="")),line=1.5,side=2,cex=0.85)
#	mtext(aV[4],adj=-.4,line=1)	
	dev.copy(pdf,paste(figurefilestem,"SinkBreakdown.pdf",sep=""),width=wd/26,height=ht)
	dev.off()
	}#end Figure 4
	
	##################################################	
	{#Figure 5.  Influence of q Figure 
	p.cex = .5
	data = peWUE
	data = data[data$closedcanopy==1,]
	data = data[data$r_cUP > 0,]
	data = data[data$G_cUP > 0,]
	data = data[data$rain<=1600,]
	X11(width=100/26,height=5)
	par(mfrow=c(2,2),mar = c(2,3.5,3,.2)+0.2, oma = c(6,0,0,0),tcl=0.5,ps=10)
	
	y = data$qUP-data$q
	ylab1 = "Change in proportion of time"
	ylab2 = expression(paste(plain("without water limitation ("),italic(q),plain(")"),sep=""))

	xbuffer = (max(data$rain,na.rm=TRUE) - min(data$rain,na.rm=TRUE))*0.1
	ybuffer = (max(y,na.rm=TRUE) - min(y,na.rm=TRUE))*0.1
	
	plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain,na.rm=TRUE)-xbuffer,max(data$rain,na.rm=TRUE)+xbuffer),ylim=c(min(y,na.rm=TRUE)-ybuffer,max(y,na.rm=TRUE)+ybuffer))
	mtext(ylab1,line=2.7,side=2,cex=0.85)
	mtext(ylab2,line=1.6,side=2,cex=0.85)
	axis(1,padj=-1)
	axis(2,padj=1)
	plot3d(y,data,p.cex=p.cex)
	mtext(aV[1],adj=0.05,line=-1)
	mtext(rain.lab,side=1,line=2,adj=0.1,cex=.9)	
		
	x = data$sink_qUP
	y = data$sink
	ylab1 = "Change in carbon storage"
	ylab2 = expression(paste(plain("(kgC "),plain(m^-2),plain(")"),sep=""))

	xbuffer = (max(x,na.rm=TRUE) - min(x,na.rm=TRUE))*0.1
	ybuffer = (max(y,na.rm=TRUE) - min(y,na.rm=TRUE))*0.1
	
	par(mfg=c(2,1))
	plot(x,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(x,na.rm=TRUE)-xbuffer,max(x,na.rm=TRUE)+xbuffer),ylim=c(min(y,na.rm=TRUE)-ybuffer,max(y,na.rm=TRUE)+ybuffer))
	mtext(ylab1,line=2.7,side=2,cex=0.85)
	mtext(ylab2,line=1.6,side=2,cex=0.85)
	axis(1,padj=-1)
	axis(2,padj=1)
	
	max = 0.9
	low = 5.5
	high = 7.5
	data$l_c = round(data$l_c)
	j = 0
	for(i in sort(unique(data$l_c,na.rm=TRUE))){
#		points(x[data$l_c==i],y[data$l_c==i],col=gray(j),pch=19,cex=p.cex)  #using this instead of the line below gives the full gradation of l_c (not rounding l_c) but it's so visually clear
		if(i>4) col=i-4
		if(i==4) col=1
		points(x[data$l_c==i],y[data$l_c==i],col=col,pch=19,cex=.8)
	}
	
	mtext("Change in proportion of time",side=1,line=2,cex=0.9)
	mtext(expression(paste("without water limitation (",italic(q),")",sep="")),side=1,line=3.3,cex=.9)
	mtext(aV[3],adj=0.05,line=-1)

	par(mfg=c(1,2))
	plot3d(data = data.frame(lambda=unique(data$lambda)),legend=TRUE)
	mtext(expression(paste(lambda,plain(" ( "),plain(day^-1),plain(")"),sep="")),line=1.5,side=2,cex=0.85)	
	mtext(aV[2],adj=-.35,line=-1)
	
	par(mfg=c(2,2))
	plot(seq(1,10),seq(1,10),type="n",xaxt="n",yaxt="n",bty="n",ylab="",xlab="")
	legend(0,9,c("5","6","7","8"),col=seq(1,4),pch=19,bty="n",cex=1.2)
	mtext(aV[4],adj=-.35,line=-1)
	mtext(expression(paste(italic(l),scriptstyle(c),plain(" ("),plain(m^2),plain(" "), plain(m^-2),plain(")"),sep="")),side=3,adj=0,line=-1.5)
	
	dev.copy(pdf,paste(figurefilestem,"deltaq.pdf",sep=""),width=100/26,height=5)
	dev.off()
	}#end Figure 5

	##################################################	
	{#Figure S1.1: lambda demo - the influence of timing of rainfall on soil moisture (demonstration of Ignacio Rodriguez-Iturbe et al. model)
	X11(width=80/25.4,height=6)
	par(mfrow=c(3,1),mar=c(.2,1.2,0,0),oma=c(5,4,4,2)+0.2,tcl=.5)
	r_c = 10
	l_c = 7
	rainV = c(600,1200,1800)
	lambdaV = seq(0.1,.9,by=0.2)
	soil = "Loam"
	iii = 0
	for(rain in rainV){
		iii = iii + 1
		for(lambda in lambdaV){
			alpha = rain/lambda/365
			a = get_pVTV(r_c,l_c,lambda,alpha,soil)
			j = (1-lambda)*0.9
			xaxt="n"
			if(iii==3) xaxt="t"
			if(lambda==0.1) plot(a$sV,a$pV,
				ylim=c(0,25),
				pch=19,col = gray(j),cex=0.7,
				xlim = c(0,1),xaxt=xaxt,type="n")
			lines(a$sV,a$pV,pch=19,col = gray(j),lwd=3)
		}
		x = 0.4
		if(rain>1500) x = 0
		y = 20
		text(x,y,paste(rain,"mm/yr"),adj=0,cex=1.2,font=1) 
	}
	#xlab
	mtext(expression(paste(plain("Soil moisture ("),italic("s"),plain(")"),sep="")),side=1,line=3)
	#ylab
	mtext(expression(paste(plain("Probability density function of soil moisture, "),italic("s"),plain("  (p("),
		italic("s"),plain("))"),sep="")),side=2,line=2,outer=TRUE)
	dev.copy(pdf,paste(figurefilestem,"lambda_demo.pdf",sep=""),width=80/25.4,height=6)
	dev.off()
	}#end Figure S1.1

	##################################################	
	{#Figure S1.2: plant demo - plant allocation can change soil moisture quite a bit
	rain = 900
	lambda = 0.35

	X11(width=169/25.4,height=3)
	par(mfrow=c(2,4),mar=c(.5,2,.25,0),oma=c(5,4,4,2),tcl=0.5)
	i = 0
	for(l_c in c(2,7)){
		for(r_c in c(2,10)){
			i = i + 1

			Ts.lab = expression(paste(plain("T("),italic("s"),plain(") (mm/day)"),sep=""))
			ps.lab = expression(paste(plain("p("),italic("s"),plain(")"),sep=""))
			s.lab = expression(paste(plain("Soil moisture ("),italic("s"),plain(")"),sep=""))

			alpha = rain/lambda/365

			paramsV = getParams(r_c,l_c,lambda,alpha,soil)
			s_star = paramsV$s_star
			a = get_pVTV(r_c,l_c,lambda,alpha,soil)
			yaxt = "l"
			if(i>1) yaxt = "n"
			par(mfg=c(1,i))
			plot(a$sV,a$TV,xlim=c(0,1),
				ylim = c(0,4),
				type="n",xaxt="n",yaxt=yaxt,
				ylab = Ts.lab
				)
			lines(a$sV,a$TV,lwd=1)
			lines(c(s_star,s_star),c(-100,100),lty=2,lwd=1)
			if(i==1) mtext(Ts.lab,side=2,line=3)
			main = paste(l_c," LAI, ",r_c," RAI",sep="")
			mtext(main,side=3,line=1,cex=.8)

			par(mfg=c(2,i))
			plot(a$sV,a$pV,ylim=c(0,14),xlim=c(0,1),
				pch=19,cex=1,type="n",yaxt=yaxt,
				ylab = ps.lab
				)
			lines(a$sV,a$pV,lwd=1)
			lines(c(s_star,s_star),c(-100,100),lty=2,lwd=1)
			if(i==1) mtext(ps.lab,side=2,line=3)
		}
	}
	mtext(s.lab,side=1,line=3,outer=TRUE)
	dev.copy(pdf,paste(figurefilestem,"plant_demo.pdf",sep=""),width=169/25.4,height=3)
	dev.off()
	}#end Figure S1.2
		
	##################################################	
	{#Figure S3.1:  Translation of rain, lambda to q and Rdry (with and without ESS allocation)
	p.cex = .5
	data = peWUE
	data = data[data$closedcanopy==1,]
	data = data[data$rain<=1600,]
	X11(width=149/26,height=5)
	par(mfrow=c(2,3),mar = c(0,5,3,.2)+0.2, oma = c(6,0,0,0),tcl=0.5,ps=10)
	
	y = data$q
	
	ylab1 = expression(italic(q))
	ylab2 = "Proportion of time without water limitation"

	xbuffer = (max(data$rain) - min(data$rain))*0.1
	ybuffer = (max(y,na.rm=TRUE) - min(y,na.rm=TRUE))*0.1
	
	par(mfg=c(1,1))	
	plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain)-xbuffer,max(data$rain)+xbuffer),ylim=c(min(y,na.rm=TRUE)-ybuffer,max(y,na.rm=TRUE)+ybuffer))
	mtext(ylab1,line=2.7,side=2,cex=0.85)
	mtext(ylab2,line=1.6,side=2,cex=0.6)
	axis(1,padj=-1)
	axis(2,padj=1)
	plot3d(y,data,p.cex=p.cex)
	mtext(aV[1],adj=0.05,line=-1.5)
	mtext(rain.lab,side=1,line=3,outer=TRUE,adj=0.3,cex=.9)	
		
	y = data$Rd

	x = data$rain
	xbuffer = (max(x,na.rm=TRUE) - min(x,na.rm=TRUE))*0.1
	ybuffer = (max(y,na.rm=TRUE) - min(y,na.rm=TRUE))*0.1
	
	
	ylab1 =expression(paste(italic(R),scriptstyle(dry),plain(" (mm "),plain(day^-1),plain(")"),sep="")) 
	ylab2 = "Mean water-limited transpiration rate"
	
	par(mfg=c(1,2))
	plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(x,na.rm=TRUE)-xbuffer,max(x,na.rm=TRUE)+xbuffer),ylim=c(min(y,na.rm=TRUE)-ybuffer,max(y,na.rm=TRUE)+ybuffer))
	mtext(ylab1,line=2.7,side=2,cex=0.85)
	mtext(ylab2,line=1.6,side=2,cex=0.6)
	axis(1,padj=-1)
	axis(2,padj=1)
	plot3d(y,data,p.cex=p.cex)
	mtext(aV[2],adj=0.05,line=-1.5)

	par(mfg=c(2,3))
	plot3d(data = data.frame(lambda=unique(data$lambda)),legend=TRUE)
	mtext(expression(paste(lambda,plain(" ( "),plain(day^-1),plain(")"),sep="")),line=1.5,side=2,cex=0.85)	

	#Rd and q when lc and rc are constant.  
	runconstlcrc = TRUE	
	if(runconstlcrc){
		l_c = 5
		r_c = 5
		constlcrc=NULL
		rainV = unique(data$rain)
		lambdaV = unique(data$lambda)
		for(rain in rainV){
			for(lambda in lambdaV){
				alpha = rain/lambda/365
				a = try(run_plants(lambda,alpha,soil,N=l_c*rho_l,r_c,l_c))
				if(is.numeric(a[[1]])){
					constlcrc = rbind(constlcrc,a)
				}
			}
		}
		constlcrc$rain = constlcrc$lambda*constlcrc$alpha*365
		
		write.table(constlcrc,"constlcrc.txt",sep="\t",row.names=FALSE)
	}
	constlcrc = read.table("constlcrc.txt",sep="\t",header=TRUE)
		

	p.cex = .5
	data = constlcrc
	data = data[data$closedcanopy==1,]
	
	y = data$q
	
	ylab1 = expression(italic(q))
	ylab2 = "Proportion of time without water limitation"

	xbuffer = (max(data$rain) - min(data$rain))*0.1
	ybuffer = (max(y,na.rm=TRUE) - min(y,na.rm=TRUE))*0.1
	
	par(mfg=c(2,1))	
	plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain)-xbuffer,max(data$rain)+xbuffer),ylim=c(min(y,na.rm=TRUE)-ybuffer,max(y,na.rm=TRUE)+ybuffer))
	mtext(ylab1,line=2.7,side=2,cex=0.85)
	mtext(ylab2,line=1.6,side=2,cex=0.6)
	axis(1,padj=-1)
	axis(2,padj=1)
	plot3d(y,data,p.cex=p.cex)
	mtext(aV[3],adj=0.05,line=-1.5)
	
	y = data$Rd
	x = data$rain
	xbuffer = (max(x,na.rm=TRUE) - min(x,na.rm=TRUE))*0.1
	ybuffer = (max(y,na.rm=TRUE) - min(y,na.rm=TRUE))*0.1
	
	ylab1 =expression(paste(italic(R),scriptstyle(dry),plain(" (mm "),plain(day^-1),plain(")"),sep="")) 
	ylab2 = "Mean water-limited transpiration rate"
	
	par(mfg=c(2,2))
	plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(x,na.rm=TRUE)-xbuffer,max(x,na.rm=TRUE)+xbuffer),ylim=c(min(y,na.rm=TRUE)-ybuffer,max(y,na.rm=TRUE)+ybuffer))
	mtext(ylab1,line=2.7,side=2,cex=0.85)
	mtext(ylab2,line=1.6,side=2,cex=0.6)
	axis(1,padj=-1)
	axis(2,padj=1)
	plot3d(y,data,p.cex=p.cex)
	mtext(aV[4],adj=0.05,line=-1.5)
	
	dev.copy(pdf,paste(figurefilestem,"translation.pdf",sep=""),width=149/26,height=5)
	dev.off()	
	}#end Figure S3.1

	##################################################	
	{#Figure S3.2:  q and q_n - points not on the line do not have a good translation from rain, lambda, to q Rdry in the last figure (not many)
	data = peWUE
	data = data[data$closedcanopy==1,]
	data = data[data$rain<=1600,]

	X11(width=3.5,height=3.5)
	par(mar=c(0.2,6,0,2)+0.2,oma=c(6,0,3,0),tcl=0.5,ps=10)
	y = data$q_n
	plot(data$q,y,pch=19,lwd=1,cex=1.5,type="n",xlab="",
		ylab="",xlim=c(0,1),ylim=c(0,1))
	plot3d(y,data,p.cex=p.cex,column=10)
	mtext(expression(italic(q)),side=1,line=3,cex=1.2)
#	mtext("Proportion of time without water limitation",side=1,line=3,cex=.8)
	mtext(expression(paste(italic(q),scriptstyle(n),sep="")),side=2,line=3,cex=1.2)
#	mtext("Proportion of time in water limitation",side=2,line=3,cex=.8)
	abline(1,-1,col="red",lwd=2)
	dev.copy(pdf,paste(figurefilestem,"qn.pdf",sep=""),width=3.5,height=3.5)
	dev.off()
	}#end Figure S3.2

	##################################################	
	{#Figure S4.1: Allocation that gives the highest community-level productivity 
	p.cex = .5
	data = O.WUEonly
	data = data[order(data$rain),]
	data = data[data$closedcanopy==1,]
	data = data[data$rain<=1600,]
	X11(width=189/26,height=2)
	par(mfrow=c(1,4),mar = c(0,3.5,2,.2)+0.2, oma = c(4,0,0,0),tcl=0.5,ps=10)
	
	for(i in seq(1,3)){
		if(i==1){
			y = leafmass_tot(L=L_not,l=data$l_c)
			ylab = expression(paste(plain("Leaf NPP (gC  "),plain(m^-2),plain(" "),plain(yr^-1),plain(")"),sep=""))
		}
		if(i==2){
			y = data$r_c*RMA/r_l
			ylab = expression(paste(plain("Fine-root NPP (gC  "),plain(m^-2),plain(" "),plain(yr^-1),plain(")"),sep=""))
		}
		if(i==3){
			y = data$Anc/(1+c_b.g)
			ylab = expression(paste(plain("Wood NPP (gC  "),plain(m^-2),plain(" "),plain(yr^-1),plain(")"),sep=""))		
		}

		xbuffer = (max(data$rain) - min(data$rain))*0.1
		ybuffer = (max(y,na.rm=TRUE) - min(y,na.rm=TRUE))*0.1

		par(mfg=c(1,i)) 	
		plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain)-xbuffer,max(data$rain)+xbuffer),ylim=c(min(y,na.rm=TRUE)-ybuffer,max(y,na.rm=TRUE+ybuffer)))
		mtext(ylab,line=1.6,side=2,cex=.8)
		axis(1,padj=-1)
		axis(2,padj=1)
		mtext(aV[i],adj=0.05,line=-1.5)	
		plot3d(y,data,p.cex=p.cex)
	}
	mtext(rain.lab,side=1,line=2,outer=TRUE,adj=0.35,cex=.9)
	
	par(mfg = c(1,4))
	plot3d(data = data.frame(lambda=unique(data$lambda)),legend=TRUE,p.cex=p.cex)
	mtext(expression(paste(lambda,plain(" ( "),plain(day^-1),plain(")"),sep="")),line=1.5,side=2,cex=0.85)
	mtext(aV[4],adj=-.4,line=-1.5)
	
	dev.copy(pdf,paste(figurefilestem,"O_alloc_fig.pdf",sep=""),width=189/26,height=2)
	dev.off()
	}#end Figure S4.1

	##################################################	
	{#Figure S4.2: Carbon storage lost because of competitive optimization
	p.cex = .5
	Optim = data.frame(rain=O.WUEonly$rain,lambda = O.WUEonly$lambda,O.cStore = O.WUEonly$cStore)
	data = merge(WUEonly,Optim)
	data = data[order(data$rain),]
	data = data[data$rain<=1600,]
	y = (data$O.cStore - data$cStore)/data$O.cStore
	ylab1="Fraction carbon storage lost"
	ylab2="to competitive fine roots"
	data = data[data$closedcanopy==1,]
	X11(width=100/26,height=2.25)
	par(mfrow=c(1,2),mar = c(0,3.5,2,.2)+0.2, oma = c(4,0,0,0),tcl=0.5,ps=10)	
	xbuffer = (max(data$rain) - min(data$rain))*0.1
	plot(data$rain,y,xlab="",ylab="",type="n",main="",xaxt="n",yaxt="n",xlim=c(min(data$rain)-xbuffer,max(data$rain)+xbuffer))
	mtext(ylab1,line=2.3,side=2,cex=0.8)
	mtext(ylab2,line=1.5,side=2,cex=0.8)
	axis(1,padj=-1)
	axis(2,padj=1)
	plot3d(y,data,p.cex=p.cex)
	mtext(rain.lab,side=1,line=2,outer=TRUE,adj=0.25,cex=.9)
	
	plot3d(data = data.frame(lambda=unique(data$lambda)),legend=TRUE,p.cex=p.cex)
	mtext(expression(paste(lambda,plain(" ( "),plain(day^-1),plain(")"),sep="")),line=1.5,side=2,cex=0.8)
	
	dev.copy(pdf,paste(figurefilestem,"cStoreloss.pdf",sep=""),width=100/26,height=2.25)
	dev.off()
	}#end Figure S4.2

	##################################################	
	{#Appendix S5 Figures: Pairwise invasibility plots
	nums = 60 
	pred = WUEonly
	pred = pred[pred$closedcanopy==1,]

	i = 0
	for(rain in c(600,1200,1600)){
		rdata = pred[pred$rain==rain,]
		for(lambda in c(.1,.4,.7,1)){
			ldata = rdata[rdata$lambda==lambda,]
			if(dim(ldata)[1]>0){
				i = i+1 
				pairwise_inv(ldata)
				dev.copy(pdf,paste(figurefilestem,"PIP",i,".pdf",sep=""),width=169/25.4,height=3)
				dev.off()
			}
		}
	}
	}#end Appendix S5

}

##################################################	
FindESSs = function(filestemAdd = "",N=maxLAI*10*rho_l){	
##################################################	
# Solves for ESS allocation under all combinations of rainfall and lambda in rainV, lambdaV
##################################################	
	
	#with competition for water
	run_rainV_lambdaV(filestem=filestemAdd,rainV=rainV,lambdaV=lambdaV,comp=TRUE,N=N)

	#without competition for water
	run_rainV_lambdaV(filestem=filestemAdd,rainV=rainV,lambdaV=lambdaV,comp=FALSE,N=N)

}

##################################################	
run_ControlFert = function(filestem,comp=TRUE){
##################################################	
# organizational function
# merges ambient and added CO2 results for specific rainfall conditions for ease of plotting
##################################################	

	control = read.table(paste(filestem,"_Control.txt",sep=""),sep="\t",header=TRUE)
  control$N = 12

	WUEup = read.table(paste(filestem,"_wueFert.txt",sep=""),sep="\t",header=TRUE)
  WUEup$N = 12
	WUEonly = mergeCF(control,WUEup,WUEboost=WUEboost,alphaboost=1,Vboost=1)
	write.table(WUEonly,paste(filestem,"_WUEonly.txt",sep=""),sep="\t",row.names=FALSE)

	peup = read.table(paste(filestem,"_peFert.txt",sep=""),sep="\t",header=TRUE)
  peup$N = 12
	peonly = mergeCF(control,peup,WUEboost=1,alphaboost=alphaboost,Vboost=Vboost)
	write.table(peonly,paste(filestem,"_peonly.txt",sep=""),sep="\t",row.names=FALSE)	

	up = read.table(paste(filestem,"_Fert.txt",sep=""),sep="\t",header=TRUE)
  up$N = 12
	peWUE = mergeCF(control,up,WUEboost=WUEboost,alphaboost=alphaboost,Vboost=Vboost)
	write.table(peWUE,paste(filestem,"_peWUE.txt",sep=""),sep="\t",row.names=FALSE)
}

##################################################	
mergeCF = function(control,fert,WUEboost=WUEboost,alphaboost=alphaboost,Vboost=Vboost,comp=TRUE,goprop=FALSE){	
##################################################	
# organizational function
# merges a control and fertilized data frame with the same rainfall conditions
##################################################	

	data1 = data.frame(rain=control$rain,lambda=control$lambda,alpha=control$alpha,soil=control$soil,N=control$N,closedcanopy=control$closedcanopy,	
		r_c = control$r_c,l_c = control$l_c,G_c=control$G_c,Anc=control$Anc,
		q=control$q,
		Rd=control$Rd,
		q_n = control$q_n)
	data2 = data.frame(rain=fert$rain,lambda=fert$lambda,alpha=fert$alpha,soil=fert$soil,N=fert$N,
		r_cUP = fert$r_c,
		l_cUP = fert$l_c,
		G_cUP = fert$G_c,
		Ancup = fert$Anc,
		qUP=fert$q,
		RdUP=fert$Rd,
		q_nUP = fert$q_n)
	data = merge(data1,data2,all.x=TRUE)	
	if(comp==TRUE){
		data = data[data$closedcanopy==1,]
		data = data[data$r_c>0,]
		data = data[data$l_c>l_tilda,]
	}

	out = NULL
	WUE <<- WUEdefault*WUEboost
	alpha_f <<- alpha_fdefault*alphaboost
	V <<- Vdefault*Vboost
	for(i in seq(1,dim(data)[1])){
		line = data[i,]
		if(!is.na(line$soil)){
			AncV = getAnc(r_c=line$r_c,l_c=line$l_c,lambda=line$lambda,alpha=line$alpha,soil=line$soil)
			Anc = AncV$Anc
			G_c = Anc/(1+c_b.g)*(alpha_w/(alpha_s*(gama+1)))	
			if(goprop){
				NPP = line$r_c*RMA/r_l + LMA*line$l_c + line$Anc
				p_r = line$r_c*RMA/r_l/NPP
				p_l = line$l_c*LMA/NPP
				AncProp = line$Anc
				l_cProp = line$l_c
				r_cProp = line$r_c
				for(jj in seq(1,100)){
					r_cProp = (AncProp + LMA*l_cProp)*p_r/(RMA/r_l*(1-p_r))
					l_cProp = (AncProp + RMA*r_cProp/r_l)*p_l/(LMA*(1-p_l))
					AncProp = max(0,getAnc(r_c = r_cProp,l_c=l_cProp,lambda=line$lambda,alpha=line$alpha,soil=line$soil)[[1]])
					if(AncProp==0) AncProp=line$Anc
				}
				NPPprop = AncProp+r_cProp*RMA + l_cProp*LMA
				print(paste(round(l_cProp*LMA/NPPprop,2),round(p_l,2)))
				print(paste(round(r_cProp*RMA/r_l/NPPprop,2),round(p_r,2)))			
			}
			if(!goprop){ l_cProp=NaN;r_cProp=NaN;AncProp=NaN }
		}
		else{ Anc=NaN; G_c = NaN;l_cProp=AncProp=r_cProp=NaN}
		out = rbind(out,data.frame(line,G_cNULL=G_c,AncNULL=Anc,qNULL=AncV$q,RdNULL=AncV$Rd,l_cProp=l_cProp,r_cProp=r_cProp,AncProp=AncProp))
	}
	WUE <<- WUEdefault
	V <<- Vdefault
	alpha_f <<-alpha_fdefault
	
	out$cStore = ((out$Anc)/mu_c + leafmass_tot(L_not,out$l_c) + out$r_c*RMA)/1000
	out$cStoreUP = ((out$Ancup)/mu_c + leafmass_tot(L_not,out$l_cUP) + out$r_cUP*RMA)/1000
	out$cStoreNULL = ((out$AncNULL)/mu_c + leafmass_tot(L_not,out$l_c) + out$r_c*RMA)/1000
	out$cStorePROP = ((out$AncProp)/mu_c + leafmass_tot(L_not,out$l_cProp) + out$r_c*RMA)/1000

	out$sinkNULL = out$cStoreNULL-out$cStore
	out$sink = out$cStoreUP-out$cStore
	out$difference = out$sink-out$sinkNULL
	out$sink_qNULL = out$qNULL - out$q
	out$sink_qUP = out$qUP - out$q
	out$sink_qDiff = out$sink_qUP - out$sink_qNULL
	out$sink_RdNULL = out$RdNULL - out$Rd
	out$sink_RdUP = out$RdUP - out$Rd
	out$sink_RdDiff = out$sink_RdNULL - out$sink_RdUP

	return(out)
}

##################################################	
run_rainV_lambdaV = function(rainV = seq(800,2000,by=100),lambdaV = seq(0.1,1,by=0.1),filestem="Loam",comp=TRUE,goup=TRUE,N=N){
##################################################	
# execute function 
# runs the ESS allocation for ambient and CO2 elevated conditions 
##################################################	

	if(comp==FALSE) filestem=paste(filestem,"_Optim",sep="")

	control_file = paste(filestem,"_Control.txt",sep="") #ambient conditions
	WUEup_file = paste(filestem,"_wueFert.txt",sep="") #elevated WUE
	peup_file = paste(filestem,"_peFert.txt",sep="") #elevated alapha_f and Vmax
	up_file = paste(filestem,"_Fert.txt",sep="") #elevated WUE, alpha_f and Vmax
	
	WUE.up <<- WUEdefault*WUEboost
	V.up <<- Vdefault*Vboost
	alpha_f.up <<- alpha_fdefault*alphaboost

	i = 0
	for(rain in rainV){
		for(lambda in lambdaV){
			alpha = rain/lambda/365
			WUE <<- WUEdefault; alpha_f <<- alpha_fdefault; V <<- Vdefault
			if(comp) a = try(run_plants(lambda,alpha,soil=soil,N),silent=TRUE)
			if(comp==FALSE) a = try(climbOptim(lambda,alpha,soil,N),silent=TRUE)
			if(is.numeric(a[[1]])){
				i = i+1
				write.table(data.frame(rain=rain,lambda=lambda,alpha=alpha,soil=soil,N=N,a),control_file,row.names=FALSE,sep="\t",col.names=i==1,append=i!=1)
				if(goup){
					#water use efficiency increases
					WUE <<- WUE.up
					if(comp) aUp = try(run_plants(lambda,alpha,soil=soil,N),silent=TRUE)
					if(comp==FALSE) aUp = try(climbOptim(lambda,alpha,soil,N),silent=TRUE)
					if(!is.numeric(aUp[[1]])) aUp = a*NaN
					write.table(data.frame(rain=rain,lambda=lambda,alpha=alpha,soil=soil,N=N,aUp),WUEup_file, row.names=FALSE,sep="\t",col.names=i==1,append=i!=1)

					#photosynthetic efficiency increases
					WUE <<- WUEdefault; alpha_f <<- alpha_f.up; V <<- V.up
					if(comp) aUp2 = try(run_plants(lambda,alpha,soil=soil,N),silent=TRUE)
					if(comp==FALSE) aUp2 = try(climbOptim(lambda,alpha,soil,N),silent=TRUE)
					if(!is.numeric(aUp2[[1]])) aUp2 = a*NaN
					write.table(data.frame(rain=rain,lambda=lambda,alpha=alpha,soil=soil,N=N,aUp2),peup_file, row.names=FALSE,sep="\t",col.names=i==1,append=i!=1)

					#all CO2 up
					WUE <<- WUE.up; alpha_f <<- alpha_f.up; V <<- V.up
					if(comp==TRUE) aUp3 = try(run_plants(lambda,alpha,soil=soil,N),silent=TRUE)
					if(comp==FALSE) aUp3 = try(climbOptim(lambda,alpha,soil,N),silent=TRUE)
					if(!is.numeric(aUp3[[1]])) aUp3 = a*NaN
					write.table(data.frame(rain=rain,lambda=lambda,alpha=alpha,soil=soil,N=N,aUp3),up_file, row.names=FALSE,sep="\t",col.names=i==1,append=i!=1)
				}
			}
			print(paste(filestem,rain,lambda))
		}
	}	
}#end run_rainV_lambdaV

##################################################	
gridOptim = function(lambda,alpha,soil,l_cV,r_cV){
##################################################	
# returns the l_c, r_c combination that gives the highest value of Photosynthesis  minus costs of leaves and roots for a monoculture
##################################################	
	AV = NULL
	for(l_c in l_cV){
		for(r_c in r_cV){
			Anc = getAnc(r_c,l_c,lambda,alpha,soil)
			AV = rbind(AV,data.frame(l_c,r_c,Anc))
		}
	}	
	optim = AV[AV$Anc==max(AV$Anc,na.rm=TRUE),]
	
	return(optim)
}#end gridOptim

##################################################	
climbOptim = function(lambda,alpha,soil,N){
##################################################	
#This function returns the optimal strategy for monocultures (i.e. if group was the level of selection)
##################################################	
	optim = gridOptim(lambda,alpha,soil,seq(0.1,min(maxLAI,N/rho_l),by=1),
		seq(.1,maxRAI,by=1))
	
	r_cOLD = optim$r_c
	l_cOLD = optim$l_c 
	AncOLD = optim$Anc
	AnV = NULL
	ttt = 0
	j = 0
	while(ttt<80){
		rand.r_c = runif(1)
		if(rand.r_c<0.3) r_cNEW = r_cOLD-r_cfine
		if(rand.r_c>0.6) r_cNEW = r_cOLD+r_cfine
		if(rand.r_c>.3) if(rand.r_c<0.6) r_cNEW = r_cOLD
		if(r_cNEW<0) r_cNEW=0
		rand.l_c = runif(1)
		if(rand.l_c<0.3) l_cNEW = l_cOLD-l_cfine
		if(rand.l_c>0.6) l_cNEW = l_cOLD+l_cfine
		if(rand.l_c>.3) if(rand.l_c<0.6) l_cNEW = l_cOLD
		if(l_cNEW<0) l_cNEW=0
		l_cNEW = min(N/rho_l,l_cNEW)
		AncV = try(getAnc(r_cNEW,l_cNEW,lambda,alpha,soil))
		if(is.numeric(AncV[[1]])){
			if(AncV$Anc>AncOLD){
				r_cOLD = r_cNEW
				l_cOLD = l_cNEW
				AncOLD = AncV$Anc
			}
			AnV = c(AnV,AncOLD)
			if(length(AnV)>1) if(AnV[length(AnV)-1]==AncOLD) ttt = ttt + 1
			else ttt = 0
			j = 0
			}
		if(!is.numeric(AncV[[1]])) j = j + 1
		if(j == 10) ttt = 80
	}
	if(j!=10){
		plot(AnV)
		AncV = getAnc(r_cOLD,l_cOLD,lambda,alpha,soil)
		G_c = AncV$Anc/(1+c_b.g)*(alpha_w/(alpha_s*(gama+1)))

		return(data.frame(l_c=l_cOLD,r_c=r_cOLD,AncV,l_c.ess=NaN,r_c.ess=NaN,closedcanopy=1,G_c=G_c))
		}		
	if(j==10) return(data.frame(l_c=NaN,r_c=NaN,Anc=NaN,q=NaN,Rd=NaN,s_star=NaN,nq=NaN,l_c.ess=NaN,r_c.ess=NaN,closedcanopy=NaN,G_c=NaN))
}#end climbOptim

##################################################	
plot3d = function(y,data,column=NaN,legend=FALSE,max=0.9,p.cex=1.5,legx=1.04){
##################################################	
# function for plotting and coding the darkness of a point as the lambda value
##################################################	

	low = 0.1
	high = 0.5

	if(!legend){
		for(i in sort(unique(data$lambda))){
			j = (1-i)*max
			x.d = data[data$lambda==i,]
			lambda.d = y[data$lambda==i]
			if(is.na(column)) points(x.d$rain,lambda.d,col=gray(j),pch=19,cex=p.cex)
			if(!is.na(column)) points(x.d[,column],lambda.d,col=gray(j),pch=19,cex=p.cex)
		}
		if(is.na(column)) lines(data[data$lambda==low,]$rain,y[data$lambda==low],col="green",lwd=2)
		if(!is.na(column)) lines(data[data$lambda==low,][,column],y[data$lambda==low],col="green",lwd=2)
		if(is.na(column)) lines(data[data$lambda==high,]$rain,y[data$lambda==high],col="blue",lwd=2)
		if(!is.na(column)) lines(data[data$lambda==high,][,column],y[data$lambda==high],col="blue",lwd=2)
	}
	if(legend){
#		y.lb = expression(paste(plain(lambda),plain(" ("),plain(day^-1),plain(")"),sep=""))
		y.lb = ""			
		plot(seq(1,dim(data)[1])*0+1.04,data$lambda,type="n",xlim=c(1,1.4),
			xaxt="n",yaxt="n",bty="n",ylim=c(0,1),xlab="",ylab=y.lb)
		for(i in sort(unique(data$lambda))){
			j = (1-i)*max
			points(legx,i,col=gray(j),pch=19,cex=p.cex)
		}
		axis(2,padj=1)
		lines(c(legx-.02,legx+.02),c(low,low),col="green",lwd=2)
		lines(c(legx-.02,legx+.02),c(high,high),col="blue",lwd=2)
	}

}#end plot3d

##################################################	
run_plants = function(lambda,alpha,soil,N,r_c=NaN,l_c=NaN){
##################################################	
# runs the numerical estimation of the ESS allocation and computes resulting stand level properties
##################################################	

	if(is.na(r_c)){
		a = newton2d(lambda=lambda,alpha=alpha,soil=soil,N=N)
		if(a==-999){
			start = get_rclc_grid(lambda,alpha,soil,N,nums=nums)
			a = try(newton2d(P=c(start[2],start[1]),lambda,alpha,soil,N))
		}
		if(is.na(a[1])) a$r_c = a$l_c = -999
		if(a[1]>10){
			start = get_rclc_grid(lambda,alpha,soil,N,nums=nums)
			a = try(newton2d(P=c(start[2],start[1]),lambda,alpha,soil,N))
		}
		if(a!=-999){
			r_c = a$r_c
			l_c = a$l_c
		}
		if(a==-999){
			r_c = 999
			l_c = 999
		}
	}	

	paramsV = getParams(r_c,l_c,lambda,alpha,soil)
	s_star = paramsV$s_star
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	Wsat = paramsV$Wsat

	if(abs(s_w-s_star)<0.01){
		AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
		CC = integrate(p.3,s_star,s_fc,paramsV)[[1]]
		DD = integrate(p.4,s_fc,1,paramsV)[[1]]
		J = 1/(AA+CC+DD)
		q = J*(CC+DD)
		q_n = 0
		Rd = NaN
	}
	
	if(Wsat){
		AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
		BB = integrate(p.2,s_w,s_star,paramsV)[[1]]
		CC = integrate(p.3,s_star,s_fc,paramsV)[[1]]
		DD = integrate(p.4,s_fc,1,paramsV)[[1]]

		J = 1/(AA + BB + CC + DD)

		q = J*(CC +DD) 
		q_n = J*BB
		Rd = J*K_p*integrate(p.2.WP,s_w,s_star,paramsV)[[1]]/q_n*r_c
	}
	if(!Wsat){
		AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
		BB = integrate(p.2,s_w,s_fc,paramsV)[[1]]
		DD = integrate(p.4,s_fc,1,paramsV)[[1]]

		J = 1/(AA + BB + DD)

		q = J*(DD) 
		q_n = J*BB
		Rd = J*K_p*integrate(p.2.WP,s_w,s_fc,paramsV)[[1]]/q_n*r_c
	}

	l_c.ess = 1/k*log((t*q*alpha_f*L_not + k*cl_exp)/cl_lin)
	r_c.ess = t*q_n*Rd*WUE/c_r
	Anc = getAnc(r_c=r_c,l_c=l_c,lambda=lambda,alpha=alpha,soil=soil)
	Anc = Anc$Anc

	G_c = Anc/(1+c_b.g)*(alpha_w/(alpha_s*(gama+1)))
	if(G_c<=0) G_c = 0

	closedcanopy = 1
	if(1 > F_c*alpha_w*G_c^gama/(mu_c^(gama+1))*gamma(gama+1)) closedcanopy = 0 
	if(G_c<=0) closedcanopy = 0

	return(data.frame(lambda,alpha,r_c,l_c,G_c,
		q,Rd,q_n,Anc,l_c.ess,r_c.ess,
		closedcanopy,s_star))
}#end run_plants 

##################################################	
newton2d = function(P=c(2,7),lambda,alpha,soil,N){
##################################################	
# newton method for finding the l_c and r_c that satisfy the ESS condition
##################################################	
	F = NULL 
	
	r_cT = NULL
	l_cT = NULL
	deltaP = 10
	thresh = 0.0005
	repyes=0
	library(MASS)
	jj = 0 
	while(sum(deltaP^2)>thresh){
		jj = jj + 1
		if(floor(jj/500)==jj/500) P = c((max(r_cT)+min(r_cT))/2,(max(l_cT)+min(l_cT))/2)
		F[1] = get_zero.rc(P[1],P[2],lambda,alpha,soil,N)
		F[2] = get_zero.lc(P[1],P[2],lambda,alpha,soil,N)
		J = matrix(nrow=2,ncol=2)
		J[1,1] = 
		(get_zero.rc(P[1]+bit.r,P[2],lambda,alpha,soil,N) - 
		 get_zero.rc(P[1]-bit.r,P[2],lambda,alpha,soil,N))/bit.r
		J[1,2] = 
		(get_zero.rc(P[1],P[2]+bit.l,lambda,alpha,soil,N) - 
		 get_zero.rc(P[1],P[2]-bit.l,lambda,alpha,soil,N))/bit.l
		J[2,1] = 
		(get_zero.lc(P[1]+bit.r,P[2],lambda,alpha,soil,N) - 
		 get_zero.lc(P[1]-bit.r,P[2],lambda,alpha,soil,N))/bit.r
		J[2,2] = 
		(get_zero.lc(P[1],P[2]+bit.l,lambda,alpha,soil,N) - 
		 get_zero.lc(P[1],P[2]-bit.l,lambda,alpha,soil,N))/bit.l
		 if(!is.na(sum(J))){ 
			 if(sum(J)!=0){
				xxx = try(solve(J),silent=TRUE)
				if(is.numeric(xxx[[1]])) deltaP = xxx%*%-F
				else return(-999)
			 }
			 if(sum(J)==0) deltaP=0
		 }
		else deltaP = (runif(2)-0.7)*50*c(bit.r,bit.l)
		P0 = P
		P = P + deltaP
		if(P[1]==P0[1]) {
				deltaP = 0
				repyes = 1
		}
		P = abs(P)
		r_cT = c(r_cT,P[1]) 
		l_cT = c(l_cT,P[2])
		if(jj>6000){
			deltaP=c(0,0)
			P = c(NaN,NaN)
			print(paste(lambda,alpha,"NoConvergence",sep=" "))
		}
	}

	r_c = P[1]
	l_c = P[2]

	return(data.frame(r_c,l_c,repyes=repyes))
}#end newton2d


##################################################	
get_rclc_grid = function(lambda,alpha,soil,N,nums=10){
##################################################	
# finds the r_c and l_c that are closest to solving the ESS condition for an grid of values. Used on a coarse grid to find the r_c and l_c to start the Newton method.
##################################################	

	grid = NULL
	zero = 100
	jj=0
	while(zero>10){
		jj = jj+1
		grid = matrix(ncol=3,nrow=(nums+1)^2)
		ii = 0
		for(r_c in seq(.1,maxRAI,by=(maxRAI-.1)/nums)){
			for(l_c in seq(2,maxLAI,by=(maxLAI-2)/nums)){
				ii = ii+1
				l_c = min(l_c,N/rho_l)
				grid[ii,] = c(l_c,r_c,abs(get_zero.lc(r_c,l_c,lambda,alpha,soil,N))+abs(get_zero.rc(r_c,l_c,lambda,alpha,soil,N)))
			}
		}
		good = grid[grid[,3]==min(grid[,3],na.rm=TRUE),]
		if(is.matrix(good)) good=good[1,]
		nums = 2*nums
		if(jj == 3) good[3]=0
		zero=good[3]
	}
	nums = min(nums,15)
	rstep = (maxRAI-.1)/nums
	lstep = (maxLAI-2)/nums
	if(good[3]==0) good[3]=10
	jj = 0
	while(good[3]>3){
		jj = jj+1
		r_cMid = good[2]
		l_cMid = good[1]
		rstepN = rstep/nums*2
		lstepN = lstep/nums*2
		grid = matrix(ncol=3,nrow=(nums+1)^2+2)
		ii = 0
		for(r_c in seq(r_cMid-rstep,r_cMid+rstep,
			by=rstepN)){
			for(l_c in seq(l_cMid-lstep,l_cMid+lstep,
				by=lstepN)){
				ii = ii + 1
				l_c = min(l_c,N/rho_l)
				grid[ii,] = c(l_c,r_c,abs(get_zero.lc(r_c,l_c,lambda,alpha,soil,N))+abs(get_zero.rc(r_c,l_c,lambda,alpha,soil,N)))
			}
		}
		good = grid[grid[,3]==min(grid[,3],na.rm=TRUE),]
		if(is.matrix(good)){
			good = good[!is.na(good[,3]),]
			if(is.matrix(good)) good = good[1,]
		}
		rstep=rstepN; lstep=lstepN
		if(jj==10) good[3]=0
	}
	return(good)
}#end get_rclc_grid

##################################################	
get_zero.rc = function(r_c,l_c,lambda,alpha,soil,N){
##################################################	
# returns the value of dLRS/dr_i evaluated at r_i = r_c for the r_c and l_c passed
# By definition the value is zero if r_c and l_c are ESS values
##################################################	

	Nlim = FALSE
	if(N/rho_l <= l_c){l_c = N/rho_l; Nlim=TRUE}
	zero.rc = NA
	paramsV = getParams(r_c,l_c,lambda,alpha,soil)
	s_star = paramsV$s_star
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	
	if(s_star > s_w+0.01 && s_star < s_fc){  
		AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
		BB = integrate(p.2,s_w,s_star,paramsV)[[1]]
		CC = integrate(p.3,s_star,s_fc,paramsV)[[1]]
		DD = integrate(p.4,s_fc,1,paramsV)[[1]]

		J = 1/(AA + BB + CC + DD)

		bug = integrate(p.2.WP,s_w,s_star,paramsV)[[1]]
		if(!Nlim) zero.rc = t*WUE*K_p*(J*bug) - c_r
		if(Nlim) zero.rc = (t*WUE*K_p*(J*bug) + t*J*(CC+DD)*alpha_f*L_not*l_c*exp(-k*l_c)/r_c - cl_lin*l_c/r_c + cl_exp*k*l_c*exp(-k*l_c)/r_c) - c_r		
		
	}

	if(s_star > s_w+0.01 && s_star > s_fc) {  
		AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
		BB = integrate(p.2,s_w,s_fc,paramsV)[[1]]
		DD = integrate(p.4,s_fc,1,paramsV)[[1]]

		J = 1/(AA + BB + DD)

		bug = integrate(p.2.WP,s_w,s_fc,paramsV)[[1]]
		zero.rc = t*WUE*K_p*(J*bug) + t*WUE*K_p*(s_fc-s_w)*DD*J - c_r
	}

	return(zero.rc)
}#end get_zero.rc

##################################################	
get_zero.lc = function(r_c,l_c,lambda,alpha,soil,N){
##################################################	
# returns the value of dLRS/dl_i evaluated at l_i = l_c for the r_c and l_c passed
# By definition the value is zero if r_c and l_c are ESS values
##################################################	
	
	if(N/rho_l <= l_c){l_c = N/rho_l; zero.lc=0}
	else{
		paramsV = getParams(r_c,l_c,lambda,alpha,soil)
		zero.lc=NA
		s_star = paramsV$s_star
		s_h = paramsV$s_h
		s_w = paramsV$s_w
		s_fc = paramsV$s_fc

		if(s_star > s_w+0.01 && s_star < s_fc){
			AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
			BB = integrate(p.2,s_w,s_star,paramsV)[[1]]
			CC = integrate(p.3,s_star,s_fc,paramsV)[[1]]
			DD = integrate(p.4,s_fc,1,paramsV)[[1]]

			J = 1/(AA + BB + CC + DD)

			if(l_c > l_tilda) zero.lc = (t*(J*CC + J*DD)*alpha_f*L_not+cl_exp*k)*exp(-k*l_c) - cl_lin
			if(l_c <= l_tilda) zero.lc = t*V*(J*CC + J*DD) + cl_exp*k*exp(-k*l_c) - cl_lin
		}

		if(s_star > s_w+0.01 && s_star > s_fc){
		#we assume if s_star > s_fc, they are just water limited all the time.  
		#there is no expression for p.4 for always water-limited plants
			if(l_c > l_tilda) zero.lc = cl_exp*k*exp(-k*l_c) - cl_lin
			if(l_c <= l_tilda) zero.lc =  cl_exp*k*exp(-k*l_c) - cl_lin
		}
	}
	return(zero.lc)
}#end get_zero.lc

##################################################	
getAnc = function(r_c,l_c,lambda,alpha,soil){
##################################################	
# returns yearly photosynthesis minus costs of leaves and fine roots 
# because in this paper we only competing individuals with varying allocation strategies; # the ess criterion using LRS reduces to a criterion for Anc  
##################################################	

	paramsV = getParams(r_c,l_c,lambda,alpha,soil)
	s_star = paramsV$s_star
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	Wsat = paramsV$Wsat

	if(Wsat){
		AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
		BB = integrate(p.2,s_w,s_star,paramsV)[[1]]
		CC = integrate(p.3,s_star,s_fc,paramsV)[[1]]
		DD = integrate(p.4,s_fc,1,paramsV)[[1]]

		J = 1/(AA + BB + CC + DD)

		Anc = t*J*(r_c*K_p*WUE*integrate(p.2.WP,s_w,s_star,paramsV)[[1]] + 
			(CC + DD)*getAL(L=L_not,l=l_c)) - c_l.tot(L=L_not,l=l_c) - c_r*r_c - c_f*F_c
		q = J*(CC + DD)
		q_n = J*BB
		Rd = J*integrate(p.2.WP,s_w,s_star,paramsV)[[1]]*r_c*K_p/q_n/t #this is the average rate of transpiration during water limitation (mm/day)		
	}
	if(!Wsat){
		AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
		BB = integrate(p.2,s_w,s_fc,paramsV)[[1]]
		DD = integrate(p.4,s_fc,1,paramsV)[[1]]

		J = 1/(AA + BB + DD)
		#assume that any s > s_fc runs off before the plants can take it up.
		Anc = t*J*(r_c*K_p*WUE*integrate(p.2.WP,s_w,s_fc,paramsV)[[1]] + 
			DD*r_c*K_p*WUE*(s_fc-s_w)) - c_l.tot(L=L_not,l=l_c) - c_r*r_c - c_f*F_c
		q = 0
		q_n = J*BB  #water limited period, when roots have some control over the water
		Rd = J*(integrate(p.2.WP,s_w,s_fc,paramsV)[[1]])*r_c*K_p/q_n/t
	}

	return(data.frame(Anc,q,Rd,s_star,q_n))
}#end getAnc

##################################################	
getAL = function(l,L){
##################################################	
# Light-limited photosynthesis, a simplified Farquhar  
##################################################	
	l_tilda = 1/k*log(alpha_f*L/V)
	if(l<=l_tilda) AL = V*l 
	if(l_tilda<l){
		if(l_tilda > 0) AL = V/k*(1+log(alpha_f*L/V) - alpha_f*L/V*exp(-k*l))
	}
	if(l_tilda<0) AL = alpha_f * L /k *(1-exp(-k*l))
	return(AL)
}#end getAL

##################################################	
pairwise_inv = function(essline){
##################################################	
# makes a pairwise invasion plot for the rainfall regime values passed through by the line from the ESS dataframe 
##################################################	

	lambda = essline$lambda
	alpha = essline$alpha
  N = essline$N 
	
	P.r_c = essline$r_c
	P.l_c = essline$l_c
  if(P.l_c*rho_l==N) pairwise_inv_Nlim(essline)
  if(P.l_c*rho_l>N){

  	r_cV = seq(0.2,P.r_c*2,by=(P.r_c*2-0.2)/nums)
  
  	l_cV = seq(0.2,P.l_c*2,by=(P.l_c*2-0.2)/nums)
  
  	X11(width=169/25.4,height=4)
  	par(mfrow=c(1,2))
  
  	r.data = NULL
  	for(r_cRES in r_cV){
  		res = NaN
  		res = try(getRes(r_cRES,P.l_c,lambda,alpha),silent=TRUE)
  		if(is.numeric(res[[1]])){
  			if(res$Anc>0){
  				for(r_cINV in r_cV){
  					inv=0
  					inv = try(getAncI(r_cINV,P.l_c,L_not,res),silent=TRUE)
  					if(is.numeric(inv[[1]])){
  						go = 0
  						if(inv$Anc>res$Anc) go = 1
  						r.data = rbind(r.data,data.frame(r_cRES,r_cINV,go))
  					}
  				}
  			}
  		}
  	}
  	godata = r.data[r.data$go==1,]
  	RR = essline$rain
  	plot(godata$r_cRES,godata$r_cINV,pch=19,
  		xlab=expression(paste(plain("Resident "),italic("r"),scriptstyle("c"),sep="")),
  		ylab=expression(paste(plain("Invader  "),italic("r"),scriptstyle("c"),sep="")),
  		main=paste("R = ",RR,", lambda = ",lambda,sep=""))
  	points(P.r_c,P.r_c,col="yellow",pch=19,cex=1)
  
  	l.data = NULL
  	for(l_cRES in l_cV){
  		res = NaN
  		res = try(getRes(P.r_c,l_cRES,lambda,alpha),silent=TRUE)
  		if(is.numeric(res[[1]])){
  			for(l_cINV in l_cV){
  				inv = 0
  				inv = try(getAncI(P.r_c,l_cINV,L_not,res),silent=TRUE)
  				go = 0
  				if(inv$Anc>res$Anc) go = 1
  				l.data = rbind(l.data,data.frame(l_cRES,l_cINV,go))
  			}
  		}
  	}
  	godata = l.data[l.data$go==1,]
  	plot(godata$l_cRES,godata$l_cINV,pch=19,
  		xlab=expression(paste(plain("Resident "),italic("l"),scriptstyle("c"),sep="")),
  		ylab=expression(paste(plain("Invader  "),italic("l"),scriptstyle("c"),sep="")))
  	points(P.l_c,P.l_c,col="yellow",pch=19,cex=1)
  }

}#end pairwise_inv

##################################################  
pairwise_inv_Nlim = function(essline){
##################################################	
# makes a pairwise invasion plot for the rainfall regime values passed through by the line from the ESS dataframe 
# if trees are N limited
##################################################	
  
  lambda = essline$lambda
  alpha = essline$alpha
  N = essline$N
  
  P.r_c = essline$r_c
  P.l_c = essline$l_c
  
  r_cV = seq(0.2,P.r_c*2,by=(P.r_c*2-0.2)/nums)
  
  X11(width=169/25.4,height=4)
  par(mfrow=c(1,2))
  
  r.data = NULL
  for(r_cRES in r_cV){
    res = NaN
    res = try(getRes(r_cRES,P.l_c,lambda,alpha),silent=TRUE)
    if(is.numeric(res[[1]])){
      if(res$Anc>0){
        for(r_cINV in r_cV){
          l_cINV = r_cINV/r_cRES*N/rho_l
          inv=0
          inv = try(getAncI(r_cINV,l_cINV,L_not,res),silent=TRUE)
          if(is.numeric(inv[[1]])){
            go = 0
            if(inv$Anc>res$Anc) go = 1
            r.data = rbind(r.data,data.frame(r_cRES,r_cINV,go))
          }
        }
      }
    }
  }
  godata = r.data[r.data$go==1,]
  RR = essline$rain
  plot(godata$r_cRES,godata$r_cINV,pch=19,
       xlab=expression(paste(plain("Resident "),italic("r"),scriptstyle("c"),sep="")),
       ylab=expression(paste(plain("Invader  "),italic("r"),scriptstyle("c"),sep="")),
       main=paste("R = ",RR,", lambda = ",lambda,sep=""))
  points(P.r_c,P.r_c,col="yellow",pch=19,cex=1)
  
}#end pairwise_inv_Nlim


##################################################	
getRes = function(r_c,l_c,lambda,alpha){
##################################################	
# returns values of the environment for a monoculture of individuals of the passed r_c, l_c 
##################################################	

	paramsV = getParams(r_c,l_c,lambda,alpha,soil)
	s_star = paramsV$s_star
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	
	AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
	BB = integrate(p.2,s_w,s_star,paramsV)[[1]]
	CC = integrate(p.3,s_star,s_fc,paramsV)[[1]]
	DD = integrate(p.4,s_fc,1,paramsV)[[1]]

	J = 1/(AA + BB + CC + DD)

	Lu = L_not*exp(-k*p*l_c)

	Anc = getAncI(r_c,l_c,L=L_not,paramsV)

	return(data.frame(Anc,paramsV))
}#end getRes

##################################################	
getAncI = function(r,l,L,paramsV){
##################################################	
# returns the value of yearly photosynthesis minus costs of allocation leaves and fine roots for the environment set by a different resident
##################################################	

	s_starR = paramsV$s_star
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	
	s_starI = getAL(L=L,l=l)/(r*K_p*WUE) + s_w 
	Wsat = TRUE
	if(s_starI > s_fc) Wsat = FALSE
		
	AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
	BB = integrate(p.2,s_w,s_starR,paramsV)[[1]]
	CC = integrate(p.3,s_starR,s_fc,paramsV)[[1]]
	DD = integrate(p.4,s_fc,1,paramsV)[[1]]
			
	if(Wsat){
		if(s_starI>=s_starR){
			AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
			BB = integrate(p.2,s_w,s_starR,paramsV)[[1]]
			CC.1 = integrate(p.3,s_starR,s_starI,paramsV)[[1]]
			CC.2 = integrate(p.3,s_starI,s_fc,paramsV)[[1]]
			DD = integrate(p.4,s_fc,1,paramsV)[[1]]
			J = 1/(AA + BB + CC.1 + CC.2 + DD)

			Anc = t*J*(r*K_p*WUE*integrate(p.2.WP,s_w,s_starR,paramsV)[[1]] + 
				r*K_p*WUE*integrate(p.3.WP,s_starR,s_starI,paramsV)[[1]] +
				(CC.2 + DD)*getAL(L=L,l=l)) - c_l.tot(L=L,l=l) - c_r*r - c_f*F_c
			if(L<L_not) Anc = Anc + c_f*F_c
			q = J*(CC.2 + DD)
			qn = J*(BB + CC.1)
			Rd = J*(integrate(p.2.WP,s_w,s_starR,paramsV)[[1]]+integrate(p.3.WP,s_starR,s_starI,paramsV)[[1]])*r*K_p/qn/t #this is the average rate of transpiration during water limitation (mm/day)
		}
		if(s_starI<s_starR){
			AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
			BB.1 = integrate(p.2,s_w,s_starI,paramsV)[[1]]
			BB.2 = integrate(p.2,s_starI,s_starR,paramsV)[[1]]
			CC = integrate(p.3,s_starR,s_fc,paramsV)[[1]]
			DD = integrate(p.4,s_fc,1,paramsV)[[1]]
			J = 1/(AA + BB.1 + BB.2 + CC + DD)

			Anc = t*J*(r*K_p*WUE*integrate(p.2.WP,s_w,s_starI,paramsV)[[1]] + 
				+ (BB.2 + CC + DD)*getAL(L=L,l=l)) - c_l.tot(L=L,l=l) - c_r*r - c_f*F_c		
			if(L<L_not) Anc = Anc + c_f*F_c
			q = J*(BB.2 + CC + DD)
			qn = J*(BB.1)
			Rd = J*(integrate(p.2.WP,s_w,s_starI,paramsV)[[1]])*r*K_p/qn/t #this is the average rate of transpiration during water limitation (mm/day)
		}
	}
	if(!Wsat){
		AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
		BB = integrate(p.2,s_w,s_fc,paramsV)[[1]]
		DD = integrate(p.4,s_fc,1,paramsV)[[1]]

		J = 1/(AA + BB + DD)
		Anc = t*J*(r*K_p*WUE*integrate(p.2.WP,s_w,s_fc,paramsV)[[1]] + 
			DD*r*K_p*WUE*(s_fc-s_w)) - c_l.tot(L=L,l=l) - c_r*r - c_f*F_c
		if(L<L_not) Anc = Anc+c_f*F_c
		q = 0
		qn = J*BB  #this is the water-limited period when roots have control over their water somewhat
		Rd = J*(integrate(p.2.WP,s_w,s_fc,paramsV)[[1]])*r*K_p/qn/t
	}

	return(data.frame(Anc,q,Rd,s_starI,qn))
}#end getAncI

##################################################	
get_pVTV = function(r_c,l_c,lambda,alpha,soil){
##################################################	
# returns the probability density and transpiration rates for the range of soil moisture values 
##################################################	

	Tmax = getAL(L=L_not,l=l_c)/WUE
	paramsV = getParams(r_c,l_c,lambda,alpha,soil)

	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_star = paramsV$s_star
	s_fc = paramsV$s_fc
	
	TV = NULL
	pV = NULL
	sV = seq(0,1,by=0.001)

	AA = integrate(p.1,s_h,s_w,paramsV)[[1]]
	BB = integrate(p.2,s_w,s_star,paramsV)[[1]]
	CC = integrate(p.3,s_star,s_fc,paramsV)[[1]]
	DD = integrate(p.4,s_fc,1,paramsV)[[1]]

	J = 1/(AA + BB + CC + DD)

	for(s in sV){
		if(s < s_h){ 
			T = 0 
			p = 0
		}
		if(s > s_h && s < s_w){
			T = 0
			p = p.1(s,paramsV)
		}
		if(s > s_w && s < s_star){
			T = Tmax*(s-s_w)/(s_star - s_w)
			p = p.2(s,paramsV)
		}
		if(s > s_star && s < s_fc){
			T = Tmax
			p = p.3(s,paramsV)
		}
		if(s > s_fc){
			T = Tmax
			p = p.4(s,paramsV)
		}
		TV = c(TV,T)
		pV = c(pV,p)
	}

	pV = J*pV

	a = data.frame(sV,TV,pV)

	return(a)
}#end get_pVTV

##################################################	
# The following functions: p.1,p.2,p.3,p.4 are the probability density of soil moisture for s in the following ranges respectively: (0,s_w);(s_w,s_star);(s_star,s_fc);(s_fc,1)
#(lacking the normalization constant)
#functions with .WP are needed for integrating over the photosynthesis at all different soil moistures. 
##################################################	
p.1 = function(s,paramsV){
	eta = paramsV$eta
	eta_w = paramsV$eta_w
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	lambda.prime = paramsV$lambda.prime
	gamma = paramsV$gamma
	beta = paramsV$beta
	s_star = paramsV$s_star
	m = paramsV$m
	1/eta_w*((s-s_h)/(s_w - s_h))^(lambda.prime*(s_w-s_h)/eta_w -1)*exp(-gamma*s)
}#end p.1

p.2 = function(s,paramsV){
	eta=paramsV$eta
	eta_w = paramsV$eta_w
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	lambda.prime=paramsV$lambda.prime
	gamma = paramsV$gamma
	beta = paramsV$beta
	s_star = paramsV$s_star
	Wsat = paramsV$Wsat
	m = paramsV$m
	if(s_star>s_fc) s_star = s_fc
	if(s_star<s_w) s_star = s_w
	1/eta_w*(1+(eta/eta_w-1)*((s-s_w)/(s_star-s_w)))^(lambda.prime*
		(s_star-s_w)/(eta-eta_w) -1)*exp(-gamma*s)
}#end p.2 

p.2.WP = function(s,paramsV){
	eta=paramsV$eta
	eta_w = paramsV$eta_w
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	lambda.prime=paramsV$lambda.prime
	gamma = paramsV$gamma
	beta = paramsV$beta
	s_star = paramsV$s_star
	Wsat = paramsV$Wsat
	m = paramsV$m
	if(s_star>s_fc) s_star = s_fc
	if(s_star<s_w) s_star = s_w
	(s-s_w)*1/eta_w*(1+(eta/eta_w-1)*((s-s_w)/(s_star-s_w)))^(lambda.prime*
		(s_star-s_w)/(eta-eta_w) -1)*exp(-gamma*s)
}#end p.2.WP

p.3 = function(s,paramsV){
	eta=paramsV$eta
	eta_w = paramsV$eta_w
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	lambda.prime=paramsV$lambda.prime
	gamma = paramsV$gamma
	beta = paramsV$beta
	s_star = paramsV$s_star
	m = paramsV$m
	1/eta*exp(-gamma*s + lambda.prime/eta*(s-s_star))*
		(eta/eta_w)^(lambda.prime*(s_star-s_w)/(eta-eta_w))
}#end p.3

p.3.WP = function(s,paramsV){
	eta=paramsV$eta
	eta_w = paramsV$eta_w
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	lambda.prime=paramsV$lambda.prime
	gamma = paramsV$gamma
	beta = paramsV$beta
	s_star = paramsV$s_star
	m = paramsV$m
	(s-s_w)*1/eta*exp(-gamma*s + lambda.prime/eta*(s-s_star))*
		(eta/eta_w)^(lambda.prime*(s_star-s_w)/(eta-eta_w))
}#end p.3.WP

p.4 = function(s,paramsV){
	eta=paramsV$eta
	eta_w = paramsV$eta_w
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	lambda.prime=paramsV$lambda.prime
	gamma = paramsV$gamma
	beta = paramsV$beta
	s_star = paramsV$s_star
	Wsat = paramsV$Wsat
	m = paramsV$m
	if(s_star>s_fc) s_star = s_fc
	if(s_star<s_w) s_star = s_w
	1/eta*exp(-(beta+gamma)*s+beta*s_fc)*(eta*exp(beta*s)/((eta-m)*
		exp(beta*s_fc) + m*exp(beta*s)))^(lambda.prime/(beta*(eta-m))+1)*
		(eta/eta_w)^(lambda.prime*(s_star-s_w)/(eta-eta_w))*
		exp(lambda.prime/eta*(s_fc-s_star))
}#end p.4

p.4.WP = function(s,paramsV){
	eta=paramsV$eta
	eta_w = paramsV$eta_w
	s_h = paramsV$s_h
	s_w = paramsV$s_w
	s_fc = paramsV$s_fc
	lambda.prime=paramsV$lambda.prime
	gamma = paramsV$gamma
	beta = paramsV$beta
	s_star = paramsV$s_star
	Wsat = paramsV$Wsat
	m = paramsV$m
	if(s_star>s_fc) s_star = s_fc
	if(s_star<s_w) s_star = s_w
	(s-s_w)*1/eta*exp(-(beta+gamma)*s+beta*s_fc)*(eta*exp(beta*s)/((eta-m)*
		exp(beta*s_fc) + m*exp(beta*s)))^(lambda.prime/(beta*(eta-m))+1)*
		(eta/eta_w)^(lambda.prime*(s_star-s_w)/(eta-eta_w))*
		exp(lambda.prime/eta*(s_fc-s_star))
}#end p.4.WP

##################################################	
getParams = function(r_c,l_c,lambda,alpha,soil){
##################################################	
# organizational function
# this returns useful parameters needed for calculation of the pdf of soil moisture  
##################################################	
	if(soil=="Sand")	text = soil_text[soil_text$Soil=="Sand",]
	if(soil=="Loamy Sand") text = soil_text[soil_text$Soil=="Loamy Sand",]
	if(soil=="Sandy Loam") text = soil_text[soil_text$Soil=="Sandy Loam",]
	if(soil=="Loam") text = soil_text[soil_text$Soil=="Loam",]
	if(soil=="Clay") text = soil_text[soil_text$Soil=="Clay",]
	
	Ks = text$Ks
	b = text$b
	n = text$n
	beta = text$Beta
	s_h = text$s_h
	s_w = text$s_w
	s_fc = text$s_fc	

	s_star = min(getAL(L=L_not,l=l_c)/(r_c*K_p*WUE) + s_w,1)
	Emax = getAL(L=L_not,l=l_c)/WUE + Ew  #mm/day
	Wsat = TRUE
	if(s_star > s_fc){
		Emax = r_c*K_p*(s_fc-s_w) + Ew
		Wsat = FALSE
	}
	eta_w = Ew/(n*Zr)
	eta = Emax/(n*Zr)
	gamma = n*Zr/alpha
	lambda.prime = lambda*exp(-delta/alpha) 
	m = Ks/(n*Zr*(exp(beta*(1-s_fc)-1))) 

	return(data.frame(eta,eta_w,s_h,s_w,s_fc,lambda.prime,gamma,beta,s_star,m,Wsat))
}#end getParams
#############################
