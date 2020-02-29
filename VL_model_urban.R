##
#  R. Retkute et al. 2020
##

par.file<-"VL_par.txt"

##  Model parameters
pher.cov<-0.7  #  Pheromone coverage
no.pher.lures<-10  # Number of pheromones u.e. in each lure
new_dog_inf<- 0.13 # Proportion of newly imported dogs which were infected
pi_high<-0.37 # Proportion of infectious dogs that are highly infectious.
pi_never<-0.55 #   Proportion of exposed dogs that are never infectious.
pi_low<-1 - pi_high - pi_never  
repl_avg<-140 #  Average time (days) for deceased dog to be replaced
new_dog_inf<-0*0.13 # Probability of a newly introduced dog being exposed
nu<- 0.0055 # Per capita rate of progression of dogs from latently exposed to a further state (days^-1). 1/nu is the average duration of the latent period (days). 
delta = 0.0011; # Per capita mortality rate for dogs (days^-1).
delta_i = 0.0018; # Per capita mortality rate for dogs (days^-1). 
a_D = 0.333 #Biting rate of sand  flies (per day) (Number of times one sand  fly would want to bite a host per unit time, if hosts were freely available).
a_H = 0.125;
mu = 0.42 # sand fly mortality rate
tau=7 # latent period of sand flies
p_v_hi = 0.39 # probability of high infectiveness dogs transmitting to a sand fly
p_v_li = 0.017 # probability of low infectiveness dogs transmitting to a sand fly
p_D = 0.321 # probability of an exposed sand fly transmitting to a dog

no.inf0<-1 # Seed cases at t=0
yd=365 # Number of in a year
maxtime = 20*yd # When to stop simulations
hinftime=5*yd # When to start human cases
phertime=10*yd #When to start pheromone control
scale.flies.urban<-0.0289293 # Scaling factor in urban environment

#  Attraction parameters: scaling for D the same as for C
p1<-0.02;  p2.H<-0.0095; p3.H<-8; p2.D<-0.091; p3.D<-1; p2.C<-0.184; p3.C<-0.199; p2.P<-1; p3.P<-0.18;   p3.D<-p3.C;
par<-c(p1, p2.H, p3.H, p2.D, p3.D, p2.C, p3.C, p2.P, p3.P)
  
# Household distribution:
house<-read.table("Households_urban.txt",header=T)
# Convert coordinates into meters  
house$x<-1000*house$x; house$y<-1000*house$y;

nn<-nrow(house)
house$dogs[house$dogs==0]<-1
house$housenum<-seq(1,nn)
house$infectedhumans<-rep(0,nn)
house$pher<-rep(0,nn)
house0<-house

house.blocks<-unique(house$block)

# Calculate the distances between houses
xy <- cbind(house$x, house$y)
dist.houses<- as.matrix(dist(xy, method = 'euclidean', diag = TRUE, upper = TRUE))

H <- function(x) as.numeric(x>0)

# Sesonality in no of sandflies
SeasSc<-function(t){
	a<-270.56; b <- 0.653576; c<- 372.395
	return(a*(1+b*cos(2*pi*(t-c)/365)))
}  

Eucl.dist<-function(p1,p2){
	dist(cbind(p1,p2), method = 'euclidean', diag = TRUE, upper = TRUE) }

# Dispersal kernel
K<-function(d){
	exp(-0.16*sqrt(d))
} 

# Attraction profile:  sum over all inhabitants
A.prof<-function(d, no.H, no.D, no.C, p1, p2.H, p3.H, p2.D, p3.D, p2.C, p3.C){
	exp(-d*p1)*(p2.H*(1 - exp(-no.H*p3.H))+p2.D*(1 - exp(-no.D*p3.D))+p2.C*(1 - exp(-no.C*p3.C)))
} # 

A.host<-function(d, no, p1, p2, p3){
	exp(-d*p1)*(p2*(1 - exp(-no*p3)))
} # 

# Hight of profile times dispersion kernel
fun.1D<- function(x, d, no, par) {
	K(abs(x))* A.prof(abs(d-x), no[1], no[2], no[3], par[1], par[2], par[3], par[4], par[5], par[6], par[7])
}  

fun.1D.hst<- function(x, d, no, par) {
	K(abs(x))* A.host(abs(d-x), no, par[1], par[2], par[3])
}  

# Integrate along x
int.fun.1D<-function(no, par){
	d<-300; dt<-0.01
	DT<-seq(0,d,dt)
	if(length(no)==1) {
		sm<-sum(fun.1D.hst(DT, d, no, par))
	} else {
		sm<-sum(fun.1D(DT, d, no, par))
	}
	return(dt*sm)
}

calc.frac.fls<-function(i, par){
	no<-as.numeric(house[i,5:7])
	fr<-int.fun.1D(no, par)
	return(fr)
}

calc.frac.fls.house<-function(i){
	no<-as.numeric(house[i,5])
	fr<-int.fun.1D(no, c(p1, p2.H, p3.H))
	return(fr)
}

calc.frac.fls.outside<-function(i){
	no<-as.numeric(house[i,6:7])
	fr<-int.fun.1D(no[1], c(p1, p2.D, p3.D))+ int.fun.1D(no[2], c(p1, p2.C, p3.C))
	return(fr)
} 

calc.frac.fls.pher<-function(i, t){
	no<-as.numeric(house$pher[i])
	if(t>= phertime){
		fr<-int.fun.1D(no, c(p1, p2.P, p3.P))
	} else {
		fr<-0
	}
	return(fr)
}

get.norm.frac.fl<-function(t){
	frc.fl.H<-rep(0,nn)
	for(i in 1:nn){ frc.fl.H[i]<-calc.frac.fls.house(i)}
	frc.fl.O<-rep(0,nn)
	for(i in 1:nn){ frc.fl.O[i]<-calc.frac.fls.outside(i) }
	frc.fl.P<-rep(0,nn)
	for(i in 1:nn){ frc.fl.P[i]<-calc.frac.fls.pher(i, t)}
	a<-frc.fl.H+frc.fl.O
	frc.fl.H<-a*frc.fl.H/(a+frc.fl.P)
	frc.fl.O<-a*frc.fl.O/(a+frc.fl.P)
	frc.fl.H<-frc.fl.H/sum(a)
	frc.fl.O<-frc.fl.O/sum(a)
	frc.fl.P<-frc.fl.P
 return(list(frc.fl.H, frc.fl.O, frc.fl.P))
}

test_frac_flies<- function(){
	frc.fl<-rep(0,nn)
	for(i in 1:nn){
		frc.fl[i]<-calc.frac.fls(i,par)
	}
	frc.fl.H<-rep(0,nn)
	for(i in 1:nn){
		frc.fl.H[i]<-calc.frac.fls.house(i)
	}
	frc.fl.O<-rep(0,nn)
	for(i in 1:nn){
		frc.fl.O[i]<-calc.frac.fls.outside(i)
	}
	hist(frc.fl.H+frc.fl.O-frc.fl)
}
	
getDogs<-function(){
	dogs<-data.frame(IDs=c(),inf_stat=c(),house=c(), introduced=c(), exposed=c(),infectious=c(), removed=c(), replaced=c(), block=c())
	id<-1
	for(i in 1:nn){ 
		if(house$dogs[i]>0){
			for(j in 1:house$dogs[i]) {	 
				infstat<-sample(c(1,2,3),1,prob=c(pi_never, pi_low, pi_high))
				   dogs<- data.frame(IDs=c(dogs$IDs, id),inf_stat=c(dogs$inf_stat,infstat), house=c(dogs$house, house$housenum[i]), introduced=c(dogs$introduced, 0), exposed=c(dogs$exposed, -1), infectious=c(dogs$infectious, -1), removed=c(dogs$removed, -1), replaced=c(dogs$replaced,0) , block=c(dogs$block, house$block[i]))
				id<-id+1  }}}  
	return(dogs)
}

#  Run model until equilibrium
run_model_br_equilibrium<-function(scal_fl){
	dogs<-getDogs()
	for (i in 1:length(house.blocks)){
		inf0<-sample(dogs$ID[dogs$inf_stat>1 & dogs$block==house.blocks[i]], no.inf0)
		dogs$exposed[inf0]<-1
		dogs$infectious[inf0]<-1
	}  
	inf.bites.H<-c()
	inf.bites.D<-c()
	inc.rate.D<-c()
	inf.bites.CH<-matrix(0,nrow=nn, ncol=maxtime-hinftime)
	FR<-get.norm.frac.fl(1)
	frc.fl.H<-FR[[1]]; frc.fl.O<-FR[[2]]; 
	frc.fl<- frc.fl.H+frc.fl.O	
	sim.D<-c()
	K.mn<- rep(0,nn);
	for (t in 2:(phertime-1)){ 
		if(t>2) K.mn<-mean(p_v_hi* K.HI  +  p_v_li* K.LI)
		new.inf.D<-0
		total.fl<- scal_fl*nn*SeasSc(t)  
        prevt.D<-100*(length(dogs$ID[dogs$infectious>0 & dogs$removed<0]))/length(dogs$ID[dogs$removed<0])
        # S -> E 
        infHI<-which(dogs$inf_stat==3 & dogs$infectious>0 & dogs$infectious<t & dogs$removed<0)
        infLI<-which(dogs$inf_stat==2 & dogs$infectious>0 & dogs$infectious<t & dogs$removed<0)
        	K.HI<- rep(0,nn); K.LI<- rep(0,nn);
        if(length(infHI)+length(infLI)>0){ 
        	if(length(infHI)==1)
        		K.HI<-as.numeric(K1[,dogs$house[infHI]])
        	if(length(infHI)>1)
        		K.HI<-as.numeric(rowSums(K1[,dogs$house[infHI]]) )
        	if(length(infLI)==1)
        		K.LI<-K1[,dogs$house[infLI]]
        	if(length(infLI)>1)
        		K.LI<-as.numeric(rowSums(K1[,dogs$house[infLI]]) )
            # vectorial capacity @ each house sum(frc.fl)
       		VC_D<- (total.fl*frc.fl.O/pmax(house$dogs, 1))*a_D^2 *exp(-mu *tau)/mu
			lambda <- H(house$dogs)*((p_D * VC_D)/pmax(house$dogs, 1))*(p_v_hi* K.HI  +  p_v_li* K.LI + K.NB* K.mn)
			pht<-(1-exp(-lambda)) # plot(seq(0,1, by=0.01),(1-exp(-seq(0,1, by=0.01))), type='l')
			u <- runif(nrow(dogs))   # head(dogs)
			no_exp<-which(u<pht[dogs$house] & dogs$exposed<0 & dogs$removed<0)
			if(length(no_exp)>0){
				dogs$exposed[no_exp]<-t
			}
		}   
		# E -> I
		wh<-which(dogs$exposed>0 & dogs$exposed<t & dogs$infectious<0 & dogs$removed<0)
		if(length(wh)>0){
			inf_prob<-1-exp(-rep(nu,length(wh)))
			u <- runif(length(wh))   
			no_inf<-which(u<inf_prob)
			new.inf.D<-new.inf.D+length(no_inf)
			if(length(no_inf)>0){
				dogs$infectious[wh[no_inf]]<-t
			}
		}
		# Remove susceptable/exposed dogs 
		whS<-which(dogs$removed<0 & dogs$infectious<0)
		if(length(whS)>0){
			rem_prob<-1-exp(-rep(delta,length(whS)))
			u <- runif(length(whS)) 
			no_rem<-which(u<rem_prob)
			if(length(no_rem)>0){
				dogs$removed[whS[no_rem]]<-t
				house$dogs[dogs$house[whS[no_rem]]]<-house$dogs[dogs$house[whS[no_rem]]]-1
			}
		}
		# Remove infectious dogs 
		whI<-which(dogs$removed<0 & dogs$infectious>-1)
		if(length(whI)>0){
			rem_prob<-1-exp(-rep(delta_i,length(whI)))
			u <- runif(length(whI)) 
			no_rem<-which(u<rem_prob)
			if(length(no_rem)>0){
				dogs$removed[whI[no_rem]]<-t
				house$dogs[dogs$house[whI[no_rem]]]<-house$dogs[dogs$house[whI[no_rem]]]-1
			}
		}
		if(length(whS)+ length(whI)>0){
			FR<-get.norm.frac.fl(t)
			frc.fl.H<-FR[[1]]; frc.fl.O<-FR[[2]]; frc.fl<- frc.fl.H+frc.fl.O
		}
		# Introduce new dogs 
		wh<-which(dogs$replaced==0 & dogs$removed >0 & dogs$removed + repl_avg ==t)
		if (length(wh)>0){
			for(i in 1:length(wh)){
				id<-max(dogs$IDs)+1
				dogs$replaced[wh[i]]<-id
				infstat<-sample(c(1,2,3),1,prob=c(pi_never, pi_low, pi_high))
				latper<-rpois(1,1/nu)
				infectiondate<- sample(c(t,-1), 1, prob=c(new_dog_inf, 1-new_dog_inf))
				dogs<- data.frame(IDs=c(dogs$IDs, id),inf_stat=c(dogs$inf_stat,infstat), house=c(dogs$house, dogs$house[wh[i]]), introduced=c(dogs$introduced, t), exposed=c(dogs$exposed, infectiondate), infectious=c(dogs$infectious, infectiondate), removed=c(dogs$removed, -1), replaced=c(dogs$replaced,0))
			}
			house$dogs[dogs$house[wh]]<-house$dogs[dogs$house[wh]]+1
			FR<-get.norm.frac.fl(t)
			frc.fl.H<-FR[[1]]; frc.fl.O<-FR[[2]]; frc.fl<- frc.fl.H+frc.fl.O
		}
	# track infected bites
	 inf.bites.D<-c(inf.bites.D, sum(H(house$dogs)*(total.fl*frc.fl.O/pmax(house$dogs, 1))*a_D^2* (K.HI  +  K.LI) * exp(-mu *tau)/mu))
	# humans 
	if(t>hinftime){
		inf.bites.CH[,t-hinftime]<-H(house$humans)*(total.fl*frc.fl.H/pmax(house$humans, 1))*a_H*a_D* (K.HI  +  K.LI + K.NB* K.mn) * exp(-mu *tau)/mu
		
	}
	prevt.D<-100*(length(dogs$ID[dogs$infectious>0 & dogs$removed<0]))/length(dogs$ID[dogs$removed<0])
	sim.D<-c(sim.D, prevt.D);
	inc.rate.D<-c(inc.rate.D, new.inf.D)
 } # t
	return(list(sim.D, inc.rate.D, inf.bites.CH, inf.bites.H, inf.bites.D, house, dogs, FR, K.mn))
}

#  Run model with control
run_model_br_pher<-function(scal_fl, simulEq, housePh){
	sim.D<-simulEq[[1]]; inc.rate.D<-simulEq[[2]]; inf.bites.CH<-simulEq[[3]]; inf.bites.H<-simulEq[[4]]; inf.bites.D<-simulEq[[5]]; house<-housePh; dogs<-simulEq[[7]]; FR<-simulEq[[8]]; K.mn<-simulEq[[9]]
	frc.fl.H<-FR[[1]]; frc.fl.O<-FR[[2]]; 
	frc.fl<- frc.fl.H+frc.fl.O	
	for (t in phertime:maxtime){ 
		if(t>phertime) K.mn<-mean(p_v_hi* K.HI  +  p_v_li* K.LI)
		new.inf.D<-0
		# is it time for pheromones?
		if(t==phertime){
			frc.fl.H<-FR[[1]]; frc.fl.O<-FR[[2]]; 
			frc.fl<- frc.fl.H+frc.fl.O
		}		
		
		total.fl<- scal_fl*nn*SeasSc(t)  
		
        prevt.D<-100*(length(dogs$ID[dogs$infectious>0 & dogs$removed<0]))/length(dogs$ID[dogs$removed<0])
        # S -> E 
        infHI<-which(dogs$inf_stat==3 & dogs$infectious>0 & dogs$infectious<t & dogs$removed<0)
        infLI<-which(dogs$inf_stat==2 & dogs$infectious>0 & dogs$infectious<t & dogs$removed<0)
        	K.HI<- rep(0,nn); K.LI<- rep(0,nn);
        if(length(infHI)+length(infLI)>0){ 
        	if(length(infHI)==1)
        		K.HI<-as.numeric(K1[,dogs$house[infHI]])
        	if(length(infHI)>1)
        		K.HI<-as.numeric(rowSums(K1[,dogs$house[infHI]]) )
        	if(length(infLI)==1)
        		K.LI<-K1[,dogs$house[infLI]]
        	if(length(infLI)>1)
        		K.LI<-as.numeric(rowSums(K1[,dogs$house[infLI]]) )
            # vectorial capacity 
       		VC_D<- (total.fl*frc.fl.O/pmax(house$dogs, 1))*a_D^2 *exp(-mu *tau)/mu
			lambda <- H(house$dogs)*((p_D * VC_D)/pmax(house$dogs, 1))*(p_v_hi* K.HI  +  p_v_li* K.LI + K.NB* K.mn)
			pht<-(1-exp(-lambda)) 
			u <- runif(nrow(dogs))   
			no_exp<-which(u<pht[dogs$house] & dogs$exposed<0 & dogs$removed<0)
			if(length(no_exp)>0){
				dogs$exposed[no_exp]<-t
			}
		}   
		# E -> I
		wh<-which(dogs$exposed>0 & dogs$exposed<t & dogs$infectious<0 & dogs$removed<0)
		if(length(wh)>0){
			inf_prob<-1-exp(-rep(nu,length(wh)))
			u <- runif(length(wh))   
			no_inf<-which(u<inf_prob)
			new.inf.D<-new.inf.D+length(no_inf)
			if(length(no_inf)>0){
				dogs$infectious[wh[no_inf]]<-t
			}
		}
		# Remove susceptable/exposed dogs 
		whS<-which(dogs$removed<0 & dogs$infectious<0)
		if(length(whS)>0){
			rem_prob<-1-exp(-rep(delta,length(whS)))
			u <- runif(length(whS)) 
			no_rem<-which(u<rem_prob)
			if(length(no_rem)>0){
				dogs$removed[whS[no_rem]]<-t
				house$dogs[dogs$house[whS[no_rem]]]<-house$dogs[dogs$house[whS[no_rem]]]-1
			}
		}
		# Remove infectious dogs 
		whI<-which(dogs$removed<0 & dogs$infectious>-1)
		if(length(whI)>0){
			rem_prob<-1-exp(-rep(delta_i,length(whI)))
			u <- runif(length(whI)) 
			no_rem<-which(u<rem_prob)
			if(length(no_rem)>0){
				dogs$removed[whI[no_rem]]<-t
				house$dogs[dogs$house[whI[no_rem]]]<-house$dogs[dogs$house[whI[no_rem]]]-1
			}
		}
		if(length(whS)+ length(whI)>0){
			FR<-get.norm.frac.fl(t)
			frc.fl.H<-FR[[1]]; frc.fl.O<-FR[[2]]; frc.fl<- frc.fl.H+frc.fl.O
		}
		
		# Introduce new dogs 
		wh<-which(dogs$replaced==0 & dogs$removed >0 & dogs$removed + repl_avg ==t)
		if (length(wh)>0){
			for(i in 1:length(wh)){
				id<-max(dogs$IDs)+1
				dogs$replaced[wh[i]]<-id
				infstat<-sample(c(1,2,3),1,prob=c(pi_never, pi_low, pi_high))
				latper<-rpois(1,1/nu)
				infectiondate<- sample(c(t,-1), 1, prob=c(new_dog_inf, 1-new_dog_inf))
				dogs<- data.frame(IDs=c(dogs$IDs, id),inf_stat=c(dogs$inf_stat,infstat), house=c(dogs$house, dogs$house[wh[i]]), introduced=c(dogs$introduced, t), exposed=c(dogs$exposed, infectiondate), infectious=c(dogs$infectious, infectiondate), removed=c(dogs$removed, -1), replaced=c(dogs$replaced,0))
			}
			house$dogs[dogs$house[wh]]<-house$dogs[dogs$house[wh]]+1
			FR<-get.norm.frac.fl(t)
			frc.fl.H<-FR[[1]]; frc.fl.O<-FR[[2]]; frc.fl<- frc.fl.H+frc.fl.O
		}
	# track infected bites
	 inf.bites.D<-c(inf.bites.D, sum(H(house$dogs)*(total.fl*frc.fl.O/pmax(house$dogs, 1))*a_D^2* (K.HI  +  K.LI) * exp(-mu *tau)/mu))
	# humans 
	if(t>hinftime){
		inf.bites.H<-c(inf.bites.H, sum(H(house$humans)*(total.fl*frc.fl.H/pmax(house$humans, 1))*a_H*a_D* (K.HI  +  K.LI) * exp(-mu *tau)/mu))
	}
	prevt.D<-100*(length(dogs$ID[dogs$infectious>0 & dogs$removed<0]))/length(dogs$ID[dogs$removed<0])
	sim.D<-c(sim.D, prevt.D);

 } # t
	return(list(inf.bites.D, inf.bites.H, prevt.D, sim.D))
}


##  Run simulation
	
	K.NB<-K(house0$dist) 
	K1<-K(dist.houses) 
	#  Reset to initial configuration
	house<-house0
	par.all<-scan(par.file)
	p<-sample(1:length(par.all),1)
	scal_fl<-scale.flies.urban*par.all[p]
	#  Run to equilibrium
	simulEq<-run_model_br_equilibrium(scal_fl) 
	
	#  Control
	add.pher<-sample(1:nn, floor(pher.cov*nn))
	house<-simulEq[[6]]
	house$pher[add.pher]<-no.pher.lures
	simul<-run_model_br_pher(scal_fl, simulEq, house)  
	
	# Visualise results
	plot(seq(1:length(simul[[4]]))/yd, simul[[4]], xlab="Time", ylab="VL prevalence in dogs", type="l", ylim=c(0,100))
			