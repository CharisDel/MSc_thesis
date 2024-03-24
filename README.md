# MSc Thesis

## title: Study of signal to background discrimination and Cross section measurement on production of $ZZ\rightaroow llvv$ in proton-proton collisions at 13 TeV

## Description

This work focuses on estimating the signal and background contributions, and calculating the fiducial cross-section $σ^{fid}(pp\rightarrow ZZ\rightarrow 2l2ν)$ for the production of two Z bosons in proton-proton collisions at $\sqrt{s} = 13$ TeV. Using data collected by the ATLAS experiment at CERN, we analyze the full Run 2 sample with an integrated luminosity of 140.1 $fb^{−1}$. 

* To basiko arxeio einai to vol0.C (kakos xamos)	
	Edw upologizei olous tous scaling factors (direct etimation (per Region)) enan enan analoga me to 
	purity tis perioxis (top(emCR_B)->WZ(3lCR)->Zjets2->Zjets1->Zjets0->WW(emCR_A)) kai telos gia tin SR
	to signal strength μ_s
	Telos kanw kati sximata gia tin SR
	
* Brute Force Optimization (BruteForce.C)

* Genetic Algorithm Optimization (GA_0.C)

* Genetic Algorithm Graphs (GeneticResults.C)
	sximata apla gia tin max signiicance kai tin Average significance gia Tournament selection kai 
	RouletteWheel selection operators 
	
* Simultaneous Fit (Simultane.C)
	Gia to simultaneous "fit" xrisimopoiwntas 7 mono likelihood functions (dinei kala apotelesmata)
	Scaling factors + statistiko sfalma + sustimatiko sfalma + pull distributions gia stat kai syst 
	
* Simultaneous Fit (Simultane_new.C)
	Kanonika gia na kaneis fit prepei na afiseis bathmous eleutherias. Arithmos eksiswsewn - arithmo 
	parametrwn = Degrees of freedom. Gi auto spame tin individual likeliood tis kathe region se 2 
	epimerous-> 14 eksiswseis - 7 parameters eisai OK
	Paraskeui to doulepsa kai Deutera parousiaza. Den to oloklirwsa...
	
* Validation Regions (ValidationRegion.C)
	Ta cuts pou bazw ta exw sto keimeno an ginei kai gia tis duo
	Mono gia mia validation region (mikri lwrida) einai i Validation Region B ousiastika (keimeno)
	
* (ScalingFactorsDisplay.C)
	Pws ebgala ekino to sximataki apla gia na parousiasw tous scaling factors
	
* (Zjets_splitter.C)
	Pws spaei to directory-samples Zjets se epimerous analoga me ton jet multiplicity
	
Gia opoiadipote aporia as epikoinwniseis mazi mou!!! ( Δελιγκάς Χάρης :)
