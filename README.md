# MSc Thesis

## Title 

Study of signal to background discrimination and Cross section measurement on production of $ZZ\rightarrow llvv$ in proton-proton collisions at 13 TeV.

## Description

This work focuses on estimating the signal and background contributions, and calculating the fiducial cross-section $σ^{fid}(pp\rightarrow ZZ\rightarrow 2l2ν)$ for the production of two Z bosons in proton-proton collisions at $\sqrt{s} = 13$ TeV, where the two Z bosons decay into 2 leptons (2l) and 2 neutrinos
(2ν). Using data collected by the ATLAS experiment at CERN, we analyze the full Run 2 sample, corresponding to an integrated luminosity of $\mathcal{L} = 140.1$ $fb^{−1}$. 

## Main Parts of the Analysis

* We define 5 datasets corresponding to 5 distinct regions. Every dataset is constructed by imposing specific selection criteria in such way that either signal events (Signal Region (SR)) are enhanced or background events (Control Regions (CRs)). All of them contain real data and Monte Carlo simulations.<n> 
Four are our primary goals:<n>
1. Optimization of SR<n>
   The SR is quite loose and full of background contamination. Thus it requires refinement. We do so by identifying such kinematic variables that are suitable for the best possible discrimination between signal and background. Subsequently, optimization is performed (using a Brute Force method and a Genetic Algorithm) on these specific kinematic variables in order to introduce further constraints on their values, always with the aim of reduc-ing the phase space while achieving the best possible signal-to-background ratio.

2. Scaling factors<n>
   Predictions from Monte Carlo simulations for background processes are used and corrected by appropriate scaling factors, such that the contribution of backgrounds explains the observed event count in the CRs. The determination of these factors is done in two ways:
   (i) Depending on how "pure" the Control Region in the background it contains is, the scaling factor of Monte Carlo events for the background is adjusted. In each subsequent CR, the calculated scaling factors are used. Finally, the corresponding scaling factor for the signal events is determined, which quantifies the compatibility of the number of measured signal events with the prediction of the Standard Model (expressed by the prediction from simula-
   tions).
   (ii) By simultaneous calculation of all scaling factors, including the signal scaling factor, by constructing the total -Log-likelihood function and finding the best combination of free parameters (scaling factors) that minimize it.

3. Validation Region (VR)<n>
   One more region is defined, VR, in order to evaluate the simultaneous application of the obtained scaling factors.

4. Cross Section Calculation<n>
   In order to translate our observations from the SR into measurements of cross section, an analysis of files corresponding to Fiducial Regions is performed. These regions are constructed in such a way as to simulate events of the specific interaction channel as they occur in nature before passing through the detector.

## Additional info about codes

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
	
## Version

This project was developed using ROOT (CERN) version 6.28/04

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

If you have any questions or feedback, feel free to contact me at [deligkasxaris@gmail.com].