# MSc Thesis

## Title 

Study of signal to background discrimination and Cross section measurement on production of $ZZ\rightarrow llvv$ in proton-proton collisions at 13 TeV.

## Description

This work focuses on estimating the signal and background contributions, and calculating the fiducial cross-section $σ^{fid}(pp\rightarrow ZZ\rightarrow 2l2ν)$ for the production of two Z bosons in proton-proton collisions at $\sqrt{s} = 13$ TeV, where the two Z bosons decay into 2 leptons (2l) and 2 neutrinos
(2ν). Using data collected by the ATLAS experiment at CERN, we analyze the full Run 2 sample, corresponding to an integrated luminosity of $\mathcal{L} = 140.1$ $fb^{−1}$. 

## Main Parts of the Analysis

We define 5 datasets corresponding to 5 distinct regions. Every dataset is constructed by imposing specific selection criteria in such way that either signal events (Signal Region (SR)) are enhanced or background events (Control Regions (CRs)). All of them contain real data and Monte Carlo (MC) simulations.<br> 
Four are our primary goals:<br>
1. Optimization of SR:<br>
   The SR is quite loose and full of background contamination. Thus it requires refinement. We do so by identifying such kinematic variables that are suitable for the best possible discrimination between signal and background. Subsequently, optimization is performed (using a Brute Force method and a Genetic Algorithm) on these specific kinematic variables in order to introduce further constraints on their values, always with the aim of reducing the phase space while achieving the best possible signal-to-background ratio.

2. Scaling factors:<br>
   Predictions from MC simulations for background processes are used and corrected by appropriate scaling factors, such that the contribution of backgrounds explains the observed event count in the CRs. The determination of these factors is done in two ways:<br>
   (i) Depending on how "pure" the CR in the background it contains is, the scaling factor of MC events for the background is adjusted. In each subsequent CR, the calculated scaling factors are used. Finally, the corresponding scaling factor for the signal events is determined, which quantifies the compatibility of the number of measured signal events with the prediction of the Standard Model (expressed by the prediction from simulations).<br>
   (ii) By simultaneous calculation of all scaling factors, including the signal scaling factor, by constructing the total -Log-likelihood function and finding the best combination of free parameters (scaling factors) that minimize it.

3. Validation Region (VR):<br>
   One more region is defined, VR, in order to evaluate the simultaneous application of the obtained scaling factors.

4. Cross Section Calculation:<br>
   In order to translate our observations from the SR into measurements of cross section, an analysis of files corresponding to Fiducial Regions is performed. These regions are constructed in such a way as to simulate events of the specific interaction channel as they occur in nature before passing through the detector.

## Info about Codes

* zjets_splitter.c<br>
  Splits one of the major background CR into three individual ones, regarding the different values that one variable can get (jets multiplicity: $n_{jets} = 0, 1, > 1$). This was done in order to help us analyze this background in the best possible way and use three scaling factors instead of one for this CR.
	
* Merging.C<br>
  This script merges two of the major background processes into one. This was done in order to test if the new concatenated region would provide us better results using only one scaling factor rather two for each one of them.

* ThreelCR.C, emCR_A.C, emCR_B.C, Zjets0CR.C, Zjets1CR.C, Zjets2CR.C, ZjetsCR.C<br>
  These scripts are used for the calculation of the scaling factors and the creation of histograms (Events vs variables) in order to evaluate the Data/MC ratio before and after the impementation of scaling factors.

* vol0.C<br>
  This is the main script that calculates all scaling factors for all regions (SR & CRs). It creates histograms (Events vs variables) where it is shown a comparison between signal and background events before the implementation of scaling factors. The respective distributions are presented as PDFs (Probability Density Functions). These histograms play pivotal role in the understanding of the so called well-behaved variables as they present in the second pad the statistical metric of Signal Significance utilizing the Integral of signal and background distributions. Also it provides the calculation of the cross section quantity using the events that are found in the Fiducial Regions.

* BruteForce.C<br>
  Performs an optimization routine through nested loops regarding the well-behaved variables that were defined by vol0.C.

* GA_0.C<br>
  Performs a more sophisticated optimization routine (Genetic Algorithm). The operations that are used here are: (i) Population Initialization. (ii) Fitness Function Calculation (Signal Significance) (iii) Selection (Tournament Selection || Roulette Wheel Selection) (iv) Crossover (Uniform Crossover) (v) Mutation.

* GeneticResults.C<br>
  This script provides graphs where a comparison is made between the max and average values of Signal Significance that are found during generations regarding the type of Selection operator that we use (Tournament Selection || Roulette Wheel Selection).
	
* Simultane.C<br>
  This script is used to assess all scaling factors utilizing a Simultaneous Fit approach. We constauct a total -Log-Likelihood function composed by seven individual ones (one for each region we have). The scaling factors are used as free parameters in our model and the minimization of the total -Log-Likelihood returns their values. Moreover, we compute separately the statistical and the systematic error for each one scaling factor we compute by imposing gaussian variations in the observed and MC simulations respectivelly. Finally, we create the pull distributions for the verification of the obtained statistical and systematic errors.
	
* Simultane_new.C<br>
  A more concrete approach about the Simultaneous Fit method. In this script we split each region into two in order to leave degrees of freedom in our system. It can be used only to extract all scaling factors. Concerning the statistical and systematic errors and their verification through the pull distributions Simultane.C is quite accurate.

* ScalingFactorsDisplay.C, ScalingFactorssimultaneous.C<br>
  They provide a visual representation of the obtained scaling factors by the direct estimation and the Simultaneous Fit method. They present their values together with their statistical errors. 
	
* ValidationRegion.C
  This script creates a validation region where we can apply and assess the behaviour of Data/MC ratio before and after we impose scaling factors. A compatible with unity value for the aforementioned ratio provides us confidence that we can apply scaling factors for the MC simulations in the SR.
	
## Version

This project was developed using ROOT (CERN) version 6.28/04

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

If you have any questions or feedback, feel free to contact me at [deligkasxaris@gmail.com](mailto:deligkasxaris@gmail.com).