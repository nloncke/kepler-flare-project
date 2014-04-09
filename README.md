Kepler Flare Project
===============

This repo contains plans + code for tagging/vetting stellar flares in Kepler data, developed by 
Nicole Loncke (Princeton '15) & Lucianne Walkowicz (Dept of Astrophysical Sciences). 


Current task list:
- NL: vet as many flagged flares as possible (preferably b/t 200-300) to create a starter training set
- LW: think/code up some simple features to calculate on flagged flares
- LW: send NL some reading material on machine learning as applied to astronomy

3 Mar 2014: NL: first set of vetted flares due (while vetting, pay close attn to flagging pitfalls)
      	    LW: first set of features for flagged flares due

12 Mar 2014: NL + LW: first attempt at writing classifier for autovetting
             NL: fill out flare_feature_calc.py over Spring Break

12 - 31 Mar 2014: Look at results of above, refine 

31 Mar 2014: NL will run flare_feature_calc on all lightcurves with flares
	     NL will write translator to take dictionaries and put into numpy array
             LW will send scikit-learn tutorials and reading to NL
    - NL thinking of using linear support vector classification based
      on handy sklearn flowchart

9 Apr 2014: Initial trial with classifier


14 Apr 2014: Report rough draft due. Should contain a basic intro + methods section, plus any initial results from above. 
      	    Figures to include: a full Kepler lightcurve showing flares, close-up of "yes", "no" and "maybe" flagged flares,
	    	       		screen shot of flareshow vetting plot window? 

22 Apr 2014: Revised draft due. As above, but with results from classifier and initial conclusions.

6 May 2014: Final Report due to Planets & Life Program
