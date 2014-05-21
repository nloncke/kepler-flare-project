Kepler Flare Project
====================

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
      - fix bugs in formatter.............DONE 
      - adjust second derivative metric...


14 Apr 2014: Report rough draft due. 
       	     Methods section can be test from workflow doc, include
       	     description of flare feature calculations and motivation
       	     (quantifying how the human eye/mind decides something is
       	     a flare), plus any initial results from above.
      	    
	    Figures to include: a full Kepler lightcurve showing flares, close-up of "yes", "no" and "maybe" flagged flares, screen shot of flareshow vetting plot window? 

21 Apr 2014: Run through iris classifier with real flare data: compare
   accuracy against test set (split 200 labelled data in half).  Write
   additional formatter to add labels to the data to feed into
   classifier.  Update often!

22 Apr 2014: Revised draft due. As above, but with results from classifier and initial conclusions.

6 May 2014: Final Report due to Planets & Life Program

Next steps: store output of classifier (even as a plain text file)
            run flareint on everything classifier says is a flare.
            plot the data
	    
	    change learning curve accuracy to account for the yeses
	    also, do subset selection as per L's email
	   
RETRIEVE:
 1. classifier output
 2. for each flare, the stats

WORKFLOW:
 1. Run flareFeatures() on all the mdwarf files.  You need to pass a
     list of the filenames and a list of their corresponding flag files.
 2. Call run_classifier() on the output of flareFeatures().  This will
     separate the labelled from the testing data, train the classifier,
      and make predictions on the test data.
