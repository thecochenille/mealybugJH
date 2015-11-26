# R analyses for Vea et al.  (submitted PloS ONE)


##Summary
This folder includes all files for the R analyses performed in the manuscript: Vea et la. Differential juvenile hormone modulation establishes extreme sexual dimorphism in scale insects (for submission in PloS ONE)

## List of files

- [sexratio.csv](https://github.com/zourloubidou/mealybugJH/blob/master/sexratio.csv): embryo counts for 5 females.
- [expressionprofile.csv](https://github.com/zourloubidou/mealybugJH/blob/master/expressionprofile.csv): data file including SDM quantity of each studied genes obtained from quantitative RT-PCR for samples of male and female mealybugs collected every 24 hour (cf. Materials and Methods, and supplementary pdf for details of experiments). List of variables as follows: 

- [Pyri5mM.csv](https://github.com/zourloubidou/mealybugJH/blob/master/Pyri5mM.csv): data file including SDM quantity of each studied gene for prepupae and pupae treated with 5mM pyriproxyfen

## List of variables for each file
###expressionprofile.csv:
- cDNA ID: unique number for cDNA
- Sample.ID: sample ID used for sample collecting
- Day.after.oviposition: number indicating the day the samples were collected for RNA extraction after they were laid.
- Stage.and.day: indicated the stage (E: embryo, L1: first instar nymph, L2: second instar nymph, L3: third instar nymph, f: female adult, m: male adult) and day within one stage.	
- Sex: male or female sample (see details in publication)

Below are the SDM values obtained from qRT-PCR

- SDM.rp49.2: SDM value for reference gene	rpL32
- SDM.JHAMT: PkJHAMT	
- SDM.Met: PkMet	
- SDM.Tai: PkTai	
- SDM.Pkkr.h1_26: PkKr-h1 core	
- SDM.Pkkr.h1A: PkKr-h1 A	
- SDM.Pkkr.h1B: PkKr-h1 B	
- SDM.Pkbr1: Pkbr1	
- SDM.Pkbr1.Z2: Pkbr1 Z2	
- SDM.Pkbr1.Z4: Pkbr1 Z4	
- SDM.Pkbr2: Pkbr2	
- SDM.Pkbr2.Z2: Pkbr2 Z2	
- SDM.Pkbr2.Z4: Pkbr2 Z4	
- SDM.Pkbr3: Pkbr3	
- SDM.Pkbr3.Z2: Pkbr3 Z2


###Pyri5mM.csv
- Sample: Sample ID
- Treatment: pyriproxyfen (5mM) or methanol	
- Stage.treated: PreD1 = prepupa (24-48 hours after molt), PuD0= pupa (0-24 hours after molt)	
- Gene	
- SDM.gene: gene SDM value from qRT-PCR	
- rp49.2: SDM rpL32 (reference gene)	

###sexratio.csv
- Mother: letter assigned to each mother
- Oviposition: day of oviposition after first oviposition (day 1)
- NumberEggs: numer of eggs oviposited
- NumberMale: number of male embryos
- NumberFemale: number of female embryos
- NumberUnknown: number of embryos that could not be sexed

##Scripts
To access sex ratio analysis click [here](https://github.com/zourloubidou/mealybugJH/blob/master/JHmealybug.md)

To access expression profiles and effect of JHM treatment click [here](https://github.com/zourloubidou/mealybugJH/blob/master/sexratio.md)






