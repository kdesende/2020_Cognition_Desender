# 2020_Cognition_Desender
Repo for "Dynamic expressions of confidence within an evidence accumulation framework", authored by Kobe Desender, Tobias H Donner, &amp; Tom Verguts. Published in Cognition.

#Raw data

The raw data from this paper can be found in the file rawData.txt

#Behavioral analysis and model predictions

The file behavioral_analysis.r contains the code needed to do all  preprocessing of the raw data, as well as generating model predictions from the hddm fit and analyzing both empirical data and model predictions. 

#HDDM fitting 

Using the datafile generated in the behavioral file above, the file "fit_data_with_hddm.py" can be used to fit these data using the HDDM module (v0.6.0) in python (v2.7).
This file saves the estimates and distributions which in turn are loaded and uses in the behavioral file above

#Help

If anything is unclear, be sure to get in contact at kobe (dot) desender (at) gmail (dot) com
