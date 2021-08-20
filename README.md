# Rt_v2

#initial date: 2021/08/19
#author: Tobias Koch 
#email: r.tobiaskoch@gmail.com or toby.koch@yale.edu
#contributing author Jessica Rothman

PURPOSE: This is the second version of the Rt calculation for covidtrackerct.org
of variants. IT has been rewritten to be reproducible by any state. Case numbers pulled from
from covidestim.org infection model and frequencies from aggregated sequencing results from state labs.

version 1: pulled frequency data directly from googlesheets and wrote outputs directly to googlesheets using googlesheets API. This exists on rtobiaskoch github account
version 2: created folders in github repository for inputs and outputs for those that do not have access to connecticut data
and wish to use this data for themselves



INPUT FILES
1. estimates.csv (covid infections by state)
    -You will need to read in cases from downloaded csv from https://covidestim.org/
    -keep name "estimate.csv" and put it into your data_input folder

2. Frequency of variants table
  -Please use the example as a template for how to format.
  -Save in data_input folder as "variant_frequencies.csv"
  -The variable names only need contain the patterns "alpha","gamma", "iota", "delta", "b.1.621". This allows for more flexible variable naming
  -(FOR CT) Frequencies pulled from Grubaugh Lab Nextstrain build.
  
  
  Nextstrain frequency data source example:
 https://covidtrackerct.com/variant-surveillance/

Frequency table example:
 [variant_frequencies.csv](https://github.com/rtobiaskoch/Rt_v2/files/7017913/variant_frequencies.csv)
 
EDITING CODE
1. STATE:
     - In section COVIDESTIM DATA CLEAN in line ~88 you would need to change your state using the full capitalized name.
2. VARIANTS:
     -The variants calculated in this script are as follows: "alpha","gamma", "iota", "delta", "b.1.621".
     -If you wish to add more, you would need to add the new variants throughout the script
     
 OUTPUT
 Code will output 3 files:
 1. a plot of the Rt values
 2. the Rt_estimate.csv
 3. reformatted Rt_estimates (this is specific to the covidtrackerCt.com website)



