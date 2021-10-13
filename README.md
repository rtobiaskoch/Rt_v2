Rt_v2

initial date: 2021/08/19
author: Tobias Koch 
email: r.tobiaskoch@gmail.com or toby.koch@yale.edu
contributing author: Jessica Rothman
email: jessica.rothman@yale.edu

PURPOSE: This is the second version of the Rt calculation for covidtrackerct.org
of variants. IT has been rewritten to be reproducible by any state. Case numbers pulled from
from covidestim.org infection model and frequencies from aggregated sequencing results from state labs.

version 1: pulled frequency data directly from googlesheets and wrote outputs directly to googlesheets using googlesheets API. This exists on rtobiaskoch github account
version 2: created folders in github repository for inputs and outputs for those that do not have access to connecticut data
and wish to use this data for themselves

GETTING STARTED:
1. Clone the repository git clone https://github.com/rtobiaskoch/Rt_v2  be sure to change your working directory to the Rt_v2 folder

INPUT FILES
1. estimates.csv (covid infections by state)
    -You will need to read in cases from downloaded csv from https://covidestim.org/
    -keep name "estimate.csv" and put it into your data_input folder

2. Frequency of variants table
  -Please use the example as a template for how to format.
  -Save in data_input folder as "variant_frequencies.csv"
  -The variable names only need contain the patterns "alpha","gamma", "iota", "delta", "b.1.621". This allows for more flexible variable naming (Note if other variants are added the entire code will need to be edited.)
  -(FOR CT) Frequencies pulled from Grubaugh Lab Nextstrain build. If you wish to use this repository you will need to collect data in the same format
  
  
  Nextstrain frequency data source example:
 https://covidtrackerct.com/variant-surveillance/

Frequency table example:
 [variant_frequencies.csv](https://github.com/rtobiaskoch/Rt_v2/files/7017913/variant_frequencies.csv)
 
EDITING CODE: Rt Estimate Ct Infections.R
1. STATE:
     - In section COVIDESTIM DATA CLEAN in line ~88 you would need to change your state using the full capitalized name.
2. VARIANTS:
     -The variants calculated in this script are as follows: "alpha","gamma", "iota", "delta", "b.1.621".
     -If you wish to add more, you would need to add the new variants throughout the script

SUMMARY CALCULATIONS BEING DONE
1. {Variant}_infections = (relative number of variant sequenced) * infections
2. {Variant}_n = (relative number of variant sequenced)*(absolute number of samples sequenced)
3. {Variant}7 = 7 day rolling average
4. {Variant}_low/upci = binomCI(two sided, method jeffreys)
5. {Variant}_Rt = estimate_R + smooth spline

 OUTPUT
 Code will output 3 files:
 1. plot of the Rt values
 2. plot of Rt values with CI
 3. the Rt_estimate.csv
 4. reformatted Rt_estimates (this is specific to the covidtrackerCt.com website)

CITATIONS
Estimate_R
Anne Cori, Neil M. Ferguson, Christophe Fraser, Simon Cauchemez, A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics, American Journal of Epidemiology, Volume 178, Issue 9, 1 November 2013, Pages 1505–1512, https://doi.org/10.1093/aje/kwt133

