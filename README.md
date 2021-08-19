# Rt_v2

#purpose: This is the second version of the Rt calculation for covidtrackerct.org
#of variants rewritten to be reproducible by any state. Case numbers pulled from
#from covidestim infection model and frequencies from aggregated sequencing results from state labs

#initial date: 2021/08/19
#author: Tobias Koch 
#email: r.tobiaskoch@gmail.com or toby.koch@yale.edu
#initial author Jessica Rothman

NECESSARY Input FILES
1. estimates.csv (covid infections by state)
You will need to read in cases from downloaded csv from covidestim.org 
keep name estimate.csv and put it into your data_input folder

2. Frequency of variants table
 please use this example as a template for how to format. Frequencies pulled from Grubaugh Lab Nextstrain build.

Nextstrain example:
 https://covidtrackerct.com/variant-surveillance/

Frequency table example:
 [variant_frequencies.csv](https://github.com/rtobiaskoch/Rt_v2/files/7017913/variant_frequencies.csv)
 
 
 OUTPUT
 Code will output 3 files:
 1. a plot of the Rt values
 2. the Rt_estimate.csv
 3. reformatted Rt_estimates (this is specific to the covidtrackerCt.com website)



