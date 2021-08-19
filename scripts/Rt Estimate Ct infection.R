#*******************************************************************************
#####HEADER#####
#*#*******************************************************************************

#purpose: This is the second version of the Rt calculation for covidtrackerct.org
#of variants in the state of connecticut by pulling case numbers 
#from covidestim infection model and sequencing results from the Grubaugh Lab

#initial date: 2021/07/28
#author: Tobias Koch 
#email: r.tobiaskoch#gmail.com and toby.koch@yale.edu
#initial author Jessica Rothman

#PACKAGE INSTALL
packages = c("readr", "tidyverse", "lubridate","EpiEstim", "xlsx","stats","zoo",
             "DescTools", "googledrive","googlesheets4")



## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

rm(packages)

set.seed(1234)

#*******************************************************************************
##### GLAB DATA IMPORT#####
#*#*******************************************************************************

#importing rt data from CT-Yale Variant results googlesheet
#if you don't have access to the sheet or can't get googlesheets to work download
#and import it
var_import <- read_sheet("12xYePgxeF3pi0YiGnDCzmBnPPqZEASuobZ1DXeWZ7QA", sheet = "Rt")

#*******************************************************************************
#####GLAB DATA CLEAN#####
#*#*******************************************************************************

#imported at beginning
var_data = var_import %>%
  rename(alpha_prop = `Freq Alpha`,
         delta_prop = `Freq Delta`,
         gamma_prop = `Freq Gamma`,
         iota_prop = `Freq Iota`,
         b1621_prop = `Freq B.1.621`,
  ) %>% #rename to match variables in rest of code
  filter(!is.na(n)) #removes blank columns

#*******************************************************************************
##### Covidestim.org#####
#*#*******************************************************************************
#will need to read in cases from downloaded csv from covidestim.org 
#that is put it your working directory
#be sure to download the latest version
infect_import = read.csv("estimates.csv")
infect = infect_import %>%
  select(state,
         date,
         infections) %>%
  filter(state == "Connecticut") %>% #filter only CT
  rename(Date = date) %>% #rename for merging with var_data
  mutate(Date = as.Date(Date)) #reformat to match var_data


#*******************************************************************************
#####MERGE DATA#####
#*#*******************************************************************************
var_merge = var_data %>%
  left_join(infect, by = "Date") %>% #merges jhop data with our data and keeps only 
  mutate(alpha_infections        = infections*alpha_prop,
         delta_infections        = infections*delta_prop,
         gamma_infections        = infections*gamma_prop,
         iota_infections    = infections*iota_prop,
         b1621_infections = infections*b1621_prop) %>%
  mutate(alpha_n        = n*alpha_prop,
         delta_n        = n*delta_prop,
         gamma_n        = n*gamma_prop,
         iota_n    = n*iota_prop,
         b1621_n = n*b1621_prop)

#*******************************************************************************
#####ROLLING 7 DAY AVG FOR RT #####
#*******************************************************************************

daily_7<- var_merge %>%
  #7 day rolling avg for samples sequenced
  mutate(alpha_n7 =         rollmean(alpha_n, k = 7, fill = NA),
         gamma_n7 =         rollmean(gamma_n, k = 7, fill = NA),
         delta_n7 =         rollmean(delta_n, k = 7, fill = NA),
         b1621_n7 =  rollmean(b1621_n, k = 7, fill = NA),
         iota_n7 =     rollmean(iota_n, k = 7, fill = NA),
         n_7 =              rollmean(n, k = 7, fill = NA)) %>%
  #7 day rolling avg for frequency(prop aka proportion)
  mutate(alpha_prop7 =         rollmean(alpha_prop, k = 7, fill = NA),
         gamma_prop7 =         rollmean(gamma_prop, k = 7, fill = NA),
         delta_prop7 =         rollmean(delta_prop, k = 7, fill = NA),
         b1621_prop7 =  rollmean(b1621_prop, k = 7, fill = NA),
         iota_prop7 =     rollmean(iota_prop, k = 7, fill = NA))%>%
  #7 day rolling avg calculate by 7 day roll avg frequency pf variant * new infections
  mutate(alpha_infections7 = alpha_prop7 * infections,
         gamma_infections7 = gamma_prop7 * infections,
         delta_infections7 = delta_prop7 * infections,
         b1621_infections7 = b1621_prop7 * infections,
         iota_infections7 = iota_prop7 * infections) %>%
  drop_na #drops future dates and first 3 days because of rollmean 


#*******************************************************************************
#CI FUNCTION
#*#******************************************************************************
ci_fun <- function(v, nn, name, c){
  out = BinomCI(x=v, 
                n=nn, 
                conf.level = 0.95, 
                sides = "two.sided",
                method = "jeffreys")
  
  
  #adds alpha prefix to output names
  cname = paste(name, colnames(out), sep = "_")
  colnames(out) = cname
  
  #binds binom output to daily_7
  out3 = cbind.data.frame(daily_7, out)
  
  #selects columns for only the variant being run for simplicity and Rt function
  out4 = out3 %>% select(Date,
                         infections,
                         contains(name)
  ) %>%
    mutate_at(vars(ends_with("est")), #searches for column with est
              funs(.*infections)
    )%>% #output is <variant_. couldnt figure out how to change
    rename_at(vars(ends_with("est")),
              funs(paste("I"))
              # use this if you need variant prefix funs(paste(name,"I",sep = "_"))
    )
}


alpha_df = ci_fun(daily_7$alpha_n7, daily_7$n_7, "alpha")
gamma_df = ci_fun(daily_7$gamma_n7, daily_7$n_7, "gamma")
delta_df = ci_fun(daily_7$delta_n7, daily_7$n_7, "delta")
b1621_df = ci_fun(daily_7$b1621_n7, daily_7$n_7, "b1621")
iota_df = ci_fun(daily_7$iota_n7, daily_7$n_7, "iota")



#*******************************************************************************
#RT FUNCTION ####
#*#******************************************************************************

#Rt Calculation
#generates Rt calculation, smooths the line then merges with 
#the other variant data in a dataframe

#run for everything -delta
rt_fun= function(df, name){
  
  non0 <- min(which(df$I > 0)) #1st day with infections of variant to start the R estimate otherise R estimate artificially high
  
  df2 = df[non0:nrow(df),] #dataframe filtered where there is the first case of variant to end of dataset
  
  #input of interval for R estimate
  t_start<-seq(2, nrow(df2)-21) 
  t_end<-t_start+21
  config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                             std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                             n1=500,n2=50,t_start=t_start, t_end=t_end)
  )
  
  
  mean_Rt = estimate_R(df2$I, #will search for column named I which was created in the ci_fun but explicitly named here
                       method="uncertain_si",
                       config = config)
  
  #adds back in days that were filtered out to match the days in the main dataframe
  mean_Rt$R$t_start = mean_Rt$R$t_start +non0 
  mean_Rt$R$t_end = mean_Rt$R$t_end + non0
  
  smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
  smooth_spline_lower_ci<-with(mean_Rt$R, smooth.spline(mean_Rt$R$t_end, mean_Rt$R$`Quantile.0.025(R)`,cv=TRUE))
  smooth_spline_upper_ci<-with(mean_Rt$R, smooth.spline(mean_Rt$R$t_end, mean_Rt$R$`Quantile.0.975(R)`,cv=TRUE))
  
  #binds them into a dataframe
  smooth_spline_mean_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y,
                                          smooth_spline_lower_ci$y,smooth_spline_upper_ci$y)
  
  # #renames smooth line Rt so it can merge and is more comprehensible
  smooth_spline_mean_df = rename(smooth_spline_mean_df,
                                 day = `smooth_spline_mean$x`,
                                 Rt = `smooth_spline_mean$y`,
                                 rtlowci = `smooth_spline_lower_ci$y`,
                                 rtupci = `smooth_spline_upper_ci$y`)
  
  #merges the Rt value with the other variant data and renames Rt to have variant suffix
  merge = df %>%
    arrange(Date)%>% #keep in date so that the day variable lines up with the first date
    mutate(day = 1:nrow(df))%>% #used to merge with the estimate_R variable output for the day
    left_join(smooth_spline_mean_df) %>%
    rename_with(.fn = ~paste0(name,"_",.), .cols = c("Rt", "rtlowci", "rtupci")) #renames the smooth_spline output to have variant prefix
}

alpha_rt = rt_fun(alpha_df, "alpha")
gamma_rt = rt_fun(gamma_df, "gamma")
b1621_rt = rt_fun(b1621_df, "b1621")
iota_rt = rt_fun(iota_df, "iota")
delta_rt = rt_fun(delta_df, "delta")



#*******************************************************************************
#RT MERGE FOR EXPORT ####
#*******************************************************************************
#me trying to extract the rt's in the join
rt_list = list(alpha_rt,
               gamma_rt,
               delta_rt,
               b1621_rt,
               iota_rt)

#new merged file that selects only the necessary variables
rt_export <- rt_list %>% 
  reduce(left_join, by = "Date") %>%
  select(Date,
         alpha_Rt, #alpha
         alpha_rtlowci,
         alpha_rtupci,
         gamma_Rt, #gamma
         gamma_rtlowci,
         gamma_rtupci,
         delta_Rt, #delta
         delta_rtlowci,
         delta_rtupci,
         iota_Rt, #iota
         iota_rtlowci,
         iota_rtupci,
         b1621_Rt, #b1621
         b1621_rtlowci,
         b1621_rtupci
  ) %>%
  filter_at(2:ncol(.), #filter all coloumns but the first column which is the date
            any_vars(!is.na(.))) #remove any rows where there is no data

#*******************************************************************************
#Rt Plot####
#*******************************************************************************
#creates a plot as a quick check to see if values are looking correct
rt_plot = rt_export %>%
  select_at(vars(Date, ends_with("_Rt"))) %>% #selects only rt values
  gather("variant", "rt", 2:ncol(.)) #puts into long format for plotting


p <- rt_plot %>%
  ggplot(aes(x = Date, y = rt, color = variant, fill = variant)) +
                        geom_line()
                        
p                       

#*******************************************************************************
#REFORMAT####
#*******************************************************************************
#creating variable for the column rename 
mutate_col = colnames(rt_export[2:length(rt_export)])
#reformat to fit formatting for the website
rt_export2 = rt_export %>%
  mutate_at(all_of(mutate_col),
            funs(round(.,2))   #rounds all variables to 2nd decimal
  ) %>%
  transmute(`First day of week` = Date,
            #ALPHA
            Alpha = alpha_Rt,
            `Alpha-CI` = paste(month(Date),"/",day(Date),
                               ": ", 
                               alpha_Rt,
                               " (",
                               alpha_rtlowci,
                               ", ",
                               alpha_rtupci,
                               ")",sep =""),
            `Alpha-low` = alpha_rtlowci,
            `Alpha-high` = alpha_rtupci,
            #DELTA
            Delta = delta_Rt,
            `Delta-CI` = paste(month(Date),"/",day(Date),
                               ": ", 
                               delta_Rt,
                               " (",
                               delta_rtlowci,
                               ", ",
                               delta_rtupci,
                               ")",sep =""),
            `Delta-low` = delta_rtlowci,
            `Delta-high` = delta_rtupci,
            #GAMMA
            Gamma = gamma_Rt, 
            `Gamma-CI` = paste(month(Date),"/",day(Date),
                               ": ", 
                               gamma_Rt,
                               " (",
                               gamma_rtlowci,
                               ", ",
                               gamma_rtupci,
                               ")",sep =""),
            `Gamma-low` = gamma_rtlowci,
            `Gamma-high` = gamma_rtupci,
            #IOTA
            Iota = iota_Rt, 
            `Iota-CI` = paste(month(Date),"/",day(Date),
                              ": ", 
                              iota_Rt,
                              " (",
                              iota_rtlowci,
                              ", ",
                              iota_rtupci,
                              ")",sep =""),
            `Iota-low` = iota_rtlowci,
            `Iota-high` = iota_rtupci,
            #OTHER
            Other = b1621_Rt, 
            `Other-CI` = paste(month(Date),"/",day(Date),
                               ": ", 
                               b1621_Rt,
                               " (",
                               b1621_rtlowci,
                               ", ",
                               b1621_rtupci,
                               ")",sep =""),
            `Other-low` = b1621_rtlowci,
            `Other-high` = b1621_rtupci,
  )
#*******************************************************************************
#EXPORT####
#*******************************************************************************
#*
#export directly back into the google sheet that we use "CT-Yale Variant results
#creates sheet name with the first and last date in dataframe
sheetname = paste("Rt_out", min(rt_export2$`First day of week`), "to", max(rt_export2$`First day of week`), sep = " ")
#exports to sheet that stores as a history
write_sheet(rt_export2, "1IskdjyQI4N6cCVO2s_hqStltEZ-q091VjAATrovfTjk", sheet = sheetname)

#exports most recent rt data to ct visualizer sheet for others to reference
write_sheet(rt_export2, "12xYePgxeF3pi0YiGnDCzmBnPPqZEASuobZ1DXeWZ7QA", sheet = "rt_R_output")




