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
#####DATA IMPORT#####
#*#*******************************************************************************

#importing rt data from CT-Yale Variant results googlesheet
#if you don't have access to the sheet or can't get googlesheets to work download
#and import it

#GOOGLESHEET DOWNLOAD FROM GRUBAUGH LABS GDRIVE
#note that column names will be slightly different
#var_import <- read_sheet("12xYePgxeF3pi0YiGnDCzmBnPPqZEASuobZ1DXeWZ7QA", sheet = "Rt")

#this is aggregated state variant frequency data
#please see readme for example format
var_import = read.csv("data_input/variant_frequencies.csv")
infect_import = read.csv("data_input/estimates.csv")

#*******************************************************************************
#####GLAB DATA CLEAN#####
#*#*******************************************************************************

#*******************************************************************************
#*******************************************************************************
# #list function 
# 
# var_list = c("B.1.617.2", "AY.3", "AY.4", "AY.25", "AY.26", "AY.33")
# 
# #experimental for loop not working
# for(i in seq_along(var_list)) {
#   x = c(print(i) )
# var_data =  var_import %>%
#     rename_at(vars(ends_with(var_list[i])),
#               funs(paste0(var_list[i],"_prop"))
#     )
# }
# 
# #experimental map function to accomplish below not working
# var_data = var_import %>%
#   nest(variants = -c("Date", "EW", "n"))
# 
# map2(colnames(var_data$variants), var_list, ~(.x %>%
#                                          rename_at(vars(ends_with(.y)),
#                                                    funs(paste0(.y,"_prop"))
#                                                   )
#                                        )
#         )
#*******************************************************************************
#*******************************************************************************


#identifies variant columns and renames them to the format of the rest of previously written code
var_data = var_import %>%
  rename_at(vars(ends_with("B.1.617.2")),
                 funs(paste0("B.1.617.2","_prop"))
            ) %>%
  rename_at(vars(ends_with("AY.3")),
            funs(paste0("AY.3","_prop"))
            ) %>%
  rename_at(vars(ends_with("AY.4")),
            funs(paste0("AY.4","_prop"))
             ) %>%
  rename_at(vars(ends_with("AY.25")),
            funs(paste0("AY.25","_prop"))
            ) %>%
  rename_at(vars(ends_with("AY.26")),
            funs(paste0("AY.26","_prop"))
             ) %>% #rename to match variables in rest of code
  filter(!is.na(n)) %>% #removes blank columns
  mutate(Date = as.Date(Date)) #converts date from string to date for merging

#*******************************************************************************
##### COVID ESTIM DATA CLEAN#####
#*#*******************************************************************************
#will need to read in cases from downloaded csv from covidestim.org 
#that is put it your working directory this is done in line 48 in the data import section
infect = infect_import %>%
  select(state,
         date,
         infections) %>%
  #**************EDIT BY STATE HERE*********
  filter(state == "Connecticut") %>% #filter only CT 
  #*****************************************
  
  rename(Date = date) %>% #rename for merging with var_data
  mutate(Date = as.Date(Date)) #reformat to match var_data


#*******************************************************************************
#####MERGE DATA#####
#*#*******************************************************************************
var_merge = var_data %>%
  left_join(infect, by = "Date") %>% #merges jhop data with our data and keeps only 
  mutate(B.1.617.2_infections        = infections*B.1.617.2_prop,
         AY.3_infections        = infections*AY.3_prop,
         AY.4_infections        = infections*AY.4_prop,
         AY.25_infections    = infections*AY.25_prop,
         AY.26_infections = infections*AY.26_prop) %>%
  mutate(B.1.617.2_n        = n*B.1.617.2_prop,
         AY.3_n        = n*AY.3_prop,
         AY.4_n        = n*AY.4_prop,
         AY.25_n    = n*AY.25_prop,
         AY.26_n = n*AY.26_prop)

#*******************************************************************************
#####ROLLING 7 DAY AVG FOR RT #####
#*******************************************************************************

daily_7<- var_merge %>%
  #7 day rolling avg for samples sequenced
  mutate(B.1.617.2_n7 =         rollmean(B.1.617.2_n, k = 7, fill = NA),
         AY.4_n7 =         rollmean(AY.4_n, k = 7, fill = NA),
         AY.3_n7 =         rollmean(AY.3_n, k = 7, fill = NA),
         AY.26_n7 =  rollmean(AY.26_n, k = 7, fill = NA),
         AY.25_n7 =     rollmean(AY.25_n, k = 7, fill = NA),
         n_7 =              rollmean(n, k = 7, fill = NA)) %>%
  #7 day rolling avg for frequency(prop aka proportion)
  mutate(B.1.617.2_prop7 =         rollmean(B.1.617.2_prop, k = 7, fill = NA),
         AY.4_prop7 =         rollmean(AY.4_prop, k = 7, fill = NA),
         AY.3_prop7 =         rollmean(AY.3_prop, k = 7, fill = NA),
         AY.26_prop7 =  rollmean(AY.26_prop, k = 7, fill = NA),
         AY.25_prop7 =     rollmean(AY.25_prop, k = 7, fill = NA))%>%
  #7 day rolling avg calculate by 7 day roll avg frequency pf variant * new infections
  mutate(B.1.617.2_infections7 = B.1.617.2_prop7 * infections,
         AY.4_infections7 = AY.4_prop7 * infections,
         AY.3_infections7 = AY.3_prop7 * infections,
         AY.26_infections7 = AY.26_prop7 * infections,
         AY.25_infections7 = AY.25_prop7 * infections) %>%
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
  
  
  #adds B.1.617.2 prefix to output names
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


B.1.617.2_df = ci_fun(daily_7$B.1.617.2_n7, daily_7$n_7, "B.1.617.2")
AY.4_df = ci_fun(daily_7$AY.4_n7, daily_7$n_7, "AY.4")
AY.3_df = ci_fun(daily_7$AY.3_n7, daily_7$n_7, "AY.3")
AY.26_df = ci_fun(daily_7$AY.26_n7, daily_7$n_7, "AY.26")
AY.25_df = ci_fun(daily_7$AY.25_n7, daily_7$n_7, "AY.25")

#*******************************************************************************
#RT FUNCTION ####
#*#******************************************************************************

#Rt Calculation
#generates Rt calculation, smooths the line then merges with 
#the other variant data in a dataframe

#run for everything -AY.3
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

B.1.617.2_rt = rt_fun(B.1.617.2_df, "B.1.617.2")
AY.4_rt = rt_fun(AY.4_df, "AY.4")
AY.26_rt = rt_fun(AY.26_df, "AY.26")
AY.25_rt = rt_fun(AY.25_df, "AY.25")
AY.3_rt = rt_fun(AY.3_df, "AY.3")



#*******************************************************************************
#RT MERGE FOR EXPORT ####
#*******************************************************************************
#me trying to extract the rt's in the join
rt_list = list(B.1.617.2_rt,
               AY.4_rt,
               AY.3_rt,
               AY.26_rt,
               AY.25_rt)

#new merged file that selects only the necessary variables
rt_export <- rt_list %>% 
  reduce(left_join, by = "Date") %>%
  select(Date,
         B.1.617.2_Rt, #B.1.617.2
         B.1.617.2_rtlowci,
         B.1.617.2_rtupci,
         AY.4_Rt, #AY.4
         AY.4_rtlowci,
         AY.4_rtupci,
         AY.3_Rt, #AY.3
         AY.3_rtlowci,
         AY.3_rtupci,
         AY.25_Rt, #AY.25
         AY.25_rtlowci,
         AY.25_rtupci,
         AY.26_Rt, #AY.26
         AY.26_rtlowci,
         AY.26_rtupci
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

rt_plot = rt_export %>%
  gather("variant", "value", -Date) %>% #puts into long format for plotting
  separate(variant, c("variant", "stat"), sep = "_") %>% #separates the variant name by the stat to spread
  spread(stat, value)


p <- rt_plot %>%
  ggplot(aes(x = Date, y = Rt, color = variant, fill = variant)) +
        geom_line()

p_ci <- p +
        geom_ribbon(aes(ymin = rtlowci, ymax = rtupci, fill = variant), 
                                               alpha=0.1, 
                                              linetype="dashed",
                                               color="grey")
                        

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
            #B.1.617.2
            B.1.617.2 = B.1.617.2_Rt,
            `B.1.617.2-CI` = paste(month(Date),"/",day(Date),
                               ": ", 
                               B.1.617.2_Rt,
                               " (",
                               B.1.617.2_rtlowci,
                               ", ",
                               B.1.617.2_rtupci,
                               ")",sep =""),
            `B.1.617.2-low` = B.1.617.2_rtlowci,
            `B.1.617.2-high` = B.1.617.2_rtupci,
            #AY.3
            AY.3 = AY.3_Rt,
            `AY.3-CI` = paste(month(Date),"/",day(Date),
                               ": ", 
                               AY.3_Rt,
                               " (",
                               AY.3_rtlowci,
                               ", ",
                               AY.3_rtupci,
                               ")",sep =""),
            `AY.3-low` = AY.3_rtlowci,
            `AY.3-high` = AY.3_rtupci,
            #AY.4
            AY.4 = AY.4_Rt, 
            `AY.4-CI` = paste(month(Date),"/",day(Date),
                               ": ", 
                               AY.4_Rt,
                               " (",
                               AY.4_rtlowci,
                               ", ",
                               AY.4_rtupci,
                               ")",sep =""),
            `AY.4-low` = AY.4_rtlowci,
            `AY.4-high` = AY.4_rtupci,
            #AY.25
            AY.25 = AY.25_Rt, 
            `AY.25-CI` = paste(month(Date),"/",day(Date),
                              ": ", 
                              AY.25_Rt,
                              " (",
                              AY.25_rtlowci,
                              ", ",
                              AY.25_rtupci,
                              ")",sep =""),
            `AY.25-low` = AY.25_rtlowci,
            `AY.25-high` = AY.25_rtupci,
            #OTHER
            Other = AY.26_Rt, 
            `Other-CI` = paste(month(Date),"/",day(Date),
                               ": ", 
                               AY.26_Rt,
                               " (",
                               AY.26_rtlowci,
                               ", ",
                               AY.26_rtupci,
                               ")",sep =""),
            `Other-low` = AY.26_rtlowci,
            `Other-high` = AY.26_rtupci,
  )
#*******************************************************************************
#EXPORT####
#*******************************************************************************
#export directly back into the google sheet that we use "CT-Yale Variant results
#creates sheet name with the first and last date in dataframe

rt_name = paste("data_output/Rt_estimates ", 
                min(rt_export2$`First day of week`), 
                " to ", max(rt_export2$`First day of week`),
                ".csv", 
                sep = "")

rt_reformat_name = paste("data_output/Rt_estimate_reformat ", 
                         min(rt_export2$`First day of week`), 
                         " to ", 
                         max(rt_export2$`First day of week`),
                         ".csv", 
                         sep = "")

#writes csv for rt estimates
write.csv(rt_export, rt_name)

#writes csv for reformatted
write.csv(rt_export2, rt_reformat_name)

#print the plot for double check
png("data_output/rt_plot.png", width = 800, height = 400)
p
dev.off()







