#current testing code####
#*************************************
variant_list = c("alpha", "delta", "gamma","iota", "b.1.621")

var_rename_fun = function(df, variant){
  df2 = df %>%
    select_at(vars(contains(variant))
    ) %>%
    rename_at(vars(contains(variant)),
              funs(paste(variant,"_prop",sep = ""))
    ) %>%
    cbind.data.frame(df, df2)
}

var_data = for(i in seq_along(variant_list)) {
  var_rename_fun(var_import, variant_list[i])
}

df %>%
  filter(!is.na(n)) %>% 
  mutate(Date = as.Date(Date)) %>%



#*************************************
#*
#*
#delta
non0 <- min(which(delta_df$delta_est > 0))

t_beg = if(non0>2){
  non0-2
}else{2}

t_start<-seq(t_beg,length(delta_df$delta_est)-18)
t_end<- 

case7_R = delta_df$delta_est[t_beg:length(delta_df$delta_est)] 

#alpha
non02 <- min(which(alpha_df$alpha_est > 0))

t_beg2 = if(non0>2){
  non02-2
}else{2}

t_start2<-seq(t_beg2,length(alpha_df$alpha_est)-18)
t_end2<-t_start2+18




config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end)
)

case7_R2 = alpha_df$alpha_est[t_beg:length(alpha_df$alpha_est)]                  

mean_Rt = estimate_R(case7_R,
                     method="uncertain_si",
                     config = config)

smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
smooth_spline_mean_alpha_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y)

#renames smooth line Rt so it can merge and is more comprehensible
smooth_spline_mean_df = rename(smooth_spline_mean_df, 
                               day =`smooth_spline_mean$x`,
                               Rt = `smooth_spline_mean$y`)

#merges the Rt value with the other variant data and renames Rt to have variant suffix
merge = delta_df %>%
  arrange(Date)%>%
  mutate(day = 1:nrow(delta_df))%>%
  left_join(smooth_spline_mean_df) %>%
  rename_with(.fn = ~paste0("alpha","_",.), .cols = Rt )


#*************************************************************************************************************************


  
  non0 <- min(which(delta_df$delta_cases7 > 0))
  
  t_beg = if(non0>2){
    non0
  }else{2}
  
  t_start<-seq(t_beg,length(delta_df$delta_cases7)-19)
  t_end<-t_start+21
  config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                             std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                             n1=500,n2=50,t_start=t_start, t_end=t_end)
                        )
                        
              case7_R = delta_df$delta_cases7[t_beg:length(delta_df$delta_cases7)]                  
                        
                        mean_Rt = estimate_R(case7_R,
                                             method="uncertain_si",
                                             config = config)
                        
                        smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
                        smooth_spline_mean_alpha_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y)
                        
                        #renames smooth line Rt so it can merge and is more comprehensible
                        smooth_spline_mean_df = rename(smooth_spline_mean_df, 
                                                       day =`smooth_spline_mean$x`,
                                                       Rt = `smooth_spline_mean$y`)
                        
                        #merges the Rt value with the other variant data and renames Rt to have variant suffix
                        merge = delta_df %>%
                          arrange(Date)%>%
                          mutate(day = 1:nrow(delta_df))%>%
                          left_join(smooth_spline_mean_df) %>%
                          rename_with(.fn = ~paste0("alpha","_",.), .cols = Rt )


#*************************************************************************************************************************
                        #works but day is off for merge
                        c7 = delta_df$delta_est
                        
                        non0 <- min(which(c7 > 0))
                        
                        c7 = c7[non0:length(c7)]
                
                      
                        t_start<-seq(2,length(c7)-21)
                        t_end<-t_start+21
                        config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                                                   std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                                                   n1=500,n2=50,t_start=t_start, t_end=t_end)
                        )
            
                        
                        mean_Rt = estimate_R(c7,
                                             method="uncertain_si",
                                             config = config)
                        
                        #trying to fix day merge
                        
                        c7 = delta_df %>% 
                          select(delta_est, Date) %>%
                          rename (I = delta_est)
                        
                        non0 <- min(which(c7[,1] > 0))
                        
                        # c7 = c7[non0:nrow(c7)]
                        
                        
                        t_start<-seq(non0, nrow(c7)-21)
                        t_end<-t_start+21
                        config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                                                   std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                                                   n1=500,n2=50,t_start=t_start, t_end=t_end)
                        )
                        
                        
                        mean_Rt = estimate_R(c7,
                                             method="uncertain_si",
                                             config = config)
                        
                        
                        #3
                        
                        c7 = delta_df %>%
                          rename_at(vars(ends_with("est")),
                                    funs(paste("I"))
                                    )
                        
                        non0 <- min(which(c7[,1] > 0))
                        
                        # c7 = c7[non0:nrow(c7)]
                        
                        
                        t_start<-seq(non0, nrow(c7)-21)
                        t_end<-t_start+21
                        config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                                                   std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                                                   n1=500,n2=50,t_start=t_start, t_end=t_end)
                        )
                        
                        
                        mean_Rt = estimate_R(c7,
                                             method="uncertain_si",
                                             config = config)