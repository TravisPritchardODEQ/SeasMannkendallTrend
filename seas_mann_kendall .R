library(tidyverse)
library(lubridate)
library(wql)
library(zoo)
library(readxl)
library(AWQMSdata)


# 2/1/2019
#Modified to use data pulled directly from AWQMS


# 4/5/2018

# This script takes 7DMADMAX data exported from volWQdb and 
# performs a seasonal mann kendall trend test.
#  Script generates a csv table with results of the seasonal mann kendall test
#  as well as graphs.
#  the slope on the graphs is generated using the reported sen's slope

# Trends are analyzed using WQL library. This analysis disregards NA values.
#  


# from https://www.epa.gov/sites/production/files/2016-05/documents/tech_notes_6_dec2013_trend.pdf :
#  The seasonal Kendall test statistic is computed by performing a Mann-Kendall calculation
#  for each season, then combining the results for each season. For monthly seasons, January
#  observations are compared only to other January observations, etc. No comparisons are
#  made across seasonal boundaries. The Seasonal Kendall test is highly robust and relatively
#  powerful, and is often the recommended method for most water quality trend monitoring.
#  The Sk statistic is computed as the sum of the S from each season

# In this script, seasons are set as July and August, since those months have the most complete
# dataset. Trends are only run on months with <1 missing day

# Trend is considered signifigant when p-value <0.10 



# Load data ---------------------------------------------------------------


#Import t_results export from database
t_results <- AWQMS_Data(startdate = '1949-09-15', 
                        enddate = '2019-02-01', 
                        station = c('30143-ORDEQ', '30147-ORDEQ', 
                                    '30154-ORDEQ', '30161-ORDEQ', 
                                    '37477-ORDEQ'), 
                        char = c('Temperature, water') , 
                        stat_base = c('7DADM') , 
                        crit_codes = TRUE)


#tmp_strd <- read.csv("Ump_temp_strd.csv")

# Load critical tau values to assess signifigance
load("tau_crit_values.Rdata")


#create empty list to acept test results
kendall_list <- list()




# Calculate monthly average 7 day average daily maximum temperature -------


#crate table with average SDADM per month
sdadm <- t_results %>%
  filter(QualifierAbbr %in% c("DQL=A", "DQL=B") | 
           is.na(QualifierAbbr)) %>%
  mutate(date = ymd(SampleStartDate),
         month = month(date),  
         yrmon = as.yearmon(date), 
         year = year(date)) %>%
  group_by(MLocID, year, month) %>%
  summarise(sdadm = mean(Result_Numeric),
            count = n()) %>%
  ungroup() %>%
  complete(MLocID, year, month) %>%
  filter(month == 7 | month == 8)

# sdadm_wide <- sdadm %>%
#   select(MLocID, year, month, sdadm) %>%
#   spread(month, sdadm)

sdadm_count <- sdadm %>%
  mutate(mo_days = days_in_month(month),
         diff = mo_days - count,
         exclude = ifelse(diff > 1 | is.na(sdadm), "yes", "no")) %>%
  filter(month != 10) %>%
  select(MLocID, year, month, exclude)

#add exclusion determination and filter out months to exclude
#format toable to fit WQdata input format
sdadm_assess <- sdadm %>%
  left_join(sdadm_count, by = c("MLocID", "year", "month"))%>%
  filter(exclude == "no") %>%
  select(MLocID, year, month, sdadm) %>%
  mutate(date = paste0(year, "-", month, "-01")) %>%
  mutate(date = ymd(date)) %>%
  mutate(depth = 1,
         variable = as.factor("TEMP")) %>%
  select(date, MLocID, depth, variable, sdadm) %>%
  rename(value = sdadm,
         site = MLocID)

# Run Seasonal Mann-Kendall Trend analysis --------------------------------


# for loop to process each MLocID
for(i in 1:length(unique(sdadm_assess$site))){
  
  #select unique site ID and arrange in chonological order
  #fill in NAs for missing months
  sdadm_seamannkenn <- sdadm_assess %>%
    filter(site == unique(sdadm_assess$site)[i]) %>%
    filter(month(date) == 7 | month(date) == 8) 
  #create wqdata format table
  wqdataframe <- wqData(sdadm_seamannkenn, 1:3, 4:5, site.order = TRUE, type = "long",
                        time.format = "%y-%m-%d")
 
  #create time series
  timeseries <- tsMake(wqdataframe, focus= gsub("-", ".", paste0("s",unique(sdadm_assess$site)[i])))
  
  res <- seaKen(timeseries)
  
  #save results as a dataframe and add MLocID and n
  df <- data.frame(as.list(unclass(res))) %>%
    mutate(MLocID = unique(sdadm_assess$site)[i],
           n = length(timeseries))
  
  
 #bind results list
   kendall_list[[i]] <- df
   
   
   

  
} #end of for loop

#bind results of seasonal mann kendall test to single dataframe
kendall_results <- bind_rows(kendall_list) 


#trend is detected when the tau > critical value and p. value (sl) is <0.10
#reorder and rename rows for better exporting
kendall_results <- kendall_results %>%
  mutate(significance = ifelse(p.value < 0.10 & sen.slope < 0, "Signifigant (-)", 
                               ifelse(p.value < 0.10 & sen.slope > 0, "Signifigant (+)",
                                      "No Trend"))) %>%
  select(MLocID, significance, p.value, sen.slope)


write.csv(kendall_results, "Umpqua ref seasonal mann kendall results.csv", row.names = FALSE)



# Set up graphing data ----------------------------------------------------


#join trend analysis to data
sdadm_trend <- sdadm %>%
  left_join(kendall_results, by = "MLocID") %>%
  left_join(sdadm_count, by = c("MLocID", "year", "month")) %>%
  filter(exclude == "no") %>%
  select(-exclude)

#set up some data for graphing
sdadm_raw_trend <-t_results %>%
  filter(QualifierAbbr %in% c("DQL=A", "DQL=B") | 
     is.na(QualifierAbbr)) %>%
  mutate(date = ymd(SampleStartDate),
         month = month(date),  
         yrmon = as.yearmon(date), 
         year = year(date),
         moname = month.name[month]) %>%
  left_join(kendall_results, by = "MLocID") %>%
  left_join(Temp_crit, by = c("FishCode" = "FishUse_code")) %>%
  filter(month(date) == 7 | month(date) == 8) %>%
  left_join(sdadm_count, by = c("MLocID", "year", "month")) %>%
  filter(exclude == "no") %>%
  select(-exclude)
  




# Graphing ----------------------------------------------------------------



for(j in 1:length(unique(sdadm_raw_trend$MLocID))){
  
#box plots
  sdadm_box_plot <- sdadm_raw_trend %>%
    filter(MLocID == unique(sdadm_raw_trend$MLocID)[j]) %>%
    filter(month == 7 | month(date) == 8)
  
  strd <- sdadm_box_plot$Temp_Criteria[[1]]
  
  
  box <- ggplot(sdadm_box_plot, aes(x = year, y = Result_Numeric, group = factor(yrmon))) +
    geom_boxplot(aes(fill = factor(moname))) +
    scale_fill_brewer(palette="Set1") +
    scale_x_continuous(breaks = seq(min(sdadm_box_plot$year), max(sdadm_box_plot$year, by = 3)), 
                       minor_breaks = 1) +
    coord_cartesian(ylim = c(12,30)) +
    geom_hline(aes(yintercept = strd), color = "red", linetype = "dashed") +
    annotate("text", label = "Water Quality Standard", 
             x = 2008,
             y = strd -  0.5,
             colour = "red", size = 3.5)+
    labs(title = "7 Day Average Daily Maximum Temperature",
         subtitle =  unique(sdadm_raw_trend$SiteDescription)[j],
         x = "Year",
         y = "Temperature (degrees C)")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle =  element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 50, hjust = 1)) +
    guides(fill=guide_legend(title = "Month"))
  
  
  
  
    
    ggsave(box, file=paste("Graphs/Box/",unique(sdadm_raw_trend$MLocID)[j], "- Box.png"), 
           width = 8, height = 5, units = c("in"))
    
   
   
    
     #point graph
    sdadm_month_average <- sdadm_trend %>%
      filter(MLocID == unique(sdadm_raw_trend$MLocID)[j]) %>%
      mutate(yearmon = as.Date(paste0(year,"-",month,"-1"))) %>%
      mutate(moname = month.name[month]) %>%
      filter(!is.na(sdadm))
      
    
    p <- ggplot(data = sdadm_month_average)+
      geom_point(aes(x = year, y = sdadm, color = factor(moname)), size = 2, position = position_dodge(0.15)) +
      scale_x_continuous(breaks = seq(min(sdadm_month_average$year), max(sdadm_month_average$year),1),
                         minor_breaks = 1) +
      scale_color_brewer(palette="Set1") +
      coord_cartesian(ylim = c(12,30)) +
      labs(title = "Average 7 Day Average Daily Maximum Temperature",
           subtitle = paste0(unique(sdadm_raw_trend$SiteDescription)[j]),
           x = element_blank(),
           y = "Temperature (degrees C)")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5, size = 12),
            plot.subtitle =  element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 50, hjust = 1)) +
      guides(color = guide_legend(title = "Month"), order = 1)
   
    #graph no trend values
    if(sdadm_month_average$significance[1] == "No Trend"){
      p = p + annotate("text", label = "No Trend", 
                           x = 2009,
                           y = 30,
                           colour = "black", size = 3.5)
      
    } #end of no trend if statement
    
    if(sdadm_month_average$significance[1] != "No Trend"){
      p = p +  annotate("text", label = "Significant Trend (p-value < 0.10)", 
                            x = 2009,
                            y = 30,
                            colour = "black", size = 3.5) 
      
      
      
      
      #plot trend line. This is taken from the coho trends process
      slope <- kendall_results[kendall_results$MLocID == unique(sdadm_raw_trend$MLocID)[j], "sen.slope"]
      x.delta <- as.numeric((max(sdadm_month_average$year) - min(sdadm_month_average$year)))/2
      SK.min <- mean(sdadm_month_average$sdadm, na.rm = TRUE) - x.delta*slope
      SK.max <- mean(sdadm_month_average$sdadm, na.rm = TRUE) + x.delta*slope
      
      p <- p + geom_segment(aes(x = min(sdadm_month_average$year), y = SK.min,
                                xend = max(sdadm_month_average$year+0.5), yend = SK.max, linetype = "Trend"),
                                size = 1.05, color = "gray49") +
        scale_linetype_manual(values=c("dashed")) +
        guides(linetype=guide_legend(title = element_blank(), order = 2))
      
      
    }
    
    
    
   

    
    ggsave(p, file=paste("Graphs/Average/",unique(sdadm_raw_trend$MLocID)[j], "- average.png"), 
           width = 8, height = 5, units = c("in"))
    

  } #end of MLocID loop
  
  
  
  




