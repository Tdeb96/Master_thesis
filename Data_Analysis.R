library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)

#Import datasets for analysis
user_data_prep <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/user_data_prep.RDS")
item_data_prep <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/item_data_prep.RDS")

#Import the uncleaned datasets
Thijs_Observations <- read_delim("~/Google Drive/Data scriptie/Thijs_Observations.csv",";", escape_double = FALSE, trim_ws = TRUE)

#Mail statistics----

#Clickrate of the emails
#Create a mail_offer object combining the mail and the offer id
Thijs_Observations[,"MailOffer"] <- Thijs_Observations %>% unite("MailOffer", c("MAILID" ,"OFFERID"), sep = "_") %>% select("MailOffer")

#Merge these with the userids to get individual mailings
Thijs_Observations[,"UserMailid"] <- Thijs_Observations %>% unite("UserOffer", c("USERID" ,"MAILID"), sep = "_") %>% select("UserOffer")
#Now get statistics for this UserOffer combination
df_raw_email_clicks <- Thijs_Observations %>% group_by(UserMailid) %>% summarize(totalclicks = sum(CLICK), number = n())

#Average amount of offers per email
mean(df_raw_email_clicks$number)
clickrate_emails <- sum(df_raw_email_clicks$totalclicks>0)/nrow(df_raw_email_clicks)

#Charts-----
#Frequency charts of opened emails per user, clicks per user and click rates
df_user_plot <- Thijs_Observations %>% group_by(USERID) %>% summarize(emails = length(unique(MAILID)), clicks = sum(CLICK), click_rate = 100*sum(CLICK)/n())

#Opened emails chart
average_emails <- mean(df_user_plot$emails)
Mails <- ggplot(data = df_user_plot, aes(emails)) + geom_histogram(bins = 30, color = "black", fill = "grey") + 
  geom_vline(aes(xintercept=average_emails), linetype = "dashed") + scale_y_log10(breaks = c(0,10,100,1000,10000,100000)) + xlab("Number of opened emails") + ylab("Frequency")+
  theme_bw()
Mails

#Clicks per user chart
average_clicks <- mean(df_user_plot$clicks)
Clicks <- ggplot(data = df_user_plot, aes(clicks)) + geom_histogram(bins = 30, color = "black", fill = "grey") + 
  geom_vline(aes(xintercept=average_clicks), linetype = "dashed") + scale_y_log10(breaks = c(0,10,100,1000,10000,100000)) + xlab("Number of total clicks") + ylab("Frequency")+
  theme_bw()
Clicks
#Clickrate vs total 
average_click_rate <- mean(df_user_plot$click_rate)
Click_rate <- ggplot(data = df_user_plot, aes(click_rate)) + geom_histogram(bins = 30 , color = "black", fill = "grey") + 
  geom_vline(aes(xintercept=average_click_rate), linetype = "dashed") + scale_y_log10(breaks = c(0,10,100,1000,10000,100000)) + xlab("Clickrates") + ylab("Frequency")+
   scale_x_continuous(labels = function(x) paste0(x, "%"))+
  theme_bw()
Click_rate

#Clicks based on offers 
df_item_plot <- Thijs_Observations %>% group_by(MailOffer) %>% summarize(receivers = length(unique(USERID)), clicks = sum(CLICK), click_rate = 100*sum(CLICK)/n())

#Email recipients 
average_recipients <- mean(df_item_plot$receivers)
Recipients <- ggplot(data = df_item_plot, aes(receivers)) + geom_histogram(bins = 30, color = "black", fill = "grey") + 
  geom_vline(aes(xintercept=average_recipients), linetype = "dashed") + xlab("Number of email recipients") + ylab("Frequency")+
  theme_bw()
Recipients

#Clicks per user chart
average_itemclicks <- mean(df_item_plot$clicks)
ItemClicks <- ggplot(data = df_item_plot, aes(clicks)) + geom_histogram(bins = 30, color = "black", fill = "grey") + 
  geom_vline(aes(xintercept=average_itemclicks), linetype = "dashed")  + xlab("Number of total clicks per offer") + ylab("Frequency")+
  theme_bw()
ItemClicks
#Clickrate
average_click_rate_item <- mean(df_item_plot$click_rate)
Click_rate_item <- ggplot(data = df_item_plot, aes(click_rate)) + geom_histogram(bins = 30 , color = "black", fill = "grey") + 
  geom_vline(aes(xintercept=average_click_rate_item), linetype = "dashed") + xlab("Clickrates") + ylab("Frequency")+
  scale_x_continuous(labels = function(x) paste0(x, "%"))+
  theme_bw()
Click_rate_item

#Scatterplot of items and click rates
ggplot(data = df_item_plot, aes(x = click_rate, y = receivers)) + geom_point( color = "grey")+ xlab("Clickrates") + ylab("Recipients")+
  theme_bw()

#Some item statistics
item_data_prep <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/item_data_prep.RDS")
summary(item_data_prep$STAR_RATING)

#Some user statistics
user_data_prep <- readRDS("~/Dropbox/Uni/Master_Econometrie/Thesis/Code/R/user_data_prep.RDS")

