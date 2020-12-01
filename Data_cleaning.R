library(tidyverse)
library(lubridate)


#Cleaning of the user data------
#One of the goals is to obtain one row per user
user_data <- read_delim("~/Google Drive/Data scriptie/Thijs_WebsiteClicks.csv", ";", escape_double = FALSE, trim_ws = TRUE)
Observations <- read_delim("~/Google Drive/Data scriptie/Thijs_Observations.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)
#Remove users not in the dataset, Observations is the dataset with the click data
users <- unique(Observations$USERID)
user_data <- user_data[which(user_data$USERID %in% users),]

#Sort on user_ID
user_data<- user_data[order(user_data$USERID),]

#Webpage destination - Assumption: Country is the most important
#extract the country
user_data$COUNTRY <- sapply(user_data$WEBPAGE_DESTINATION,function(x){strsplit(x, "[_]")[[1]][1]})
user_data$COUNTRY2 <- sapply(user_data$WEBPAGE_DESTINATION,function(x){strsplit(x, "[_]")[[1]][2]})
user_data$COUNTRY3 <- sapply(user_data$WEBPAGE_DESTINATION,function(x){strsplit(x, "[_]")[[1]][3]})

#Drop missing observations
user_data <- drop_na(user_data)

#Remove all the rows corresponding to columns showing null or undefined, as these correspond to errorous datapoints
user_data<- user_data[user_data$COUNTRY3!="null",]
user_data<- user_data[user_data$COUNTRY3!="undefined",]

#Manually drop remaining data errors
user_data<- user_data[user_data$COUNTRY!="41F43E44044244343343043B43844F",]
user_data<- user_data[user_data$COUNTRY!="41344043544643844F",]

original <- unique(user_data$COUNTRY)
#make a hardcode key, done with help of dput
new <- c("Griekenland", "Turkije", "Egypte", "Spanje", "Portugal", "Cyprus", 
         "Italië", "Verenigde Arabische Emiraten", "Griekenland", "Spanje", 
         "Turkije", "Egypte", "Tunesië", "Tunesië", "Bulgarije", 
         "Montenegro", "Portugal", "Marokko", "Malta", "Kaapverdië", 
         "Kroatië", "Israël", "Bulgarije", "Kroatië", "Marokko", 
         NA, "Israël", "Kaapverdië", "Cyprus", "Griekenland", "Italië", 
         "Malta", "Montenegro", "Verenigde Arabische Emiraten", NA, 
         NA, "Spanje", "Griekenland", NA, "Gambia", "Marokko", 
         NA, NA, "Griekenland", "Turkije", "Griekenland", "Egypte", 
         "Spanje", "Portugal", "Spanje", NA, "Griekenland", 
         "Spanje", "Marokko", "Portugal", "Egypte", "Griekenland", "Turkije", 
         "Egypte", "Spanje", "Griekenland", "Tunesië", "Marokko", "Cyprus", 
         "Griekenland", "Turkije", "Tunesië", "Spanje", "Italië", 
         "Bulgarije", "Griekenland", "Italië", "Kroatië")
key <- data.frame("original" = unique(user_data$COUNTRY), "new" = new)

user_data <- left_join(user_data, key, by =c("COUNTRY"="original"))
user_data$COUNTRY <- user_data$new

#Drop the few unnecessary columns
user_data <- user_data %>% select(-c("COUNTRY2", "COUNTRY3", "new"))
user_data <- drop_na(user_data)

#Recode the Webpage to a indicator of seriousness for the offer
original <- unique(user_data$WEBPAGE)
WEBPAGE_order <- c(2, 4, 3, 
         1, 5, 7, 6, 
         8, 9, 10, 
         2)
key <- data.frame("original" = original, "WEBPAGE_order" = WEBPAGE_order)
user_data <- left_join(user_data, key, by =c("WEBPAGE"="original"))

#Find the most_recent session for each user and furthest shopping for each user
max_time <- user_data %>% group_by(USERID) %>% summarize("most_recent_session" = max(SESSION_DATETIME), "Furthest_shopping" = max(WEBPAGE_order))
  #Merge with the original matrix
  user_data <- left_join(user_data, max_time)

#Now we move to the final format, in which we for each user extract the destination from the most recent session and the destination from the furthest shopping
#If there is no furthest shopping we simply extract the most_recent, furthest shopping.
  
  #First we format a final database
  user_final <- user_data %>% group_by(USERID) %>% summarize(n(), Max_session = max(SESSION_INDEX)+1) %>%select(-("n()"))

  #Final destination of the most recent session
  
  #First order the matrix first on user, then by session time
  user_data <- user_data[order(user_data$USERID, user_data$SESSION_DATETIME),]
  #extract matrix with most_recent sessions
  most_recent <- user_data[user_data$SESSION_DATETIME==user_data$most_recent_session,]
  max_index <- most_recent %>% group_by(USERID) %>% summarize("most_recent_sessionindex" = max(SESSION_INDEX))
  most_recent <- left_join(most_recent, max_index)  
  
  #select only the last action
  most_recent <- most_recent[most_recent$most_recent_sessionindex==most_recent$SESSION_INDEX,]
  
  #Some users have errorous data, multiple max clicks, so we select only the most_recent website clicktime
  max_index <- most_recent %>% group_by(USERID) %>% summarize("most_recent_clickid" = max(WEBSITECLICKID))
  most_recent <- left_join(most_recent, max_index)  
  #select only the last clickid
  most_recent <- most_recent[most_recent$most_recent_clickid==most_recent$WEBSITECLICKID,]
  
  #select relevant variables
  most_recent <- most_recent %>% select(c("USERID","SESSION_DATETIME", "SESSION_INDEX", "WEBPAGE","COUNTRY"))
  most_recent <- most_recent %>% rename( "MOST_RECENT_SESSIONDATETIME"= "SESSION_DATETIME", "MOST_RECENT_SESSION_INDEX" = "SESSION_INDEX",
                                         "MOST_RECENT_WEBPAGE"="WEBPAGE", "MOST_RECENT_COUNTRY"="COUNTRY")
  #Paste variables of interest in the final database
  user_final <- left_join(user_final, most_recent)
  
  #Furthest shopping.
  
  #select furthest clicks
  furthest <- user_data[user_data$WEBPAGE_order==user_data$Furthest_shopping,]
  #find the latest of these clicks
  max_index <- furthest %>% group_by(USERID) %>% summarize("most_recent_clickid" = max(WEBSITECLICKID))
  furthest <- left_join(furthest, max_index)  
  #select only the last clickid
  furthest <- furthest[furthest$most_recent_clickid==furthest$WEBSITECLICKID,]
  
  #select relevant variables
  furthest <- furthest %>% select(c("USERID","WEBPAGE","COUNTRY"))
  furthest <- furthest %>% rename( "FURTHEST_WEBPAGE"="WEBPAGE", "FURTHEST_COUNTRY"="COUNTRY")
  #Paste variables of interest in the final database
  user_final <- left_join(user_final, furthest)
  
#Extract variables from these items
  #Quarter of latest session
  user_final$MOST_RECENT_QUARTER <- quarter(user_final$MOST_RECENT_SESSIONDATETIME)
  user_final$MOST_RECENT_YEAR <- year(user_final$MOST_RECENT_SESSIONDATETIME)
  
  #recode session to start at 1 (for effects calculation (0*something))
  user_final$MOST_RECENT_SESSION_INDEX =user_final$MOST_RECENT_SESSION_INDEX+ 1
  
  #Recode year and quarter into factos
  user_final$MOST_RECENT_YEAR <- as.factor(user_final$MOST_RECENT_YEAR)
  user_final$MOST_RECENT_QUARTER <- as.factor(user_final$MOST_RECENT_QUARTER)
  
  #turn the relevant variables into dummies
  user_final <-  fastDummies::dummy_cols(user_final, remove_selected_columns = T)
  
  #Drop the non-relevant columns
  user_final <- user_final %>% select(-c("MOST_RECENT_SESSIONDATETIME"))

#save dataset
saveRDS(user_final, "user_final.RDS")
#Clear environment
rm(list = ls())

#cleaning of item data------
item_data <- read_delim("~/Google Drive/Data scriptie/Thijs_OfferDetails.csv", ";", escape_double = FALSE, trim_ws = TRUE)
Observations <- read_delim("~/Google Drive/Data scriptie/Thijs_Observations.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)
#Merge the mail_offer id's
item_data[,"MailOffer"] <- item_data %>% unite("MailOffer", c("MAILID" ,"OFFERID"), sep = "_") %>%
  select("MailOffer")
#reorder columns
item_data <- item_data %>% select(-c("MAILID", "OFFERID"))%>% select("MailOffer", everything())

#remove all mail_offer combinations not in observations
items <- Observations %>% unite("MailOffer", c("MAILID" ,"OFFERID"), sep = "_") %>%
  select("MailOffer") %>% unique()
items <- unlist(items)
item_data <- item_data[which(item_data$MailOffer %in% items),]
rm(items)
rm(Observations)

#Split the date of the mail into a year and a quarter
item_data$Year <- year(as.POSIXct(item_data$MAIL_DATETIME, format = "%d-%m-%Y"))
item_data$Quarter <- quarter(as.POSIXct(item_data$MAIL_DATETIME, format = "%d-%m-%Y"))
item_data <- item_data %>% select(-c("MAIL_DATETIME"))

#Replace the values of "Costa del sol" into spain and "Roda" into greece
item_data[item_data$COUNTRY_NAME=="Costa del Sol", "COUNTRY_NAME"] <- "Spanje"
item_data[item_data$COUNTRY_NAME=="Roda", "COUNTRY_NAME"] <- "Griekenland"

#Remove region, city and acconomdation name
item_data <- item_data %>% select(-c("REGION_NAME", "CITY_NAME", "ACCOMMODATION_NAME"))
item_data <- item_data %>% select(-c("USP1", "USP2", "USP3"))

#Cleaning room occupancy
ROOM_OCCUPANCY <- unique(item_data$ROOM_OCCUPANCY)
new <- c("2", "3", "3_k", 
         "2", "2", "4", "2_k", 
         "5", "2", "3_k", 
         "1", "2_k", "2_k", 
         "2_k", "3", 
         "2_k", "2_k", 
         "2_k", "2_k", "2_k", 
         "2", "3_k", "2_k", 
         "2_k", "2_k", "2_k", 
         "3_k", "3_k", "2_k", 
         "2_k", "3_k", 
         "3_k", "3_k", 
         "2_k", "2", "2_k", 
         "2_k", "2_k", 
         "2_k", "2_k", 
         "2_k", "2_k", 
         "2_k", "2_k", 
         "2_k", "2_k", "2_k", "2_k", 
         "2_k", "2_k", "2_k", 
         "2_k", "2_k", "2_k", 
         "2_k", "2_k"
)
key <- data.frame("ROOM_OCCUPANCY" = ROOM_OCCUPANCY, "maxpers" = new )
item_data <- left_join(item_data, key)

#create a dummy for children
item_data$Childdummy <- NA
item_data[,(ncol(item_data)-1):ncol(item_data)]<- str_split_fixed(item_data$maxpers, "_", 2)
item_data$maxpers <- as.numeric(item_data$maxpers)
item_data[item_data$Childdummy=="","Childdummy"] <-0
item_data[item_data$Childdummy=="k","Childdummy"] <-1
item_data$Childdummy <- as.numeric(item_data$Childdummy)
item_data <- item_data %>% select(-c("ROOM_OCCUPANCY"))

#Define extraction pattern for departure month
pattern <- "(de)|(ap)|(ja)|(no)|(Mar)|(me)|(Mar)|(Oct)|(Feb)|(Jan)|(Apr)|(Dec)|(Nov)|(Sep)|(Aug)|(May)|(Jun)|(Jul)|(januari)|(january)|(jan)|(februari)|(february)|(feb)|(maart)|(march)|(mar)|(april)|(apr)|(mei)|(may)|(june)|(juni)|(jun)|(jul)|(juli)|(july)|(aug)|(augustus)|(august)|(okt)|(oct)|(october)|(oktober)|(sep)|(september)|(november)|(nov)|(dec)|(december)" 
item_data$dep_month <- sapply(item_data$DEPARTURE_DATE, function(x){str_extract(x,pattern)})
item_data[is.na(item_data$dep_month),"DEPARTURE_DATE"]

#uniform all the various months into one single format
original <- unique(item_data$dep_month)
Dep_month <- c("mei", "jun", "jul", "jan", "apr", "feb", "jun", "mei", 
               "jul", "aug", "aug", "maa", "apr", "sep", "okt", "sep", "dec", 
               "jan", "nov", "okt", "nov", "dec", "feb", "maa", NA)
key <- data.frame("dep_month" = original, "Dep_month" = Dep_month)
item_data <- left_join(item_data, key) %>% select(-c("dep_month", "DEPARTURE_DATE"))

#clean the meal_plan column
original <- unique(item_data$MEAL_PLAN)
Meal_plan <- c("All inclusive", "All inclusive", "Logies", "Logies en ontbijt", 
               "Logies en ontbijt", "Ultra all inclusive", "Halfpension", "Ultra all inclusive", 
               "Logies", "Halfpension", "Logies en ontbijt", "Ultra all inclusive", 
               "Ultra all inclusive", "All inclusive", "Logies en ontbijt", "Ultra all inclusive", 
               "All inclusive", "Volpension", "Logies en ontbijt"
)
key <- data.frame("MEAL_PLAN" = original, "Meal_plan" = Meal_plan)
item_data <- left_join(item_data, key) %>% select(-c("MEAL_PLAN"))

#drop roomtype
item_data <- item_data %>% select(-c("ROOMTYPE"))

#normalize star rating
item_data$STAR_RATING[item_data$STAR_RATING>5] <- item_data$STAR_RATING[item_data$STAR_RATING>5]/10

#replace zero star rating and review rating by na
item_data$STAR_RATING[item_data$STAR_RATING==0] <- NA
item_data$REVIEW_RATING[item_data$REVIEW_RATING==0] <- NA

#Recode commas to points
item_data$REVIEW_RATING <- as.numeric(gsub(",", ".", gsub("\\.", "", item_data$REVIEW_RATING)))
#Drop the NA's ASSUMPTION, MOET NOG NAAR GEKEKEN WORDEN
item_data <- drop_na(item_data)

item_data$REVIEW_RATING[item_data$REVIEW_RATING>10] <- item_data$REVIEW_RATING[item_data$REVIEW_RATING>10]/10

#Drop the NA's ASSUMPTION, MOET NOG NAAR GEKEKEN WORDEN
item_data <- drop_na(item_data)

#Change variables to factors for dummy creation
item_data$PRICE <- as.numeric(item_data$PRICE)
item_data$PRICE_ORIGINAL <- as.numeric(item_data$PRICE_ORIGINAL)
item_data$OFFER_VISUALISATION <- as.factor(item_data$OFFER_VISUALISATION)
item_data$OFFER_POSITION <- as.factor(item_data$OFFER_POSITION)
item_data$COUNTRY_NAME <- as.factor(item_data$COUNTRY_NAME)
item_data$Year <- as.factor(item_data$Year)
item_data$Quarter <- as.factor(item_data$Quarter)
item_data$maxpers <- as.factor(item_data$maxpers)

#Drop the NA's ASSUMPTION, MOET NOG NAAR GEKEKEN WORDEN
item_data <- drop_na(item_data)

#Create the dummies 
item_final <- fastDummies::dummy_cols(item_data[,-1], remove_selected_columns = TRUE)
item_final$MailOffer <- item_data$MailOffer
item_final <- item_final %>% select("MailOffer", everything())

#save dataset
saveRDS(item_final, "item_final.Rds")
#Clear environment
rm(list = ls())



