###############RFID Data Analysis 2018################

library(feedr)
library(tidyverse)
library(dplyr)
library(compare)

setwd("~/OneDrive - University of Oklahoma/University of Oklahoma/Bridge Lab/EABL Anthro Noise/data_analysis")

##Loading raw RFID data##
data = load_raw("ARF_04_072018.TXT")
data = data[(which(nchar(data$animal_id) == 10)),] #removing tag reads that are less than 10 characters long
attach(data)

#df <- data[-grep(c(1:9), data$animal_id),]
####Housekeeping before using feedr####
tag_reads = data %>% count(animal_id)
write.table(tag_reads, file = "arf_4_tag_reads.csv", append = FALSE, col.names = TRUE, sep = ",") #Manually created a "problems" csv in Excel that has original tag read and corrected read

##Housekeeping##
animal_index = read.csv("ok_int_noise_animal_index.csv", header = TRUE) #importing list of all the birds banded with RFID tags (incomplete list as of now)
head(animal_index)
names(animal_index) = c("animal_id", "species","age", "sex") #renaming the column headers so it includes at least "animal_id" and "species"

comparison = semi_join(tag_reads,animal_index) #comparing the tag_reads and animal_index vectors to see which tags are found in each vector, that way other tagged EABLs do not get counted as an error when we do the error check

animal_index = read.csv("ok_int_noise_animal_index.csv", header = TRUE) #loading animal_index dataframe that now has additional error values
animal_index$animal_id = as.character(animal_index$animal_id)
animal_index = animal_index[(which(nchar(animal_index$animal_id) == 10)),] #removing tag reads that are less than 10 characters long

error_clean = check_ids(data, ids = animal_index) #checking the tags read by the reader against the list of the known RFID tags. You will get tags that do not show up in animal_index. You will need to enter these tag IDs into the animal_index csv file and label them as error or wand.

problems <- read.csv("arf_4_problems.csv")
clean_data <- check_problems(error_clean, problems = problems) #taking misreads and correcting them to the actual PIT tag id
clean_data$animal_id = as.character(clean_data$animal_id) #removing tag reads that are less than 10 characters long

v <- visits(error_clean, bw = 2, na_rm=TRUE) #between reads is 2 s. will need to discuss this. should try to base time between visits on previous bluebird data.

visit_summary <- v %>%
  group_by(animal_id) %>%
  summarize(n_visits = length(start)/2,
            n_loggers = length(unique(logger_id)),
            n_days = length(unique(date)),
            n_hour = length(unique(hour)),
            feedr_mean_visit_length = mean(as.numeric(difftime(end, start, units = "sec"))))

visit_summary$mean_visit_rate = (visit_summary$n_visits/visit_summary$n_days)
visit_summary$treatment = rep("int_noise",nrow(visit_summary))
visit_summary$box = "arf_4"
visit_summary

write.table(visit_summary, file = "arf_4_visit_summaries.csv", append = FALSE, sep = ",", col.names = TRUE)
write.table(visit_summary, file = "ok_int_visit_summaries.csv", sep = ",", append = TRUE, col.names = FALSE)

###########################################################################################################
######## Mixed Model analysis################################################################
################################################################################################
setwd("~/OneDrive - University of Oklahoma/University of Oklahoma/Bridge Lab/EABL Anthro Noise/data_analysis")

library(nlme)
library(lme4)
library(lattice)
library(bbmle)
library(ggplot2)
library(ordinal)
library(ggthemes)
library(plyr)

rfid_visit_summaries <- read_csv("~/OneDrive - University of Oklahoma/University of Oklahoma/Bridge Lab/EABL Anthro Noise/data_analysis/rfid_visit_summaries.csv", col_types = cols(julian_band_date = col_date(format = "%m/%d/%Y")))
View(rfid_visit_summaries)
attach(rfid_visit_summaries)
adult_visits = filter(rfid_visit_summaries, age != "hy")
adult_visits = filter(adult_visits, age != "n/a")
detach(rfid_visit_summaries)
attach(adult_visits)
adult_visits = na.omit(adult_visits)

#plotting for distributions
hist(mean_visit_rate)
shapiro.test(mean_visit_rate)
#poisson.test(mean_visit_rate)

#Log transforming data to make it more normal, p-value is still below 0.05...
adult_visits$mean_visit_rate = adult_visits$mean_visit_rate/adult_visits$num_nest
adult_visits$log_mean_visit_rate = log(adult_visits$mean_visit_rate)
attach(adult_visits)
shapiro.test(log_mean_visit_rate)
hist(log_mean_visit_rate)


#testing for correlations among variables
cor(log_mean_visit_rate)

####Models For Visit Rate####
m1 = lmer(log_mean_visit_rate ~ 1 + (1|band_num), data = adult_visits)
m2 = lmer(log_mean_visit_rate ~ treatment + (1|band_num), data = adult_visits)
m3 = lmer(log_mean_visit_rate ~ treatment * sex + (1|band_num), data = adult_visits)
m4 = lmer(log_mean_visit_rate ~ treatment * age + (1|band_num), data = adult_visits)
m5 = lmer(log_mean_visit_rate ~ julian_band_date + (1|band_num), data = adult_visits)
m6 = lmer(log_mean_visit_rate ~ treatment + julian_band_date + (1|band_num), data = adult_visits)

AICctab(m1,m2,m3,m4,m5,m6, nobs = 37, base=TRUE,delta=TRUE, sort=TRUE, weights=TRUE)

summary(m2)
coefs <- data.frame(coef(summary(m2))) 
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs



#Displaying average number of visits - Scatterplot
sp = ggplot(adult_visits, aes(x=treatment, y=log_mean_visit_rate, color = age))
sp  + 
  geom_point(size=5, shape=19) + 
  labs(x="Treatment", y = "Visit Rate \n(# of visits/day)") +
  theme_few(base_size = 20)  + 
  theme(axis.text=element_text(size=28), axis.title=element_text(size=30,face="bold")) 


#Graphing visitation rate - box plot
ggplot(adult_visits, aes(x = treatment, y = log_mean_visit_rate)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #rotate/shift id labels
  geom_bar(stat = "identity") + # Create bars
  theme_few(base_size = 20)  + 
  geom_errorbar(aes(ymin=log_mean_visit_rate-sd(log_mean_visit_rate), ymax=log_mean_visit_rate+sd(log_mean_visit_rate)), width=.2,
                position=position_dodge(.9)) +
  labs(x = "Treatment Group", y = "Average visit rate (visits/day)") # Make labels pretty

#######Measuring Growth Rates#######



# ###Filtering visit data by tag read so we can calculate visit length correctly
# SVS2_M_v = filter(v,animal_id == "011016E6D7")
# SVS2_M_v$visit_length = as.numeric(SVS2_M_v$end) - as.numeric(lag(SVS2_M_v$start))
# mean(SVS2_M_v$visit_length[2:531],)
# 
# SVS2_F_v = filter(v,animal_id == "011016C39D")
# SVS2_F_v$visit_length = as.numeric(SVS2_F_v$end) - as.numeric(lag(SVS2_F_v$start))
# mean(SVS2_F_v$visit_length[2:3084],)
# 
# ARF2_F_v = filter(v,animal_id == "01103F3A27")
# ARF2_F_v$visit_length = as.numeric(ARF2_F_v$end) - as.numeric(lag(ARF2_F_v$start))
# mean(ARF2_F_v$visit_length[2:3084],)
#   
# write.table(ARF2_F_v, file = "01103F3A27_reads.csv", append = FALSE, sep = ",", col.names = TRUE)



################################################################################################


p <- presence(v, bw = 5) #Gives duration of time inside nest box, time bewteen visits is 5 min
ARF2_F_p = filter(p,animal_id == "01103F3A27")

a <- activity(p, res = 2)
ARF2_F_a = filter(a,animal_id == "01103F3A27")

i <- convert_activity(clean_data)
library(activity)

# Calculate daily activity pattern for a single individual
a <- fitact(i[[1]], sample = "none")
plot(a)
