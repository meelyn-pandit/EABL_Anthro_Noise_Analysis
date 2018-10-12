###############RFID Data Analysis 2018################

library(feedr)
library(tidyverse)

##Loading raw RFID data##
data = load_raw("control_rfid_logs.TXT")
attach(data)

#df <- data[-grep(c(1:9), data$animal_id),]
####Housekeeping before using feedr####
tag_reads = data %>% count(animal_id)
write.table(tag_reads, file = "tag_reads.csv", append = TRUE, col.names = TRUE)

##Housekeeping##
animal_index = read.csv("animal_index.csv") #importing list of all the birds banded with RFID tags (incomplete list as of now)
head(animal_index)
names(animal_index) = c("animal_id", "species","age", "sex") #renaming the column headers so it includes at least "animal_id" and "species"

clean_data = check_ids(data, ids = animal_index) #checking the tags read by the reader against the list of the known RFID tags. You will get tags that do not show up in animal_index. You will need to enter these tag IDs into the animal_index csv file and label them as error or wand.

v <- visits(clean_data, bw = 3) #between reads is 5 min. will need to discuss this. should try to base time between visits on previous bluebird data.
ARF2_F_v = filter(v,animal_id == "01103F3A27")
head(v)

library(dplyr)
visit_summary <- v %>%
  group_by(animal_id) %>%
  summarize(n_visits = length(start),
            n_loggers = length(unique(logger_id)),
            n_days = length(unique(date)),
            mean_visit_length = mean(as.numeric(difftime(end, start, units = "sec"))))
visit_summary
visit_summary$mean_visit_rate = (visit_summary$n_visits/visit_summary$n_days)
visit_summary$treatment = rep("control",nrow(visit_summary))
write.table(visit_summary, file = "rfid_visit_summaries.csv", append = TRUE, sep = ",", col.names = FALSE)

library(ggplot2) #Displaying average number of visits
ggplot(data = visit_summary, aes(x = animal_id, y = mean_visit_length)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #rotate/shift id labels
  geom_bar(stat = "identity") + # Create bars
  labs(x = "Animal ID", y = "Average visit length (s)") # Make labels pretty

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
