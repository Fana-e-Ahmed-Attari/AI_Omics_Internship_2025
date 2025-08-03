install.packages("tidyverse")   # For data manipulation & cleaning
install.packages("readr")       # For reading CSV files
install.packages("here")        # For managing file paths

# Set Working Directory 
setwd("D:/R_Projects/AI_Omics_Internship_2025/Module_I")

# Create folders
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

# Load library
library(tidyverse)

# Load dataset
patient_data <- read_csv("patient_info.csv")

# View first few rows
head(patient_data)

# Structure of dataset
str(patient_data)

# Summary stats
summary(patient_data)

# Check column types
glimpse(patient_data)

# Convert age from character to numeric
patient_data$age <- as.numeric(patient_data$age)

# Assuming the column is 'smoker'
patient_data <- patient_data %>%
  mutate(smoker_binary = ifelse(smoker == "Yes", 1, 0))

# Check missing values
colSums(is.na(patient_data))

# Remove rows with too many NAs or impute
patient_data <- drop_na(patient_data)  # Or use fill()

write_csv(patient_data, "patient_info_clean.csv")

save.image("FanaeAhmed_Class_Ib_Assignment.RData")
