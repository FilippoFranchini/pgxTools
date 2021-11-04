#Matching birthdays with ids & only dna available

library(readxl)
library(writexl)

scqm.insel.pat <- read_excel("Insel SCQM Patients 2021-07-14.xlsx")
scqm.pat.bdates <- read_excel("Insel SCQM bdates.xlsx", col_names = F)
colnames(scqm.pat.bdates) <- c("SCQM ID", "bday")

scqm.pat.dat <- merge(x = scqm.insel.pat, y = scqm.pat.bdates, by = "SCQM ID")

scqm.pat.dat <- subset(scqm.pat.dat, dna_sample_available == "yes" & !Diagnose == "UA")

write_xlsx(scqm.pat.dat, path = paste0("SCQM/","scqm_insel_pat.xlsx"))




# Select patients with PsA (M07), RA (M06), axSpA (M45) ----

diagnoses <- read.table("diagnoses.csv", header = T, sep = ";")

icd10.split <- strsplit(diagnoses$icd10_code, split = ".", fixed = T)
diagnoses$ICD10 <- unlist(lapply(icd10.split, function(x) x[1]))

myICD10 <- c("M07", "M05", "M06", "M45")

diagnoses <- diagnoses[diagnoses$ICD10 %in% myICD10,]

patid.icd10 <- unique(diagnoses$pseudo_pid)


#Merge general with lab----

lab <- read.table("lab_values.csv", header = T, sep = ";")
lab <- lab[,c(1,7,9)]

general <- read.table("general_data.csv", header = T, sep = ";")

for(i in 1:length(lab[,1])){

  lab$pseudo_scqm_id[i] <- general[general$pseudo_pid %in% lab$pseudo_pid[i],]$pseudo_scqm_id

  lab$sex[i] <- general[general$pseudo_pid %in% lab$pseudo_pid[i],]$sex

  lab$year_of_birth[i] <- general[general$pseudo_pid %in% lab$pseudo_pid[i],]$year_of_birth

  lab$nsaid[i] <- general[general$pseudo_pid %in% lab$pseudo_pid[i],]$nsaid

}

date.time.split <- strsplit(lab$date_request, split = " ", fixed = T)
lab$date_request <- unlist(lapply(date.time.split, function(x) x[1]))

lab$date_request <- as.Date(lab$date_request, "%Y-%m-%d")
lab$year_of_birth <- as.Date(as.character(lab$year_of_birth), "%Y")

colnames(lab) <- c("pseudo_pid", "creatinin", "creatinindate", "pseudo_scqm_id", "sex", "birthday", "nsaid")


lab$creatininAge <- as.vector(difftime(time2 = lab$birthday,
                                        time1 = lab$creatinindate,
                                        units = "days"))/365

lab$sex <- as.factor(lab$sex)
levels(lab$sex) <- c("female", "male")

lab$creatinin_c <- lab$creatinin/10^6*113.12*10^3/10 #creatinin in mg/dl

egfrs <- c()

for(i in 1:nrow(lab)){

  egfrs[i] <- gfr_calc(scr = lab$creatinin_c[i],
                       age = lab$creatininAge[i],
                       gender = lab$sex[i])

}

lab$egfr <- egfrs





ids <- unique(lab$pseudo_pid)

change <- list()

for(i in 1:length(ids)){

  sub.data <- subset(lab, pseudo_pid == ids[i])

  oldest.date <- min(sub.data$creatinindate)

  if(sum(sub.data$egfr[sub.data$creatinindate == oldest.date] >= 90) > 0){

    sub.data$change <- abs(sub.data$egfr - sub.data$egfr[sub.data$creatinindate == oldest.date])/sub.data$egfr[sub.data$creatinindate == oldest.date]*100

  } else {

    sub.data$change <- NA

  }

  change[[i]] <- sub.data

}

data.final <- do.call(change, what = rbind)
data.final <- na.omit(data.final)

cases1st <- data.final[data.final$change >= 25,] # take only data greater than threshold

cs.sum <- group_by(cases1st, pseudo_pid) %>% summarize(n = length(change))
cs.sum <- cs.sum[cs.sum$n > 1,] #taking patients with more than 1 measurement for confirmation period

id.cs <- cs.sum$pseudo_pid

cases2nd <- data.final[data.final$pseudo_pid %in% id.cs,]

####
cases2nd <- cases2nd[cases2nd$pseudo_pid %in% patid.icd10,]

cases2nd <- cases2nd[cases2nd$pseudo_pid %in% patid.med,]

id.cs <- unique(cases2nd$pseudo_pid)
####





plot(x = cases2nd$creatinindate, y = cases2nd$pseudo_pid, type = "n")

for(i in 1:length(id.cs)){

  sub.data <- subset(cases2nd, pseudo_pid == id.cs[i])
  col <- ifelse(sub.data$change >= 25, "red", "green")

  abline(h = sub.data$pseudo_pid, col = "grey")
  points(sub.data$creatinindate, sub.data$pseudo_pid, col = col, pch = 19, cex = 0.5)


}










# Medication----

#Selection for methotrexate (ATC L01BA01)

med <- read.table("medications.csv", header = T, sep = ";")
map <- read.table("IDSC202101727_data_v2_keys_20210915.csv", header = T, sep = ";")




methotrexate <- c("L01BA01", "L04AX03")

med <- med[med$atc_code %in% methotrexate,]

patid.med <- med$pseudo_pid









