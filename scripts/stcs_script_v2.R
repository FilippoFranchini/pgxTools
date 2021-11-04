library(readxl)
library(writexl)
library(gdata)
library(dplyr)

data <- read.xls("PGX_input_29SEP21.xlsx", header = T, sheet = 1) #
data.gwas <- read.table("List_patients_STCS_with_GWAS.txt")

data <- data[!is.na(data$crea),] #removing NA creatinin NA rows

# Date formatting
data$tpxdate <- as.Date(data$tpxdate, "%Y-%m-%d")
data$birthday <- as.Date(data$birthday, "%Y-%m-%d")
data$creatinindate <- as.Date(data$creatinindate, "%Y-%m-%d")

mydata <- data[,c(1:4,8:9,11:12)] # Removing useless columns

mydata <- mydata[mydata$crea >= 0,] #remove negative values of crea

# Calculate age of the patients at creatinindate
mydata$creatininage <- as.vector(difftime(time2 = mydata$birthday,
                                        time1 = mydata$creatinindate,
                                        units = "days"))/365

mydata$sex <- as.factor(mydata$sex) #gender as factor

mydata$creatinin_c <- mydata$crea/10^6*113.12*10^3/10 #creatinin in mg/dl

#calculate EGFR and add it to mydata
egfrs <- c()

for(i in 1:length(mydata[,1])){

  egfrs[i] <- gfr_calc(scr = mydata$creatinin_c[i],
                       age = mydata$creatininage[i],
                       gender = mydata$sex[i])

}

mydata$egfr <- egfrs

mydata <- mydata[!mydata$organ == "Kidney",] #remove kidney transplant

ids <- unique(mydata$patid) #patients' ids

change <- list()

for(i in 1:length(ids)){ #this nested loop removes ids with < 3 measurements, discards ids with baseline < 90, calculates % change in egfr from baseline if baseline >= 90

  sub.data <- subset(mydata, patid == ids[i]) #subset by id

  oldest.date <- min(sub.data$creatinindate) #take earliest date available

  if(length(sub.data[,1]) < 3){ #remove ids with 1 and 2 measurements

    sub.data$change <- NA

  } else if(sum(sub.data$egfr[sub.data$creatinindate == oldest.date] >= 90) > 0){ #if baseline >= 90, then calculate change in egfr for the other values

    sub.data$change <- abs(sub.data$egfr - sub.data$egfr[sub.data$creatinindate == oldest.date])/sub.data$egfr[sub.data$creatinindate == oldest.date]*100

  } else { #if baseline < 90, then NA

    sub.data$change <- NA

  }

  change[[i]] <- sub.data

}

data.final <- do.call(change, what = rbind)
data.final <- na.omit(data.final)

#CASES----

#cases1st <- data.final[data.final$change >= 25,] # take only data greater than threshold

cs.sum <- group_by(data.final, patid) %>% summarize(n = length(change), pct = sum(change >= 25)/length(change)*100)
cs.sum <- cs.sum[cs.sum$pct >= 30,] #taking patients with at least 30% measures >= 25% drop from 90

id.cs <- cs.sum$patid

cases <- data.final[data.final$patid %in% id.cs,]

#cs.sum <- group_by(cases, organ) %>% summarize(n = length(unique(patid)))


plot(x = cases$creatinindate, y = cases$patid, type = "n")

for(i in 1:length(id.cs)){

  sub.data <- subset(cases, patid == id.cs[i])
  col <- ifelse(sub.data$change >= 25, "red", ifelse(sub.data$change >= 15, "orange", "green"))

  abline(h = sub.data$patid, col = "grey", lwd = 0.2)
  points(sub.data$creatinindate, sub.data$patid, col = col, pch = 19, cex = 0.5)


}


#CONTROLS----

data.nocases <- data.final[!data.final$patid %in% id.cs,] #remove cases from data.final

id.ct <- unique(data.nocases$patid) #ids of controls (no incidence density sampling)

plot(x = data.nocases$creatinindate, y = data.nocases$patid, type = "n")

for(i in 1:length(id.ct)){

  sub.data <- subset(data.nocases, patid == id.ct[i])
  col <- ifelse(sub.data$change >= 25, "red", ifelse(sub.data$change >= 15, "orange", "green"))

  abline(h = sub.data$patid, col = "grey", lwd = 0.2)
  points(sub.data$creatinindate, sub.data$patid, col = col, pch = 19, cex = 0.5)


}

#INCIDENCE DENSITY SAMPLING----


id.ct.ids <- rep(NA, length(id.cs)) #vector filled with NA of length of cases

for(i in 1:length(id.cs)){ #iterations over cases

  cs.sub <- data.final[data.final$patid %in% id.cs[i],] #subset specific case

  tpx.d.cs <- unique(cs.sub$tpxdate) #tpx date of the case

  cs.date1 <- cs.sub$creatinindate[cs.sub$creatinindate == min(cs.sub$creatinindate)] #baseline date
  cs.date2 <- cs.sub$creatinindate[cs.sub$change >= 25] #dates with drop >= 25%
  cs.date2 <- cs.date2[cs.date2 == min(cs.date2)] #take earliest date available with drop >= 25%

  followup.cs <- as.vector(difftime(time2 = tpx.d.cs, time1 = cs.date2, units = "days"))/365 #followup period tpx date - earliest date available with drop >= 25%

  id.ct.corr <- id.ct[!id.ct %in% id.ct.ids[!is.na(id.ct.ids)]] #removes ct that have been already selected from the ct pool

  ids.ct <- matrix(nrow = length(id.ct.corr), ncol = 2) #empty matrix 2 columns and n rows = n cts available

  for(j in 1:length(id.ct.corr)){ #loop over all available cts

    ct.sub <- data.final[data.final$patid %in% id.ct.corr[j],] #subset specific control

    tpx.d.ct <- unique(ct.sub$tpxdate) #tpx date of the control

    dt1 <- as.vector(difftime(time2 = cs.date1,
                              time1 = ct.sub$creatinindate, units = "days"))/365 #time difference between cs baseline and all dates of the ct

    dt2 <- as.vector(difftime(time2 = cs.date2,
                              time1 = ct.sub$creatinindate, units = "days"))/365 #time difference between earliest cs data with >= 25% drop and all dates of the ct

    first.date <- min(ct.sub$creatinindate[dt2 >= 0]) #earliest ct date bigger than earliest cs data with >= 25%

    followup.ct <- as.vector(difftime(time2 = tpx.d.ct, time1 = first.date, units = "days"))/365 #followup period tpx date - earliest ct date bigger than earliest cs data with >= 25%

    followup.diff <- unique(followup.ct - followup.cs) #difference between cs and ct followups

    dt1.log <- sum(abs(dt1) <= 1) #+- tollerance at baseline (>= 1 if abs(dt1) contains values <= 1 year)
    dt2.log <- sum(dt2 >= 0) #last ct date must be > first date of cs (>= 1 if dt2 contains values >= 0 years)
    followup.log <- abs(followup.diff) <= 1 #followup difference 1 year (TRUE if followup.diff contains values <= 1 years)
    sex.log <- unique(cs.sub$sex) == unique(ct.sub$sex) #TRUE if ct and cs of same sex
    organ.log <- unique(cs.sub$organ) == unique(ct.sub$organ) #TRUE if ct and cs same organ

    if(dt1.log > 0 & dt2.log > 0 & followup.log & sex.log & organ.log){ #if same organ & sex, followup.diff <= 1, ct has measurement greater than irst date of cs and within +-1 y from cs baseline --> then  take cs id together with followup.diff

      ids.ct[j,1] <- unique(ct.sub$patid) #because there are multiple lines for each patient
      ids.ct[j,2] <- followup.diff

    } else { #if the ct does not satisfy the above criteria --> next ct

      next

    }

  }

  if(sum(!is.na(ids.ct)) == 0){ #if no ct is selected for the specific case --> next case

    next

  }

  followup.min <- ids.ct[abs(ids.ct[,2]) %in% min(abs(ids.ct[,2]), na.rm = T),][1] #if multiple cts are selected for the case --> take the ct with lowes followup.diff

  id.ct.ids[i] <- followup.min #store ct id --> next case

}

na.omit(id.ct.ids)

#CHECKING IDS WITH GWAS ALREADY----

sum(id.cs %in% data.gwas[,1]) #number of cases with already gwas data

sum(na.omit(id.ct.ids) %in% data.gwas[,1]) #number of controls with already gwas data



















