library(readxl)
library(writexl)
library(gdata)
library(dplyr)

data <- df.crea_2
data.pat <- df.pat.sel_2
data.gwas <- read.table("List_patients_STCS_with_GWAS.txt")

ids <- unique(data$patid)

change <- list()

for(i in 1:length(ids)){ #this nested loop removes ids with < 3 measurements, discards ids with baseline < 90, calculates % change in egfr from baseline if baseline >= 90

  sub.data <- as.data.frame(subset(data, patid == ids[i])) #subset by id

  oldest.date <- min(sub.data$crea.date.new) #take earliest date available
  diff.check <- difftime(time2 = oldest.date , time1 = sub.data$crea.date.new, units = "days")

  if(length(sub.data[,1]) < 3 | sum(diff.check >= 30) == 0){ #remove ids with 1 and 2 measurements

    sub.data$change <- NA

  } else if(sum(sub.data$egfr[sub.data$crea.date.new == oldest.date] >= 90) > 0){ #if baseline >= 90, then calculate change in egfr for the other values

    sub.data$change <- (sub.data$egfr - 90)/90*100 #drop from 90

  } else { #if baseline < 90, then NA

    sub.data$change <- NA

  }

  change[[i]] <- sub.data

}

data.final <- do.call(change, what = rbind)
data.final <- data.final[!is.na(data.final$change),]

#CASES----

#cases1st <- data.final[data.final$change >= 25,] # take only data greater than threshold

cs.sum <- group_by(data.final, patid) %>% summarize(n = length(change), pct = sum(change <= -25)/length(change)*100)
cs.sum <- cs.sum[cs.sum$pct >= 30,] #taking patients with at least 30% measures >= 25% drop from 90

id.cs <- cs.sum$patid

cases <- data.final[data.final$patid %in% id.cs,]

#cs.sum <- group_by(cases, organ) %>% summarize(n = length(unique(patid)))

par(mfrow=c(1,2))

plot(x = cases$crea.date.new, y = cases$patid, type = "n", ylab = "Patient ID", xlab = "Year")

for(i in 1:length(id.cs)){

  sub.data <- subset(cases, patid == id.cs[i])
  col <- ifelse(sub.data$change <= -25, "red", ifelse(sub.data$change <= -15 & sub.data$change > -25, "orange", "green"))

  abline(h = sub.data$patid, col = "grey", lwd = 0.2)
  points(sub.data$crea.date.new, sub.data$patid, col = col, pch = 19, cex = 0.5)


}


#CONTROLS----

data.nocases <- data.final[!data.final$patid %in% id.cs,] #remove cases from data.final

id.ct <- unique(data.nocases$patid) #ids of controls (no incidence density sampling)

plot(x = data.nocases$crea.date.new, y = data.nocases$patid, type = "n", ylab = "Patient ID", xlab = "Year")

for(i in 1:length(id.ct)){

  sub.data <- subset(data.nocases, patid == id.ct[i])
  col <- ifelse(sub.data$change <= -25, "red", ifelse(sub.data$change <= -15 & sub.data$change > -25, "orange", "green"))

  abline(h = sub.data$patid, col = "grey", lwd = 0.2)
  points(sub.data$crea.date.new, sub.data$patid, col = col, pch = 19, cex = 0.5)


}

#INCIDENCE DENSITY SAMPLING----

set.seed(123456)

id.ct.ids <- rep(NA, length(id.cs)) #vector filled with NA of length of cases

potential.ct.ids <- list()

for(i in 1:length(id.cs)){ #iterations over cases

  cs.sub <- data.final[data.final$patid %in% id.cs[i],] #subset specific case

  cs.organ <- df.pat.sel_2[df.pat.sel_2$patid == id.cs[i],]$tpx

  tpx.d.cs <- unique(cs.sub$tpxdate) #tpx date of the case

  cs.date1 <- cs.sub$crea.date.new[cs.sub$crea.date.new == min(cs.sub$crea.date.new)] #baseline date
  cs.date2 <- cs.sub$crea.date.new[cs.sub$change <= -25] #dates with drop >= 25%
  cs.date2 <- cs.date2[cs.date2 == min(cs.date2)] #take earliest date available with drop >= 25%

  followup.cs <- as.vector(difftime(time2 = tpx.d.cs, time1 = cs.date2, units = "days"))/365 #followup period tpx date - earliest date available with drop >= 25%

  id.ct.corr <- id.ct[!id.ct %in% id.ct.ids[!is.na(id.ct.ids)]] #removes ct that have been already selected from the ct pool

  ids.ct <- rep(NA, length(id.ct.corr)) #empty matrix 2 columns and n rows = n cts available

  for(j in 1:length(id.ct.corr)){ #loop over all available cts

    #print(j)

    ct.sub <- data.final[data.final$patid %in% id.ct.corr[j],] #subset specific control

    ct.organ <- df.pat.sel_2[df.pat.sel_2$patid == id.ct.corr[j],]$tpx

    tpx.d.ct <- unique(ct.sub$tpxdate) #tpx date of the control

    tpx.diff <- as.vector(difftime(time2 = tpx.d.cs,
                                   time1 = tpx.d.ct,
                                   units = "days"))/365 #time difference between cs and ct tpx dates

    followups.ct <- as.vector(difftime(time2 = tpx.d.ct,
                                       time1 = ct.sub$crea.date.new,
                                       units = "days"))/365 #time difference between ct tpx date and all ct creatinine dates

    followup.diff <- abs(followup.cs - followups.ct) #time difference btw ct followup and cs followups

    if(min(followup.diff)*365 >= 180){ #if min followup diff >= 30 days go to next control

      next

    }

    followup.ct.date <- ct.sub$crea.date.new[followup.diff == min(followup.diff)] #take minimum time difference as ct followup upr limit


    followup.check <- ct.sub$change[ct.sub$crea.date.new <= followup.ct.date] >= -15 #check if all dates within ct followup have egfr >= -15

    after.fup.date <- min(ct.sub$crea.date.new[ct.sub$crea.date.new > followup.ct.date])


    dt1.log <- abs(tpx.diff) <= 1 #TRUE if tpx ct is +- 1 from tpx cs
    green.followup <- sum(followup.check) == length(followup.check) #TRUE if followup ct is green
    followup.next.check <- ct.sub$change[ct.sub$crea.date.new == after.fup.date] >= -15

    sex.log <- unique(cs.sub$sex) == unique(ct.sub$sex) #TRUE if ct and cs of same sex
    organ.log <- cs.organ == ct.organ #TRUE if ct and cs same organ

    if(sum(c(dt1.log, green.followup, followup.next.check, sex.log, organ.log)) == 5){ #if same organ & sex, followup.diff <= 1, ct has measurement greater than irst date of cs and within +-1 y from cs baseline --> then  take cs id together with followup.diff

      ids.ct[j] <- unique(ct.sub$patid) #because there are multiple lines for each patient

    } else { #if the ct does not satisfy the above criteria --> next ct

      next

    }

  }

  if(sum(!is.na(ids.ct)) == 0){ #if no ct is selected for the specific case --> next case

    next

  } else if(sum(!is.na(ids.ct)) == 1){

    id.ct.ids[i] <- na.omit(ids.ct)

  } else {

    id.ct.ids[i] <- sample(na.omit(ids.ct), size = 1)

  }

  potential.ct.ids[[i]] <- ids.ct

}

na.omit(id.ct.ids)

pat.sel.ids <- data.frame(cases = id.cs, controls = id.ct.ids)
usethis::use_data(pat.sel.ids, overwrite = T)

#CHECKING IDS WITH GWAS ALREADY----

sum(id.cs %in% data.gwas[,1]) #number of cases with already gwas data

sum(na.omit(id.ct.ids) %in% data.gwas[,1]) #number of controls with already gwas data



















