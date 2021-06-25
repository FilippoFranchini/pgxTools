# Methotrexate only

data.med <- d.med_final[d.med_final$m.medication_generic_drug == "methotrexate" & is.na(d.med_final$discontinuation_reason),]

id.med <- sort(unique(data.med$patient_id))

data.pat <- d.pat_final[d.pat_final$patient_id %in% id.med,]
data.pat <- data.pat[data.pat$diagnose == "PsA",] #methotrexate + PsA only

id.pat <- sort(unique(data.pat$patient_id))

data.vis <- d.vis_final[d.vis_final$patient_id %in% id.pat,]
data.vis <- data.vis[is.na(data.vis$v.present_medication_arthritis_conventional_nsaid),]

for (i in 1:length(data.vis[,1])){

  data.vis$gender[i] <- data.pat$p.gender[data.pat$patient_id %in% data.vis$patient_id[i]]

}

data.vis$kreatinine_c <- data.vis$v.kreatinin_value/10^6*113.12*10^3/10


eGFRs <- c()

for (i in 1:length(data.vis[,1])){

  eGFRs[i] <- gfr_calc(scr = data.vis$kreatinine_c[i], age = as.numeric(data.vis$v.visit_date_age[i], units="days")/365, gender = data.vis$gender[i])

}

hist(eGFRs, breaks = 40)

data.vis$eGFR <- eGFRs

data.vis <- data.vis[,c(1,6:9)]

data.vis <- na.omit(data.vis)


library(dplyr)

n <- data.frame(group_by(data.vis, patient_id) %>% summarise(c = length(gender)))

id.final <- n[n$c > 2,]$patient_id


data.vis <- data.vis[data.vis$patient_id %in% id.final,]


change <- list()

for(i in 1:length(id.final)){

  sub.data <- subset(data.vis, patient_id == id.final[i])

  if(sub.data$eGFR[1] >= 90){

    sub.data$change <- abs(sub.data$eGFR - sub.data$eGFR[1])/sub.data$eGFR[1]*100

  } else {

    sub.data$change <- NA

  }

  change[[i]] <- sub.data

}


data.final <- do.call(change,what = rbind)
data.final <- na.omit(data.final)
data.final <- data.final[!data.final$change == 0,]

cases <- data.final[data.final$change >= 25,]

cs.sum <- group_by(cases, patient_id) %>% summarise(n = length(change))

cs.sum <- cs.sum[cs.sum$n > 1,]

boxplot(change~patient_id,data=data.final)
