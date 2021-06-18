gfr_calc <- function(scr, age, gender){

  if(gender == "male"){

    kappa <- 0.9
    alpha <- -0.411
    beta <- 1

  } else if(gender == "female"){

    kappa <- 0.7
    alpha <- -0.329
    beta <- 1.018

  }

  eGFR <- 141*(min(scr/kappa, 1)^alpha)*(max(scr/kappa, 1)^-1.209)*(0.993^age)*beta

  return(eGFR)

}


stcs_select <- function(datastring, eGFR.limit){

  data <- read.xls(datastring, header = T, sheet = 1)

  data <- data[!is.na(data$creatinin),] #removing NA creatinin NA rows

  # Date formatting
  data$tpxdate <- as.Date(data$tpxdate, "%Y-%m-%d")
  data$birthday <- as.Date(data$birthday, "%Y-%m-%d")
  data$creatinindate <- as.Date(data$creatinindate, "%Y-%m-%d")

  # Removing useless columns
  data <- data[,c(3:7,11:12)]

  # Calculating age @ creatinin and eGFR

  data$creatininAge <- as.vector(difftime(time2 = data$birthday,
                                          time1 = data$creatinindate))/365

  data$sex <- as.factor(data$sex)
  levels(data$sex) <- c("female", "male")

  data$creatinin_c <- data$creatinin/10^6*113.12*10^3/10 #creatinin in mg/dl

  egfrs <- c()

  for(i in 1:length(data[,1])){

    egfrs[i] <- gfr_calc(scr = data$creatinin_c[i],
                         age = data$creatininAge[i],
                         gender = data$sex[i])

  }

  data$egfr <- egfrs

  # Calculating % change from baseline > 90 mL/min/1.73m2

  ids <- unique(data$patid)

  change <- list()

  for(i in 1:length(ids)){

    sub.data <- subset(data, patid == ids[i])

    if(sub.data$egfr[sub.data$assperiod == 0] >= eGFR.limit){

      sub.data$change <- abs(sub.data$egfr - sub.data$egfr[sub.data$assperiod == 0])/sub.data$egfr[sub.data$assperiod == 0]*100

    } else {

      sub.data$change <- NA

    }

    change[[i]] <- sub.data

  }

  data.final <- do.call(change,what = rbind)
  data.final <- na.omit(data.final) # remove patients with baseline < 90
  data.final <- data.final[!data.final$change == 0,] # remove baselines

  data.final <- data.final[data.final$change >= 25,] # take only data with change >= 25%

  sum.tab1 <- group_by(data.final, organ, patid) %>% summarize(n = length(change))

  sum.tab2 <- group_by(sum.tab, organ) %>% summarize(npat = length(patid))

  dev.new()
  hist(data$egfr, breaks = 20,
       xlab = expression(paste("eGFR (mL/min/",m^2,")")),
       main = "eGFR distribution")

  return(list(tab1 = data.frame(sum.tab1), tab2 = data.frame(sum.tab2)))


}
