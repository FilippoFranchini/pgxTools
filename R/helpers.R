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


stcs_select <- function(datastring, GWAS.data, eGFR.limit = 90, l.cs = 25,
                        l.ct = 15, tol = 1){

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
                                          time1 = data$creatinindate,
                                          units = "days"))/365

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

    if(sum(sub.data$egfr[sub.data$assperiod == 0] >= eGFR.limit) > 0){

      sub.data$change <- abs(sub.data$egfr - sub.data$egfr[sub.data$assperiod == 0])/sub.data$egfr[sub.data$assperiod == 0]*100

    } else {

      sub.data$change <- NA

    }

    change[[i]] <- sub.data

  }

  data.final <- do.call(change, what = rbind)
  data.final <- na.omit(data.final) # remove patients with baseline < 90


  #cases
  cases1st <- data.final[data.final$change >= l.cs,] # take only data greater than threshold

  cs.sum <- group_by(cases1st, organ, patid) %>% summarize(n = length(change))
  cs.sum <- cs.sum[cs.sum$n > 1,] #taking patients with more than 1 measurement for confirmation period

  id.cs <- cs.sum$patid

  cases2nd <- data.final[data.final$patid %in% id.cs,] # patients must have at least two measurements with change >= threshold

  cs.sum <- group_by(cases2nd, organ) %>% summarize(n = length(unique(patid)))


  #controls
  data.nocases <- data.final[!data.final$patid %in% id.cs,]

  controls1st <- data.nocases[data.nocases$change <= l.ct,] #1st ct selection: change less than threshold

  ct.sum <- group_by(controls1st, organ, patid) %>% summarize(n = length(change))
  ct.sum <- ct.sum[ct.sum$n > 1,]

  controls2nd <- data.nocases[data.nocases$patid %in% ct.sum$patid,] #2nd ct selection: more than two measurement must be below threshold

  id.ct <- ct.sum$patid


  #incidence density sampling

  id.ct.ids <- list()

  for(i in 1:length(id.cs)){

    cs.sub <- data.final[data.final$patid %in% id.cs[i],]
    cs.sub <- cs.sub[order(cs.sub$creatinindate),] #order from oldest to newest

    cs.date1 <- cs.sub$creatinindate[cs.sub$assperiod == 0] #baseline
    cs.date2 <- cs.sub$creatinindate[cs.sub$change >= l.cs][1] #first available measure

    ids.ct <- c()

    for(j in 1:length(id.ct)){

      ct.sub <- data.final[data.final$patid %in% id.ct[j],]

      dt1 <- as.vector(difftime(time2 = cs.date1,
                                time1 = ct.sub$creatinindate, units = "days"))/365

      dt2 <- as.vector(difftime(time2 = cs.date2,
                                time1 = ct.sub$creatinindate, units = "days"))/365

      dt1.log <- sum(dt1 <= 1) #+- tollerance at baseline
      dt2.log <- sum(dt2 >= 0) #last ct date must be > first date of cs

      if(dt1.log >= 1 & dt2.log >= 1){

        ids.ct[j] <- unique(ct.sub$patid) #because there are multiple lines for each patient

      } else {

        ids.ct[j] <- NA

      }

    }


    id.ct.ids[[i]] <- na.omit(ids.ct)

  }


  id.ct <- unique(unlist(id.ct.ids))

  controls <- controls2nd[controls2nd$patid %in% id.ct,]

  ct.sum <- group_by(controls, organ) %>% summarize(n = length(unique(patid)))


  #GWAS comparison

  GWAS.ids <- unlist(read.table(GWAS.data, header = F))
  GWAS.n <- sum(c(id.cs, id.ct) %in% GWAS.ids)


  #plots
  #dev.new()

  #par(mfrow=c(1,2))

  #hist(cases$egfr, breaks = 40,
  #     xlab = expression(paste("eGFR (mL/min/",m^2,")")),
  #     main = "CS eGFR distribution")

  #hist(controls$egfr, breaks = 40,
  #     xlab = expression(paste("eGFR (mL/min/",m^2,")")),
  #     main = "CT eGFR distribution")

  #dev.new()

  #par(mfrow=c(2,1))

  #boxplot(cases$change ~ cases$patid, ylim = c(0, 100))
  #boxplot(controls$change ~ controls$patid, ylim = c(0, 100))

  #print(paste("Controls contain", paste0(sum(id.cs %in% id.ct)), " cases."))

  return(list(tab1 = data.frame(cs.sum),
              tab2 = data.frame(ct.sum),
              ids = list(cases = id.cs, controls = id.ct),
              GWAS.n = GWAS.n
              )
         )

}
