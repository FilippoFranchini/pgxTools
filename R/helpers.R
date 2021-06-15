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
