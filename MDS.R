########################################################################################
################################## function for MDS  ###################################
########################################################################################
#scientific notation
options(scipen=0)
options(digits=22)

library(fpc)
library(clValid)

#this function is called only after the arrival
MDS <- function(k){
  if(k<2){
    cat("Number of cluster must be equal or greater than 2!!\n")
    cat("Current best k: ",2,"\n")
    return()
  }else{
    #fit classical MDS module
    fit <- cmdscale(d=sim.matrix,eig=TRUE, k=k)
    #vectorize and list data, cluster fit, used k, internal measures and 
    #stability measures
    internal <- as.numeric(attributes(clValid(fit$points,nClust = k, validation="internal"))$measures)
    stability <- as.numeric(attributes(clValid(fit$points,nClust = k, validation="stability"))$measures)
    names(internal) <- c("con","DI","silhue")
    names(stability) <- c("APN","AD","ADM","FOM")
    res <- list(k=k,int=internal,stab=stability)
    return(res)
  }
}

#function to optimize k (number of clusters in MDS)
optimize.k <- function(k.old,increase=TRUE){
  #first clusterization input
  if(is.null(increase)){
    #call MDS function with new k - 1, then k + 1
    k.minus.1 <- MDS(k=k.old-1)
    k.plus.1 <- MDS(k=k.old+1)
    k.curr <- MDS(k=k.old)
    
    #Internal measures
    #Connectivity
    aux.minus.con <- ifelse(k.curr$int["con"] > k.minus.1$int["con"], TRUE, FALSE)
    aux.plus.con <- ifelse(k.curr$int["con"] > k.plus.1$int["con"], TRUE, FALSE)
    
    #silhuette
    aux.minus.silhue <- ifelse(k.curr$int["silhue"] < k.minus.1$int["silhue"], TRUE, FALSE)
    aux.plus.silhue <- ifelse(k.curr$int["silhue"] < k.plus.1$int["silhue"], TRUE, FALSE)
    
    #DI
    aux.minus.DI <- ifelse(k.curr$int["DI"] < k.minus.1$int["DI"], TRUE, FALSE)
    aux.plus.DI <- ifelse(k.curr$int["DI"] < k.plus.1$int["DI"], TRUE, FALSE)
    
    minus <- length(which(c(aux.minus.con,aux.minus.silhue,aux.minus.DI)==TRUE))
    plus <- length(which(c(aux.plus.con,aux.plus.silhue,aux.plus.DI)==TRUE))
    
    curr.stab.sum <- sum(k.curr$stab["APN"],k.curr$stab["AD"],k.curr$stab["ADM"],k.curr$stab["FOM"])
    minus.stab.sum <- sum(k.minus.1$stab["APN"],k.minus.1$stab["AD"],k.minus.1$stab["ADM"],k.minus.1$stab["FOM"])
    plus.stab.sum <- sum(k.plus.1$stab["APN"],k.plus.1$stab["AD"],k.plus.1$stab["ADM"],k.plus.1$stab["FOM"])
    
    
    
    #if no change in current k is better
    if(minus == plus & (minus == 0 || minus == 1)){
      res.k <- k.old
      cat("Current best k: ",res.k,"\n")
      return(res.k)
    }
   
    if(minus > plus){#if k-1 better than current k
      res.k <- k.old - 1
      #continue decreasing until local minimum
      optimize.k(k.old=k.old - 1,increase = FALSE)
    }else if(minus < plus){#if k+1 better than current k
      res.k <- k.old + 1
      #continue increasing until local minimum
      optimize.k(k.old=k.old + 1,increase = TRUE)  
    }else{
      res.k <- k.old
      #stop in local minimum
      cat("Current best k: ",res.k,"\n")    
    }
    
    #if internal measures are equally better
    #use best global stability to decide to increase or decrease 
    if(minus == plus & (minus > 1)){
      
      #stability
      if(minus.stab.sum < plus.stab.sum){
        res.k <- k.old - 1
        #stop in local minimum
        cat("Current best k: ",res.k,"\n")
        return(res.k)
      }else{
        res.k <- k.old + 1
        #stop in local minimum
        cat("Current best k: ",res.k,"\n")
        return(res.k)
      }
    }
    
    
    
  }else if(increase==TRUE){
    k.plus.1 <- MDS(k=k.old)
    k.curr <- MDS(k=k.old-1)
    
    #browser()
    #Internal measures
    #Connectivity
    aux.plus.con <- ifelse(k.curr$int["con"] > k.plus.1$int["con"], TRUE, FALSE)
    
    #silhuette
    aux.plus.silhue <- ifelse(k.curr$int["silhue"] < k.plus.1$int["silhue"], TRUE, FALSE)
    
    #DI
    aux.plus.DI <- ifelse(k.curr$int["DI"] < k.plus.1$int["DI"], TRUE, FALSE)
    
    plus <- length(which(c(aux.plus.con,aux.plus.silhue,aux.plus.DI)==TRUE))
    
    curr.stab.sum <- sum(k.curr$stab["APN"],k.curr$stab["AD"],k.curr$stab["ADM"],k.curr$stab["FOM"])
    plus.stab.sum <- sum(k.plus.1$stab["APN"],k.plus.1$stab["AD"],k.plus.1$stab["ADM"],k.plus.1$stab["FOM"])
    
    if(plus > 1){
      res.k <- k.old
      #continue increasing until local minimum
      optimize.k(k.old=k.old + 1,increase = TRUE)
    }else{
      #stop in current k
      res.k <- k.old - 1
      cat("Current best k: ",res.k,"\n")
      return(res.k)
    }
    
  }else if(increase==FALSE){
    k.minus.1 <- MDS(k=k.old)
    k.curr <- MDS(k=k.old+1)
    
    #browser()
    
    #Internal measures
    #Connectivity
    aux.minus.con <- ifelse(k.curr$int["con"] > k.minus.1$int["con"], TRUE, FALSE)
    
    #silhuette
    aux.minus.silhue <- ifelse(k.curr$int["silhue"] < k.minus.1$int["silhue"], TRUE, FALSE)
    
    #DI
    aux.minus.DI <- ifelse(k.curr$int["DI"] < k.minus.1$int["DI"], TRUE, FALSE)
    
    minus <- length(which(c(aux.minus.con,aux.minus.silhue,aux.minus.DI)==TRUE))
    
    curr.stab.sum <- sum(k.curr$stab["APN"],k.curr$stab["AD"],k.curr$stab["ADM"],k.curr$stab["FOM"])
    minus.stab.sum <- sum(k.minus.1$stab["APN"],k.minus.1$stab["AD"],k.minus.1$stab["ADM"],k.minus.1$stab["FOM"])
    
    if(minus > 1){
      res.k <- k.old
      #continue decreasing until local minimum
      optimize.k(k.old=k.old - 1,increase = FALSE)
    }else{
      #stop in current k
      res.k <- k.old + 1
      cat("Current best k: ",res.k,"\n")
      return(res.k)
    }
    
    
  }
  
}

########################################################################################
############################ function for MDS  - END ###################################
########################################################################################

#example
data(mouse)
#distance matrix of numerical data
sim.matrix <<- dist(mouse[,2:7])
#starting number of clusters
curr.k = 6
cat("Starting k: ",curr.k,"\n")
#optimize MDS (for k parameter)
optimize.k(k.old=curr.k,increase = NULL)
cat("Ended simulation Iterations \n")

cat("\n")
#2nd example
#starting number of clusters
curr.k = 8
cat("Starting k: ",curr.k,"\n")
#optimize MDS (for k parameter)
optimize.k(k.old=curr.k,increase = NULL)
cat("Ended simulation Iterations \n")

cat("\n")
#3rd example
#starting number of clusters
curr.k = 4
cat("Starting k: ",curr.k,"\n")
#optimize MDS (for k parameter)
optimize.k(k.old=curr.k,increase = NULL)
cat("Ended simulation Iterations \n")