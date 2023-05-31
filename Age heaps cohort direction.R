library(rdhs)
library(demogsurv)
library(data.table)
library(survival)

pyears_data <- readRDS("C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/data/pyears_data.rds")
#sib_data <- readRDS("C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/data/sib_data.rds")

kernel.0 <- function(x, smooth.width = 2, dist.width = 2){
  cat(as.character(x$survey_id[1]), " ", as.character(x$mm1[1]), "\n")
  if(!is.numeric(x$tips)) {x$tips <- as.numeric(levels(x$tips))[x$tips]}
  y <- smooth.width + 1
  w <- dist.width
  x$pyears2 <- x$pyears
  x$adjusted <- x$expected <- x$weights <- x$event
  heap.ages <- seq(5, 65, by = 5)
  ind<-which(x$agegr %in% heap.ages)
  
  for(p in (1:length(x$event))[-ind]){
    tar.x<-x[p,]
    cohort.neighbors.ind<-c(p)
    for(p2 in 1:smooth.width){
      cohort.neighbors.ind<-if(!(tar.x$agegr - p2 %in% heap.ages)) c(which(x$mm1==tar.x$mm1 & x$period==tar.x$period-p2 & x$tips==as.numeric(tar.x$tips)+p2 & x$agegr==tar.x$agegr-p2),cohort.neighbors.ind)
      cohort.neighbors.ind<-if(!(tar.x$agegr + p2 %in% heap.ages)) c(cohort.neighbors.ind,which(x$mm1==tar.x$mm1 & x$period==tar.x$period+p2 & x$tips==as.numeric(tar.x$tips)-p2 & x$agegr==tar.x$agegr+p2))
    }
    cohort.x<-x[cohort.neighbors.ind,]
    weights<-(1/y-abs(cohort.x$agegr-tar.x$agegr)/(y^2))
    x$weights[p]<-sum(cohort.x$event*(weights/sum(weights)))
  }
  
  for(i in ind){
    tar.x<-x[i,]
    cohort.neighbors.ind<-c()
    dist.ind<-c(i)
    for(j in 1:smooth.width){
      cohort.neighbors.ind<-c(which(x$mm1==tar.x$mm1 & x$period==tar.x$period-j & x$tips==as.numeric(tar.x$tips)+j & x$agegr==tar.x$agegr-j),cohort.neighbors.ind)
      cohort.neighbors.ind<-c(cohort.neighbors.ind,which(x$mm1==tar.x$mm1 & x$period==tar.x$period+j & x$tips==as.numeric(tar.x$tips)-j & x$agegr==tar.x$agegr+j))
    }
    
    cohort.x<-x[cohort.neighbors.ind,]
    weights<-(1/y-abs(cohort.x$agegr-tar.x$agegr)/(y^2))
    x$weights[i]<-x$expected[i]<-sum(cohort.x$event*(weights/sum(weights)))
    
    for(k in 1:dist.width){
      dist.ind<-c(which(x$mm1==tar.x$mm1 & x$period==tar.x$period-k & x$tips==as.numeric(tar.x$tips)+k & x$agegr==tar.x$agegr-k),dist.ind)
      dist.ind<-c(dist.ind,which(x$mm1==tar.x$mm1 & x$period==tar.x$period+k & x$tips==as.numeric(tar.x$tips)-k & x$agegr==tar.x$agegr+k))
    }
    
    if(x$event[i]>x$expected[i]){
      if(all(x$expected[dist.ind]==0)){x$adjusted[dist.ind]<-x$expected[dist.ind]+(x$event[i]-x$expected[i])/(2*w+1)} else {
        x$adjusted[dist.ind]<-x$expected[dist.ind]+(x$event[i]-x$expected[i])*x$weights[dist.ind]/sum(x$weights[dist.ind])
      }
      
      for(q in 1:length(dist.ind)){
        if(dist.ind[q]>i){
          x$pyears2[dist.ind[q]]<-x$pyears[dist.ind[q]]+sum(x$adjusted[dist.ind[q:length(dist.ind)]]-x$expected[dist.ind[q:length(dist.ind)]])/x$event[i]*x$pyears[i]
        } else{
          x$pyears2[dist.ind[q]]<-x$pyears[dist.ind[q]]-(sum(x$adjusted[dist.ind]-x$expected[dist.ind])-sum(x$adjusted[dist.ind[q:length(dist.ind)]]-x$expected[dist.ind[q:length(dist.ind)]]))/x$event[i]*x$pyears[i]
        }
      }
    }
  }
  return(x)
}

pyears.list <- pyears_data %>%
  mutate(agegr = as.numeric(levels(agegr)[agegr]),
         period = as.numeric(levels(period)[period])) %>%
  rename(mm1 = sibsex) %>%
  group_by(survey_id, survey_type, mm1) %>%
  group_split()

pyears.cohort.smooth.list <- lapply(pyears.list, kernel.0, smooth.width = 2, dist.width = 2)

aggr.mat.cohort.0 <- pyears.cohort.smooth.list %>% bind_rows()

saveRDS(aggr.mat.cohort.0, "C:/Users/ktang3/Desktop/Imperial/SSA_mort/JE/data/pyears_data_smoothed.rds")
