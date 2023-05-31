library(MASS)
library(survival)
library(TMB)
library(magic)
library(tidyverse)

#load cpp
#compile("C:/Users/ktang3/Desktop/Imperial/SSA_mort/TMB/nbinom_timesp_split_crossage_RW_dispmf_crosstime_avg_sex_rwanda.cpp")
dyn.load(dynlib("C:/Users/ktang3/Desktop/Imperial/SSA_mort/TMB/nbinom_timesp_split_crossage_RW_dispmf_crosstime_avg_sex_rwanda"))

pyears_data <- readRDS("data/pyears_data_smoothed.rds")

#basis####
age.start <- 10; age.end <- 65
year.start <- 1990; year.end <- 2020

aggr.mat.cohort.0 <- filter(pyears_data, period >= year.start) %>%
  group_by(country) %>% group_split %>%
  setNames(pyears_data$country %>% levels)

joint.countries <- names(aggr.mat.cohort.0)

knots.time <- seq(-3 * 2.5 + year.start, year.end + 3 * 2.5, by = 2.5)
no.basis.time <- length(knots.time) - 4
knots.age <- seq(-3 * 2.5 + age.start, age.end + 3 * 2.5, by = 2.5)
no.basis.age <- length(knots.age) - 4

bspline<-function (x,k,i,m=2) {
  if (m==-1) {basis<-as.numeric(x<k[i+1] & x>=k[i])} else {
    z0<-(x-k[i])/(k[i+m+1]-k[i])
    z1<-(k[i+m+2]-x)/(k[i+m+2]-k[i+1])
    basis<-z0*bspline(x,k,i,m-1)+z1*bspline(x,k,i+1,m-1) }
  basis
}

A.age<-c()
for(j in 1:no.basis.age) {
  A.age<-cbind(A.age, bspline(age.start:age.end, knots.age, j))
}

A.year<-c()
for(j in 1:no.basis.time) {
  A.year<-cbind(A.year,bspline(year.start:year.end, knots.time, j))
}

te.spline <- A.year %x% A.age


#geojson
library(geojsonR); library(sp); library(geojsonio); library(rgeos); library(stringr)
bound <- geojson_read("custom.geo.json", what = "sp") 

cen <- coords <- list(); cen.names <- c()
for(i in 1:length(bound)){
  if(bound[i,]$name %in% joint.countries){
    cen[[length(cen)+1]] = c(gCentroid(bound[i,])$x, gCentroid(bound[i,])$y)
    coords <- append(coords, list(bound[i,]@polygons[[1]]@Polygons[[1]]@coords))
    cen.names <- c(cen.names, bound[i,]$name)
  }
}

names(coords) <- names(cen) <- cen.names
cen <- c(cen, list("Congo Democratic Republic" = c(gCentroid(bound[which(bound$name == "Dem. Rep. Congo"),])$x, gCentroid(bound[which(bound$name == "Dem. Rep. Congo"),])$y)))
cen <- c(cen, list("Cote d'Ivoire" = c(gCentroid(bound[which(bound$name == "Côte d'Ivoire"),])$x, gCentroid(bound[which(bound$name == "Côte d'Ivoire"),])$y)))
cen <- c(cen, list("Eswatini" = c(gCentroid(bound[which(bound$name == "Swaziland"),])$x, gCentroid(bound[which(bound$name == "Swaziland"),])$y)))
cen <- c(cen, list("Central African Republic" = c(gCentroid(bound[which(bound$name == "Central African Rep."),])$x, gCentroid(bound[which(bound$name == "Central African Rep."),])$y)))

coords[["Congo Democratic Republic"]] <- bound[which(bound$name == "Dem. Rep. Congo"),]@polygons[[1]]@Polygons[[1]]@coords
coords[["Cote d'Ivoire"]] <- bound[which(bound$name == "Côte d'Ivoire"),]@polygons[[1]]@Polygons[[1]]@coords
coords[["Eswatini"]] <- bound[which(bound$name == "Swaziland"),]@polygons[[1]]@Polygons[[1]]@coords
coords[["Central African Republic"]] <- bound[which(bound$name == "Central African Rep."),]@polygons[[1]]@Polygons[[1]]@coords

cen <- cen[order(names(cen))]
coords <- coords[order(names(coords))]

coords.neighbors <- setNames(lapply(joint.countries, function(x){
  ph <- c()
  if(!is.null(coords[[x]])){
    for(j in joint.countries[!joint.countries %in% x]){
      if(!is.null(coords[[j]])) {
        #if(!all(is.na(match(data.frame(t(coords[[x]])),data.frame(t(coords[[j]])))))) {
        if(min(apply(coords[[x]], 1, function(i){apply(coords[[j]], 1, function(y){sum(sqrt((y-i)^2))})})) < 1e-3){
          ph <- c(ph, j)
        }
      }
    }
    ph
  }}), joint.countries)

coords.neighbors$Lesotho <- c("South Africa")
coords.neighbors$"South Africa" <- c(coords.neighbors$`South Africa`, "Lesotho")

no.polygon <- length(which(sapply(coords.neighbors, is.null)))

dist.from.south.africa <- sapply(cen, function(x){sqrt((x[1] - cen$`South Africa`[1]) ^ 2 + (x[2] - cen$`South Africa`[2]) ^ 2)})
dist.from.south.africa["Eswatini"] <- 1
dist.from.south.africa[names(which(sapply(coords.neighbors, is.null)))] <- Inf
dist.order <- order(dist.from.south.africa[joint.countries])

#create data and paras for TMB ####
dm.all <- df.all <- Em.all <- Ef.all <- c()
DX.spline.m.all <- DX.spline.f.all <-
  DX.spline.age.m.all <- DX.spline.age.f.all <- 
  DX.spline.time.m.all <- DX.spline.time.f.all <- 
  tp.mat.m.all <- tp.mat.f.all <- 
  intercept.mat.m.all <- intercept.mat.f.all <-NA

penal.age.all <- penal.time.all <- halfpenal.time.all1 <- halfpenal.time.all2 <- 
  penal.age.2d.all <- penal.time.2d.all <- penal.age.time.2d.all <- NA

null.age <- null.time <- null <- NA

avg.m.all <- avg.f.all <- tp.params.all <- 
  spline.params.age.m.all <- spline.params.age.f.all <- 
  spline.params.time.m.all <- spline.params.time.f.all <- 
  spline.params.2d.m.all <- spline.params.2d.f.all <- c()

tp.mat.m.all.RW <- tp.mat.f.all.RW <- NA
tp.mat.m.all.common.RW <- tp.mat.f.all.common.RW <- c()
col.split <- 5
rwanda.mat.m <- rwanda.mat.f <- t(rep(NA, 3))

for (i in 1:length(joint.countries)){
  cat(i, "\n")
  aggr.mat.reduced <- aggr.mat.cohort.0[[joint.countries[i]]]
  
  data.mat.m <- aggr.mat.reduced %>%
    filter(mm1 == "male",
           agegr %in% age.start:age.end,
           period %in% year.start:year.end) %>%
    dplyr::select(adjusted, pyears2, agegr, period, tips)
  
  data.mat.f <- aggr.mat.reduced %>%
    filter(mm1 == "female",
           agegr %in% age.start:age.end,
           period %in% year.start:year.end) %>%
    dplyr::select(adjusted, pyears2, agegr, period, tips)
  
  DX.spline.m <- te.spline[(age.end - age.start + 1) * (data.mat.m$period - year.start) + data.mat.m$agegr - age.start + 1,]
  DX.spline.f <- te.spline[(age.end - age.start + 1) * (data.mat.f$period - year.start) + data.mat.f$agegr - age.start + 1,]
  
  DX.spline.age.m <- A.age[data.mat.m$agegr - age.start + 1,]
  DX.spline.age.f <- A.age[data.mat.f$agegr - age.start+1,]
  DX.spline.time.m <- A.year[data.mat.m$period - year.start+1,]
  DX.spline.time.f <- A.year[data.mat.f$period - year.start+1,]
  
  P.age <- diff(diag(no.basis.age), differences = 2)
  P.time <- diff(diag(no.basis.time), differences = 2)
  P.age.2d <- diag(no.basis.time) %x% crossprod(P.age)
  P.time.2d <- crossprod(P.time) %x% diag(no.basis.age)
  P.age.time.2d <- crossprod(diff(diag(no.basis.time)) %x% diff(diag(no.basis.age)))
  #  C.constraints<-rbind(diag(no.basis)%x%t(apply(A.age,2,sum)),t(apply(A.year,2,sum))%x%diag(no.basis))
  
  C.constraints<-rbind(diag(no.basis.time) %x% t(apply(A.age[min(data.mat.m$agegr):max(data.mat.m$agegr) - age.start + 1,], 2, sum)),
                       t(apply(A.year[min(data.mat.m$period):max(data.mat.m$period) - year.start + 1,], 2, sum)) %x% diag(no.basis.age))
  
  dm.all <- c(dm.all, data.mat.m$adjusted); df.all <- c(df.all, data.mat.f$adjusted)
  Em.all <- c(Em.all, data.mat.m$pyears2); Ef.all <- c(Ef.all, data.mat.f$pyears2)
  
  DX.spline.m.all <- as(adiag(as.matrix(DX.spline.m.all), DX.spline.m), "sparseMatrix")
  DX.spline.f.all <- as(adiag(as.matrix(DX.spline.f.all), DX.spline.f), "sparseMatrix")
  DX.spline.age.m.all <- as(adiag(as.matrix(DX.spline.age.m.all), DX.spline.age.m), "sparseMatrix")
  DX.spline.age.f.all <- as(adiag(as.matrix(DX.spline.age.f.all), DX.spline.age.f), "sparseMatrix")
  DX.spline.time.m.all <- as(adiag(as.matrix(DX.spline.time.m.all), DX.spline.time.m), "sparseMatrix")
  DX.spline.time.f.all <- as(adiag(as.matrix(DX.spline.time.f.all), DX.spline.time.f), "sparseMatrix")
  
  penal.age.all <- as(adiag(as.matrix(penal.age.all), crossprod(P.age)), "sparseMatrix")
  penal.time.all <- as(adiag(as.matrix(penal.time.all), crossprod(P.time)), "sparseMatrix")
  penal.age.2d.all <- as(adiag(as.matrix(penal.age.2d.all), P.age.2d), "sparseMatrix")
  penal.time.2d.all <- as(adiag(as.matrix(penal.time.2d.all), P.time.2d), "sparseMatrix")
  penal.age.time.2d.all <- as(adiag(as.matrix(penal.age.time.2d.all), P.age.time.2d), "sparseMatrix")
  halfpenal.time.all1 <- as(adiag(as.matrix(halfpenal.time.all1), P.time[1:(col.split - 1),]), "sparseMatrix")
  halfpenal.time.all2 <- as(adiag(as.matrix(halfpenal.time.all2), P.time[-(1:(col.split - 1)),]),"sparseMatrix")
  #  null.age<-as(adiag(as.matrix(null.age),1e-3*diag(ncol(P.age))+exp(15)*tcrossprod(apply(A.age,2,sum))),"sparseMatrix")
  #  null.time<-as(adiag(as.matrix(null.time),1e-3*diag(ncol(P.time))+exp(15)*tcrossprod(apply(A.year,2,sum))),"sparseMatrix")
  
  null.age <- as(adiag(as.matrix(null.age), 1e-3 * diag(ncol(P.age)) + exp(15) * tcrossprod(apply(A.age[min(data.mat.m$agegr):max(data.mat.m$agegr) - age.start + 1,], 2, sum))), "sparseMatrix")
  null.time <- as(adiag(as.matrix(null.time), 1e-3 * diag(ncol(P.time)) + exp(15) * tcrossprod(apply(A.year[min(data.mat.m$period):max(data.mat.m$period) - year.start + 1,], 2, sum))), "sparseMatrix")
  null <- as(adiag(as.matrix(null), 1e-3 * diag(ncol(P.age.2d)) + exp(15) * crossprod(C.constraints)), "sparseMatrix")
  
  intercept.mat.m.all <- as(adiag(as.matrix(intercept.mat.m.all), cbind(rep(1, nrow(DX.spline.age.m)))), "sparseMatrix")
  intercept.mat.f.all <- as(adiag(as.matrix(intercept.mat.f.all), cbind(rep(1, nrow(DX.spline.age.f)))), "sparseMatrix")
  tp.mat.m.all.RW <- adiag(tp.mat.m.all.RW, model.matrix( ~ factor(data.mat.m$tips, levels = 0:14) - 1))
  tp.mat.f.all.RW <- adiag(tp.mat.f.all.RW, model.matrix( ~ factor(data.mat.f$tips, levels = 0:14) - 1))
  tp.mat.m.all.common.RW <- rbind(tp.mat.m.all.common.RW, model.matrix( ~ factor(data.mat.m$tips, levels = 0:14) - 1))
  tp.mat.f.all.common.RW <- rbind(tp.mat.f.all.common.RW, model.matrix( ~ factor(data.mat.f$tips, levels = 0:14) - 1))
  
  if(joint.countries[i] == "Rwanda"){
    rwanda.mat.m <- as(rbind(as.matrix(rwanda.mat.m), cbind(as.numeric(data.mat.m$period == 1993),
                                                            as.numeric(data.mat.m$period == 1994),
                                                            as.numeric(data.mat.m$period == 1995))), "sparseMatrix")
    rwanda.mat.f <- as(rbind(as.matrix(rwanda.mat.f), cbind(as.numeric(data.mat.f$period == 1993),
                                                            as.numeric(data.mat.f$period == 1994),
                                                            as.numeric(data.mat.f$period == 1995))), "sparseMatrix")
  } else{
    rwanda.mat.m <- as(rbind(as.matrix(rwanda.mat.m), matrix(0, nrow(data.mat.m), 3)), "sparseMatrix")
    rwanda.mat.f <- as(rbind(as.matrix(rwanda.mat.f), matrix(0, nrow(data.mat.f), 3)), "sparseMatrix")
  }
}
