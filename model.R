##setwd("N:/BALTCOAST DATA/")
##setwd("C:/MV/GeoPop/geopop/")

## PLOT_TO_FILE <- TRUE
## ABOVE <- TRUE
## FOLDER <- "."
## if(!FALSE) {
##     SPECIES <- "Pleuronectes platessa"
##     MINDSTEMAAL <- 26
## }
## if(FALSE) {
##     SPECIES <- "Gadus morhua"
##     MINDSTEMAAL <- 30
## }

## install.packages("devtools")
## devtools::install_github("kaskr/gridConstruct",subdir="gridConstruct")
## install.packages("TMB")
## install.packages("sp")

library(DATRAS)
library(gridConstruct)
load("d.tot.RData2",verbose =TRUE)
## Build grid
##d <- subset(d,lon!=0 & lon<30 & lat>-38 & lat>-10-1.4*lon)
d.tot <- subset(d.tot, lon<13)

gr <- gridConstruct(d.tot[["HH"]],
                    km=10,
                    scale=1.2,
                    center=data.frame(lon=11.5,lat=56),
                    nearestObs=40)

##plot(gr)

## Attach grid cell ID to each haulid:
##d.tot[[2]]$gf <- gridFactor(d.tot,gr)
d.tot$gf <- gridFactor(d.tot, gr)

## Change invalid TMB name:
##d.tot[[2]]$haulId <- d.tot[[2]]$haul.id
d.tot$haulId <- d.tot$haul.id

## Log( haul_duration ) standardized to 30 minutes
d.tot$logHaulDur <- log( d.tot$HaulDur / 30 )

## Select species
d.tot <- subset(d.tot, Species == SPECIES)

## Response variables
group <- ABOVE + 1  ## = 1 or 2
d.tot <- addSpectrum(d.tot, cm.breaks = c(0, MINDSTEMAAL, Inf))
response <- d.tot$N[,group]
d.tot[[2]]$response <- response

## Output filename
FILENAME <- paste0(FOLDER,"/",SPECIES," ",colnames(d.tot$N)[group])

## Factor defining our time unit
time <- factor(d.tot$Year)
##d.tot[[2]]$time <- time
d.tot$time <- time

## FIXME:
d.tot$Depth[is.na(d.tot$Depth)] <- mean(d.tot$Depth, na.rm=TRUE)

## Sparse matrices for GMRF: Q = Q0+delta*I
Q0 <- -attr(gr,"pattern")
diag(Q0) <- 0
diag(Q0) <- -rowSums(Q0)
I <- .symDiagonal(nrow(Q0))

## Load TMB
library(TMB)
##compile("model.cpp","-O0 -g")
compile("model.cpp")
dyn.load(dynlib("model"))


## TMB Data
data <- list(
  Q0=Q0,
  I=I,
  time=d.tot[[2]]$time ,
  gf=d.tot[[2]]$gf  ,
  haulId=d.tot[[2]]$haulId  ,
  response=d.tot[[2]]$response,
  gear = d.tot[[2]]$Gear,
  depth = d.tot[[2]]$Depth,
  logHaulDur = d.tot[[2]]$logHaulDur,
  ## Asbjoern:
  tbott = d.tot[[2]]$tbott,
  sbott = d.tot[[2]]$sbott,
  tsurf = d.tot[[2]]$tsurf,
  ssurf = d.tot[[2]]$ssurf,
  LevitusT = d.tot[[2]]$LevitusT,
  LevitusR = d.tot[[2]]$LevitusR,
  stratidx = d.tot[[2]]$stratidx,
  oxy_max= d.tot[[2]]$oxy_max,
  oxy_bot= d.tot[[2]]$oxy_bot
)
data0 <- data[sapply(data,length) == length(data$gf)] ## Covariates for model matrix
## How independent are AB's predictors ?:
## qw <- data[c("tbott","sbott","ssurf","LevitusT","LevitusR","stratidx",
##              "oxy_max","oxy_bot")]
## cor(as.data.frame(qw))

fit_model <- function(data, with_static_field=FALSE) {
    ## Stuff for prediction: Dummy dataset that matches the space time grid:
    ## Xpredict: design matrix for prediction
    ## Apredict: Area factor for prediction
    DFpredict <- expand.grid(gf=levels(data$gf), time=levels(data$time))
    ## FIXME: We should include depth and covariates here !
    ##        But that requires depth on the entire grid...
    Xpredict <- model.matrix(~time, data=DFpredict) ## <-- gear removed
    stopifnot( all( colnames(Xpredict) %in% colnames(data$X) ) ) ## Validity check
    tmp <- matrix(0, nrow(Xpredict), ncol(data$X))
    tmp[,match(colnames(Xpredict), colnames(data$X))] <- Xpredict
    data$Xpredict <- tmp
    data$Apredict <- factor(icesSquare(gr))
    data$doPredict <- 1
    ## Perhaps we want measure all indices relative to a fixed reference square:
    ##   plot(d.tot, plot.response = FALSE)
    ## "42G1" seems appropriate
    ## data$refindex <- which(levels(data$Apredict) == "42G1")
    data$refindex <- 0 ## <-- Disable
    
    parameters <- list(
        eta_density = matrix(0,nrow(Q0),nlevels(time)) ,
        eta_nugget = rep(0, nrow(d.tot[[2]])),
        logdelta=-4.85,
        logscale=0.4,
        logsd_nugget=0.22,
        time_corr=2.57,
        beta = rep(0, ncol(data$X))
    )
    parameters$eta_static <- rep(0, nlevels(data$gf) * with_static_field )
    parameters$logdelta_static <- rep(0, 1 * with_static_field )
    parameters$logscale_static <- rep(0, 1 * with_static_field )
    obj <- MakeADFun(data, parameters, random=c("eta_density","eta_nugget","eta_static"),
                     profile = "beta")
    fit <- nlminb(obj$par,obj$fn,obj$gr)
    rep <- obj$report(obj$env$last.par.best)
    rownames(rep$logindex) <- levels(data$Apredict)
    colnames(rep$logindex) <- levels(data$time)
    sdr <- sdreport(obj)
    s <- summary(sdr)
    est <- s[rownames(s) != "eta_density" & rownames(s) != "eta_nugget",]
    s <- summary(sdr,p.value=TRUE)
    s1 <- s[rownames(s)=="beta",]
    rownames(s1)<-colnames(data$X)
    s1
    s2 <- s[rownames(s)=="eta_density",]
    return(environment())
}


#################### Model configurations
## Base model (no effects)
X0 <- model.matrix(~time + gear + logHaulDur,
                   data = data0)
## Include depth
X1 <- model.matrix(~time + gear + logHaulDur + poly(depth,2),
                   data = data0)
## Include asbjoern
X2 <- model.matrix(~time + gear + logHaulDur + poly(depth,2) +
                       tbott + sbott + ssurf +
                       LevitusT + LevitusR + stratidx +
                       oxy_max + oxy_bot,
                   data = data0)
## We can add a 2nd order polynomial for each period if we like:
data0$time_cut <- data0$time
levels(data0$time_cut) <-
    as.character(cut(as.numeric(levels(data0$time)),4))
X3 <- model.matrix(~time + gear + time_cut:poly(depth,2),
                   data = data0)
if(FALSE){
    X4 <- model.matrix(~time + gear + I(factor(time%in%1991:2007)):poly(depth,2),data=data[c("time","gear","depth")])
}

################### Fit models
data$X <- X0
env0 <- fit_model(data)

data$X <- X1
env1 <- fit_model(data)

data$X <- X2
env2 <- fit_model(data)

## data$X <- X3
## env3 <- fit_model(data)

if(PLOT_TO_FILE) pdf(paste0(FILENAME,".pdf"))

## Full area overview
layout(1)
plot(d.tot, plot.response = !TRUE, col="red", pch=16, cex=.5)

## Legend pos
legendPos <- ifelse(
    mean(diff(colMeans(env0$rep$logindex))) > 0 ,
    "topleft",
    "bottomleft"
)

with(env0,
{
    ## Plot results
    layout(1)
    ## Grand mean
    matplot(as.numeric(levels(data$time)),
            as.matrix(colMeans(rep$logindex)),
            type="l",xlab="Year",ylab="log(index)")
    legend(legendPos,legend="Grand mean",col=1,lty=1)
})

## Sub Area overview
my_squares <- outer(43:41, c("G0","G1","G2"), paste0)

with(env0,
{
    ## Plot results
    layout(matrix(1:4,2))
    plot(subset(d.tot, (StatRec %in% my_squares) & lat<57.5),
         plot.response = !TRUE, col="red", pch=16, cex=.5)
    ## West -> Center -> East
    for(i in 1:3) {
        matplot(as.numeric(levels(data$time)),
                t(rep$logindex[my_squares[,i],]),
                type="l",xlab="Year",ylab="log(index)")
        legend(legendPos,legend=my_squares[,i],col=1:3,lty=1:3)
    }
})

source("utils.R")
modelList <- list(model2=env2, model1=env1, model0=env0)
tab <- do.call("LRtest",modelList)
getExtraOutput <- function(x) {
    c(SpatialSD = x$rep$SpatialSD,
      NuggetSD = as.numeric(exp(x$fit$par["logsd_nugget"])),
      AvgSDlogindex = mean(with(x, est[rownames(est)=="logindex",2])),
      ForecastSD = x$rep$ForecastSD
      )
}
tab <- cbind(tab,
             t(sapply(modelList, getExtraOutput)))

## Color maps - entire area
heatpal <- colorRampPalette(c("white", "yellow","orange" ,"red"))
pl <- as.list(env0$sdr,"Estimate")
layout(1)
for(i in 1:nlevels(d.tot$time)){
  plot(gr,type="n",xlab="Longitude",ylab="Latitude",las=1)
  image(gr, concTransform(pl$eta_density[,i]),
        col=heatpal(18)[-(1:2)],
        map=quote(map("worldHires",add=TRUE, fill=TRUE, col="grey")),
        add=TRUE)
  title(levels(d.tot$time)[i])
}

if(PLOT_TO_FILE) {
    dev.off()
    options(width=100); sink(paste0(FILENAME,".txt"));
    print(tab);
    cat("\n")
    print(env2$s1)
    sink()
}
