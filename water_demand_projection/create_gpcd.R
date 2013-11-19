# Use the GPCD regression from WaterDemandModel.R to create a table projecting
# water demand for all climate models listed in model_selection.R.
# Predicted demand amounts are stored in projected.demand

source('WaterDemandModel.R')
source('load_tasmax.R')
source('load_pr.R')
source('model_selection.R')

Projected.Yearly.GPCD <- function(ldf, model)
{
  tasmax.data <- ldf$tasmax[ , c('Year','Month','bounded.year')]
  tasmax.data$avg.max.temp.c <- ldf$tasmax[ , model]
  
  pr.data <- ldf$pr[ , c('Year','Month','bounded.year')]
  pr.data$rain.rate.mm <- ldf$pr[ , model]
  
  data <- merge(x=tasmax.data, y=pr.data, by=c('Year','Month','bounded.year'))
  
  projected.gpcd <- predict(gpcd.lm, data)
  projected.gpcm <- projected.gpcd * rep(days.per.month, nrow(ldf$tasmax)/12)
  projected.gpcy <- aggregate(projected.gpcm ~ data$Year, FUN=sum)
  names(projected.gpcy) <- c('Year',model)
  projected.gpcy[model] <- projected.gpcy[model] / sum(days.per.month)
  return(projected.gpcy)
}

projection = list()
projection$tasmax <- monthly.tmax[[grid.type]][ , c('Year', 'Month', model.selection)]
projection$pr <- monthly.pr[[grid.type]][ , c('Year', 'Month', model.selection)]

projection$tasmax$bounded.year <- Bounded.Year(projection$tasmax$Year)
projection$pr$bounded.year <- Bounded.Year(projection$pr$Year)

projected.demand <- Projected.Yearly.GPCD(projection, names(projection$tasmax)[3])

for(nn in names(projection$tasmax)[-c(1:3, which(names(projection$tasmax)=='bounded.year'))])
{
  temp <- Projected.Yearly.GPCD(projection, nn)
  projected.demand <- merge(x=projected.demand, y=temp, by='Year')
  rm(temp)
}

# Save with:
#  write.csv(projected.demand, file='/path/to/file', row.names=FALSE, quote=FALSE)

#------------------------------------------------------------------------------
# For plotting the fit to the historical demand
fit <- list()
fit$tasmax <- demand[ , c('Year', 'Month', 'avg.max.temp.c')]
names(fit$tasmax$avg.max.temp.c) <- 'Historic'
fit$pr <- demand[ , c('Year', 'Month', 'rain.rate.mm')]
names(fit$tasmax)[3] <- 'Fitted'
names(fit$pr)[3] <- 'Fitted'

fit$tasmax$bounded.year <- Bounded.Year(fit$tasmax$Year)
fit$pr$bounded.year <- Bounded.Year(fit$pr$Year)

historical.demand <- Projected.Yearly.GPCD(fit,names(fit$tasmax[3]))

temp <- demand$GPCD * rep(days.per.month, nrow(demand)/12)
temp <- aggregate(temp ~ demand$Year, FUN=sum)
names(temp) <- c('Year','Observed')
temp$Observed <- temp$Observed / sum(days.per.month)
historical.demand <- merge(x=historical.demand, y=temp, by='Year')
rm(temp)

long.historical <- melt(data=historical.demand, id.vars='Year', value.name='demand', variable.name='Data')

aspectratio <- 2/(1+sqrt(5))
width <- (8.5-2*1) # Paper width with 1" margins
height <- 0.75*width
my.theme <- theme_bw(base_size=16) +
  theme(aspect.ratio = aspectratio )

p <- qplot(x=Year, y=demand, data=long.historical, color=Data) +
  scale_y_continuous(name='Water Demand (GPCD)') +
  my.theme
print(p)
# End fit for historical record
#------------------------------------------------------------------------------
