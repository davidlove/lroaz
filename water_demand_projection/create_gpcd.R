# Use the GPCD regression from WaterDemandModel.R to create a table projecting
# water demand for all climate models listed in model_selection.R.
# Predicted demand amounts are stored in projected.demand

source('WaterDemandModel.R')
source('load_tasmax.R')
source('load_pr.R')
source('model_selection.R')

projection = list()
projection$tasmax <- monthly.tmax[[grid.type]][ , c('Year', 'Month', model.selection)]
projection$pr <- monthly.pr[[grid.type]][ , c('Year', 'Month', model.selection)]

# For plotting the fit to the historical demand
#projection$tasmax <- demand[ , c('Year', 'Month', 'avg.max.temp.c')]
#names(projection$tasmax$avg.max.temp.c) <- 'Historic'
#projection$pr <- demand[ , c('Year', 'Month', 'rain.rate.mm')]
#names(projection$tasmax)[3] <- 'Historic'
#names(projection$pr)[3] <- 'Historic'
# End fit for historical record

projection$tasmax$bounded.year <- Bounded.Year(projection$tasmax$Year)
projection$pr$bounded.year <- Bounded.Year(projection$pr$Year)

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

projected.demand <- Projected.Yearly.GPCD(projection, names(projection$tasmax)[3])

for(nn in names(projection$tasmax)[-c(1:3, which(names(projection$tasmax)=='bounded.year'))])
{
  temp <- Projected.Yearly.GPCD(projection, nn)
  projected.demand <- merge(x=projected.demand, y=temp, by='Year')
  rm(temp)
}

# Save with:
#  write.csv(projected.demand, file='/path/to/file', row.names=FALSE, quote=FALSE)
