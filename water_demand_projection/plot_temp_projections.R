# Plot temperature projections for each model in model_selection.R, compared to the historical data

require(reshape2)
require(ggplot2)

source('load_tasmax.R')
source('model_selection.R')

yearly.tmax <- list()
for(nn in names(monthly.tmax))
{
  yearly.tmax[[nn]] <- aggregate(. ~ Year, FUN=max, data=monthly.tmax[[nn]])
  yearly.tmax[[nn]]$Month <- NULL
}

datepaste <- function(array) 
{
  return( paste( paste(array, sep='', collapse='-'), '-1', sep='' ) )
}

for(nn in names(monthly.tmax))
{
  monthly.tmax[[nn]]$Date <- apply(monthly.tmax[[nn]][ , c('Year','Month')], 1, datepaste)
}

long.monthly <- melt(data=monthly.tmax[[grid.type]][ , !(names(monthly.tmax[[grid.type]]) %in% c('Year','Month'))],
                     id.vars=c('Date'),
                     variable.name='Model',
                     value.name='Temperature')
long.yearly <- melt(data=yearly.tmax[[grid.type]],
                     id.vars=c('Year'),
                     variable.name='Model',
                     value.name='Temperature')
long.subset <- subset(long.yearly, Model %in% model.selection)

source('WaterDemandModel.R')

historic.tmax <- aggregate(avg.max.temp.c ~ Year, FUN=max, data=demand)
names(historic.tmax) <- c('Year', 'Historic')
historic.tmax <- melt(data=historic.tmax,
                      id.vars=c('Year'),
                      variable.name='Model',
                      value.name='Temperature')
long.subset <- rbind(long.subset, historic.tmax)

print(qplot(x=Year, y=Temperature, data=long.subset, color=Model, geom='line')
      + theme_bw())
