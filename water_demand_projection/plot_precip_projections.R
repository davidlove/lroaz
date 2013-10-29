# Plot precipitation projections for each model in model_selection.R, compared to the historical data

require(reshape2)
require(ggplot2)

source('load_pr.R')
source('model_selection.R')

yearly.pr <- list()
for(nn in names(monthly.pr))
{
  yearly.pr[[nn]] <- aggregate(. ~ Year, FUN=mean, data=monthly.pr[[nn]])
  yearly.pr[[nn]]$Month <- NULL
}

datepaste <- function(array) 
{
  return( paste( paste(array, sep='', collapse='-'), '-1', sep='' ) )
}

for(nn in names(monthly.pr))
{
  monthly.pr[[nn]]$Date <- apply(monthly.pr[[nn]][ , c('Year','Month')], 1, datepaste)
}

long.monthly <- melt(data=monthly.pr[[grid.type]][ , !(names(monthly.pr[[grid.type]]) %in% c('Year','Month'))],
                     id.vars=c('Date'),
                     variable.name='Model',
                     value.name='Precipitation')
long.yearly <- melt(data=yearly.pr[[grid.type]],
                     id.vars=c('Year'),
                     variable.name='Model',
                     value.name='Precipitation')
long.subset <- subset(long.yearly, Model %in% model.selection)

source('WaterDemandModel.R')

historic.pr <- aggregate(rain.rate.mm ~ Year, FUN=mean, data=demand)
names(historic.pr) <- c('Year', 'Historic')
historic.pr <- melt(data=historic.pr,
                      id.vars=c('Year'),
                      variable.name='Model',
                      value.name='Precipitation')
long.subset <- rbind(long.subset, historic.pr)

print(qplot(x=Year, y=Precipitation, data=long.subset, color=Model, geom='line')
      + theme_bw())
