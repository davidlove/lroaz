# Compare the values of tasmax for the daily and monthly projections.
# This script shows that the tasmax variable in the monthly projections corresponds to the average daily high
# temperature from the daily projections.  For each year, we plot the highest result from each month.
# 
# Thus we end up with:
#  monthtmax: The highest monthly temperature from the monthly projection
#  avgdaytmax: The highest average daily high over the 12 months in each year
#  maxdaytmax: The highest daily high in each month, maximizing over the year.
#              Equivalent to the highest temperature projected for a given year.

require(reshape2)
require(ggplot2)

source('load_tasmax.R')

yearly.tmax <- list()
for(nn in names(monthly.tmax))
{
  yearly.tmax[[nn]] <- aggregate(. ~ Year, FUN=max, data=monthly.tmax[[nn]])
  yearly.tmax[[nn]]$Month <- NULL
}

yearly.maxtmax <- list()
yearly.avgtmax <- list()
for(nn in names(daily.tmax))
{
  yearly.maxtmax[[nn]] <- aggregate(. ~ Year, FUN=max, data=daily.tmax[[nn]])
  yearly.maxtmax[[nn]]$Month <- NULL
  
  yearly.avgtmax[[nn]] <- aggregate(. ~ Year + Month, FUN=mean, data=daily.tmax[[nn]])
  yearly.avgtmax[[nn]] <- aggregate(. ~ Year, FUN=max, data=yearly.avgtmax[[nn]])
  yearly.avgtmax[[nn]]$Month <- NULL
}

data.type <- 'bc5'

common.models <- intersect(names(yearly.tmax[[data.type]]), names(yearly.maxtmax[[data.type]]))
common.models <- common.models[common.models != 'Year']

cm <- common.models[1]

compare <- yearly.tmax[[data.type]][ , names(yearly.tmax[[data.type]]) %in% c('Year',cm)]
names(compare) <- c('Year', 'monthtmax')
compare$maxdaymax <- yearly.maxtmax[[data.type]][ , cm]
compare$avgdaymax <- yearly.avgtmax[[data.type]][ , cm]

long.compare <- melt(data=compare, id.vars=c('Year'), variable.name='Projection', value.name='Temperature')
print(qplot(x=Year, y=Temperature*9/5+32, data=long.compare, color=Projection))


