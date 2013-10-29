# Represents an early attempt to estimate water demand using data from all 12 months simultaneously.
#
# total.demand.lasso projects demand based on all 12 monthly measurements of temperature and rain, and population.
#                    It badly overestimates usage for the lower populations necessary for RESIN projections.
# percap.demand.lasso projects per-capita demand on all 12 monthly measurements of temperature and rain.
#                     This results in a constant being the best fit for the data.

require(reshape2)
require(lars)

source('WaterDemandModel.R')

demand$Month <- as.factor(demand$Month)
month.names <- list('Jan'=1, 'Feb'=2, 'Mar'=3, 'Apr'=4, 'May'=5, 'Jun'=6,
                    'Jul'=7, 'Aug'=8, 'Sep'=9, 'Oct'=10, 'Nov'=11, 'Dec'=12)
levels(demand$Month) <- month.names

yearly.demand <- yearly.demand <- aggregate(Demand ~ Year, data=demand, FUN=sum)
yearly.demand$Population <- pop$Population

Cast.Month <- function(var, name)
{
  short.temp <- dcast( demand[ , c('Year','Month',var)], Year ~ Month, value.var=var )
  nst <- names(short.temp)
  names(short.temp) <- c(nst[1], paste(name, nst[-1], sep='.'))
  rm(nst)
  return(short.temp)
}

short.temp <- Cast.Month('avg.max.temp.c', 'Temperature')
yearly.demand <- merge(yearly.demand, short.temp, by='Year')

short.temp <- Cast.Month('rain.rate.mm', 'Rain')
yearly.demand <- merge(yearly.demand, short.temp, by='Year')

rm(short.temp)

total.demand.lasso <- lars(x=as.matrix(yearly.demand[ , -(1:2)]), y=yearly.demand$Demand)
cv.total.demand <- cv.lars(x=as.matrix(yearly.demand[ , -(1:2)]), y=yearly.demand$Demand, K=7)

percap.demand.lasso <- lars(x=as.matrix(yearly.demand[ , -(1:3)]), y=yearly.demand$Demand / yearly.demand$Population * 1e6)
cv.percap.demand <- cv.lars(x=as.matrix(yearly.demand[ , -(1:3)]), y=yearly.demand$Demand / yearly.demand$Population * 1e6, K=7)

predict(total.demand.lasso, type='coefficients', mode='fraction', 
        s=cv.total.demand$index[which.min(cv.total.demand$cv + cv.total.demand$cv.error)])$coefficients
