# Create data.frame 'total' comparing average daily high temperature by month
# for the projections and the historic data

source('load_tasmax.R')
source('model_selection.R')
source('WaterDemandModel.R')

historic <- demand[ , c('Year','Month','avg.max.temp.c')]
projection <- monthly.tmax[[grid.type]][ , c('Year', 'Month', model.selection)]

Name.Months <- function(months)
{
  month.names <- list('Jan'=1, 'Feb'=2, 'Mar'=3, 'Apr'=4, 'May'=5, 'Jun'=6,
                      'Jul'=7, 'Aug'=8, 'Sep'=9, 'Oct'=10, 'Nov'=11, 'Dec'=12)
  months <- as.factor(months)
  levels(months) <- month.names
  return(months)
}

historic$Month <- Name.Months(historic$Month)
projection$Month <- Name.Months(projection$Month)

hist.temp.by.month <- aggregate(avg.max.temp.c ~ Month, FUN=mean, data=historic)
proj.temp.by.month <- aggregate(. ~ Month, FUN=mean, 
                                data=projection[ projection$Year >= 2040, names(projection) != 'Year'])

total <- merge(x=hist.temp.by.month, y=proj.temp.by.month, by='Month')
total <- total[order(total$Month), ]
