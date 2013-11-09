# Create data.frame 'total' comparing average daily high temperature by month
# for the projections and the historic data

source('load_tasmax.R')
source('model_selection.R')
source('WaterDemandModel.R')

# Only use RCP8.5 scenarios
model.selection <- model.selection[grep('rcp85$' ,model.selection)]

# Gather historic and projection data
historic <- demand[ , c('Year','Month','avg.max.temp.c')]
projection <- monthly.tmax[[grid.type]][ , c('Year', 'Month', model.selection)]

# Only look at years >= first.year in projection
first.year <- 2040
projection <- projection[projection$Year >= first.year, ]

# Simplify names of projections
names(historic)[3] <- 'Historic'
names(projection) <- gsub('\\.1\\.rcp85$', '', names(projection))

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

hist.temp.by.month <- aggregate(Historic ~ Month, FUN=mean, data=historic)
proj.temp.by.month <- aggregate(. ~ Month, FUN=mean, 
                                data=projection[, names(projection) != 'Year'])

agg.table <- merge(x=hist.temp.by.month, y=proj.temp.by.month, by='Month')
agg.table <- agg.table[order(agg.table$Month), ]
agg.table[, names(agg.table)!='Month'] <- round(agg.table[, names(agg.table)!='Month'], 1)
