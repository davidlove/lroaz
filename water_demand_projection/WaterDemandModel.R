# Create regression of per-capita water demand (in GPCD) from the historical Tucson Water data.
# Note: we have to regress total population on the service counts for Tucson
# Water to estimate Tucson's population per month.

require(gdata)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

pop <- read.xls('Book2.xlsx', sheet='Population', 
                method = 'tab', stringsAsFactors=FALSE, strip.white=TRUE)
pop$Population <- as.numeric(gsub(',', '', trim(pop$Population)))
services <- read.xls('Book2.xlsx', sheet='Services')
demand <- read.xls('Book2.xlsx', sheet='Sheet1', strip.white=TRUE)

base.month <- 7
s.months = services$Month == base.month
pp <- pp <- services[s.months, ]
pp$Population <- pop$Population

# Only single family homes reaches statistical significance
#pop.lm <- lm(Population ~ Single.Fam + Duplex.Triplex + Multifamiily + Commercial + Industrial.TUSD, data=pp)
pop.lm <- lm(Population ~ Single.Fam, data=pp)

print(summary(pop.lm))

demand$Population <- predict(pop.lm, services)
days.per.month <- c(31,28,31,30,31,30,31,31,30,31,30,31)

Bounded.Year <- function(year)
{
  min.year <- 2004
  max.year <- 2011
  return(pmin(pmax(year, min.year), max.year))
}

demand <- within(demand, 
                 {
                   GPCD <- Demand / Population * 1e6 / rep(days.per.month, nrow(demand)/12) # Gallons Per Capita per Day
                   rain.rate.mm <- Rain. * 25.4 / rep(days.per.month, nrow(demand)/12) # Convert to mm/day
                   avg.max.temp.c <- (AveMaxTemp - 32)*5/9 # Convert to celsius
                   bounded.year <- Bounded.Year(demand$Year)
                 })

# Linear regression on demand, including a linear fit in the year
gpcd.lm <- lm(GPCD ~ avg.max.temp.c + Year + rain.rate.mm, data=demand)
print(summary(gpcd.lm))

# Linear regression on demand, with upper and lower bounds on the year
gpcd.lm <- lm(GPCD ~ avg.max.temp.c + bounded.year + rain.rate.mm, data=demand)
print(summary(gpcd.lm))
