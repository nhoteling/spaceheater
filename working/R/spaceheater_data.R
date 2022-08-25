# Read data
library(dplyr)
library(stringr)

df.cities <- read.csv("../data/us_cities.csv")
df.breweries <- read.csv("../data/us_breweries.csv") %>% 
  select(byname, lon, lat, type) 
d <- lapply(seq_len(nrow(df.breweries)), function(i) {pt <- st_point(c(df.breweries$lon[i],df.breweries$lat[i]))})
df.breweries$geometry <- st_sfc(d) %>% st_set_crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')


# Select some places for testing the functions
places_of_interest <- data.frame( city = c("Portland", "San Francisco", "San Diego", "Denver", "Chicago", "New York"),
                                  state = c("Oregon", "California", "California", "Colorado", "Illinois", "New York"))


# Get lon/lat values for each city/state and convert to spatial point via sf
d <- lapply(seq_len(nrow(places_of_interest)), function(i) {
  lon <- df.cities$longitude[ df.cities$state == places_of_interest$state[i] &  df.cities$city == places_of_interest$city[i] ]
  lat <- df.cities$latitude[ df.cities$state == places_of_interest$state[i] &  df.cities$city == places_of_interest$city[i] ]
  pt <- st_point(c(lon,lat))
})
df.places <- cbind(places_of_interest, st_sfc(d) %>% st_set_crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs') )


# Compute distance from each city-of-interest
d <- lapply(seq_len(nrow(df.places)), function(i) {
  v <- as.integer(st_distance(df.places$geometry[i], df.breweries$geometry))/1000  # convert to km
})
df.tmp <- data.frame(do.call(cbind,d))
names(df.tmp) <- paste("dist_",str_replace(df.places$city, pattern=" ", ""), sep="")
df.breweries <- cbind(df.breweries, df.tmp)



# Datasets to look at later
df.pdx <- df.breweries %>% filter(dist_Portland < 150)
df.chi <- df.breweries %>% filter(dist_Chicago < 200)