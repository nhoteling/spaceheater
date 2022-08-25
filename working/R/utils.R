# utils.R


#url <- "https://dwr.virginia.gov/fishing/trout-stocking-schedule/?start_date=January+1%2C+2015&end_date=January+3%2C+2022"


build_trout_url <- function(tstart, tend) {
  # use current date for upper limit
  url0 <- "https://dwr.virginia.gov/fishing/trout-stocking-schedule"
  url1 <- paste0(url0, 
                 "/?",
                 "start_date=",format.Date(tstart,"%B"),"+",format.Date(tstart,"%d"),"%2C+",format.Date(tstart,"%Y"),
                 "&end_date=",format.Date(tend,"%B"),"+",format.Date(tend,"%d"),"%2C+",format.Date(tend,"%Y"))
  return(url1)
}


cleanup_trout_data <- function(df) {
  df %>%
    rename(date_string = Date, 
           waterbody_string = Waterbody,
           species = `Species Stocked`,
           county_string = County) %>%
    mutate(date = strptime(date_string, format="%B %d, %Y")) %>%
    filter(waterbody_string != "No Stockings Today") %>%
    mutate(waterbody1 = waterbody_string %>% str_replace_all("\\((.*?)\\)","") %>% str_replace_all("\\[(.*?)\\]","") %>% str_trim(),
           waterbody2 = stringr::str_extract_all(waterbody_string, "\\((.*?)\\)") %>% str_replace_all("[\\(\\)]",""),
           specialtype = stringr::str_extract_all(waterbody_string, "\\[(.*?)\\]") %>% str_replace_all("[\\[\\]]",""),
           waterbody2 = ifelse(waterbody2 == "character0",NA,waterbody2),
           specialtype = ifelse(specialtype == "character(0)",NA,specialtype),
           county = str_replace(county_string, "County","") %>% str_replace("City of","") %>% str_trim()) %>%
    dplyr::select(-date_string, -waterbody_string)
}


make_line <- function(lonlat) { lonlat %>% as.matrix() %>% sf::st_linestring() }

make_point <- function(lon, lat, crs=4326) {
  sf::st_point(c(lon,lat)) %>% sf::st_sfc() %>% sf::st_set_crs(crs)
}


crd_to_points <- function(lon,lat, CRS=4326) {
  d <- lapply(1:length(lon), function(i) {sf::st_point(c(lon[i],lat[i]))})
  return(d %>% sf::st_sfc(crs=CRS))
}


# NOT USED
find_close_match <- function(s, v) {
  stringdist::stringdist(s, v, method="jw")
}


# convert osmdata bbox to pgon
bbox_to_pgon <- function(v, CRS=4326) {
  x <- v %>% as.vector()
  if (!anyNA(x)) {
    names(x) <- c("xmin", "ymin","xmax","ymax")
    pgon <- x %>% sf::st_bbox() %>% sf::st_as_sfc() %>% sf::st_set_crs(CRS)
  } else {
    pgon <- NA
  }
  return(pgon) 
}
