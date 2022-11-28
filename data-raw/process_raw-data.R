# Kucera et al. 2005 MARGO Atlantic data
library("readxl")
library("tidyverse")


if (!file.exists("data-raw/MARGO_Atlantic_modern_PF.xls")) {
  download.file(
    url = "http://store.pangaea.de/Projects/MARGO/regionaldata/MARGO_Atlantic_modern_PF.xls",
    destfile = "data-raw/MARGO_Atlantic_modern_PF.xls"
  )
}

N_Atlantic <- read_excel(
  path = "data-raw/MARGO_Atlantic_modern_PF.xls",
  sheet = "North Atlantic", skip = 1
) %>%
  mutate(
    across(
      c(`Longitude (decimal, from -180 to +180)`, `Water depth (m)`),
      as.numeric
    )
  )

S_Atlantic <- read_excel(
  path = "data-raw/MARGO_Atlantic_modern_PF.xls",
  sheet = "South Atlantic", skip = 1
)


Atlantic <- bind_rows(
  N_Atlantic,
  S_Atlantic %>% filter(!Core %in% N_Atlantic$Core) # remove overlap
) %>%
  select(Core:`Total Planktics`) %>%
  select(
    -`Coring device`,
    -`Water depth (m)`,
    -(Ocean:`date of addition`),
    -(`Other ID`:`Total Planktics`)
  ) %>%
  rename(
    Longitude = `Longitude (decimal, from -180 to +180)`,
    Latitude = `Latitude (decimal, from -90 to +90)`
  ) %>%
  # remove duplicate taxa
  select(
    -`Globigerinoides ruber total`, -`Globigerinoides sacc total`,
    -`Neogloboquadrina pachyderma R`, -`Globorotalia truncatulinoides L`,
    -`Globorotalia truncatulinoides R`, -`Globorotalia menardii`,
    -`Globorotalia tumida`, -`Globorotalia menardii flexuosa`
  ) %>%
  distinct(Latitude, Longitude, .keep_all = TRUE) # remove dupllicate samples



# remove rare taxa
Atlantic <- Atlantic %>%
  gather(key = taxon, value = percent, -(Core:Longitude)) %>%
  filter(percent > 0) %>%
  group_by(taxon) %>%
  filter(n() > 50) %>%
  ungroup() %>%
  spread(key = taxon, value = percent, fill = 0)


# get 1° seasonal temperatures from WOA
if (!file.exists("data-raw/woa13_decav_t13mn01v2.csv.gz")) {
  download.file(
    url = "http://data.nodc.noaa.gov/woa/WOA13/DATAv2/temperature/csv/decav/1.00/woa13_decav_t13mn01v2.csv.gz",
    destfile = "data-raw/woa13_decav_t13mn01v2.csv.gz"
  )
}
if (!file.exists("data-raw/woa13_decav_t14mn01v2.csv.gz")) {
  download.file(
    url = "http://data.nodc.noaa.gov/woa/WOA13/DATAv2/temperature/csv/decav/1.00/woa13_decav_t14mn01v2.csv.gz",
    destfile = "data-raw/woa13_decav_t14mn01v2.csv.gz"
  )
}
if (!file.exists("data-raw/woa13_decav_t15mn01v2.csv.gz")) {
  download.file(
    url = "http://data.nodc.noaa.gov/woa/WOA13/DATAv2/temperature/csv/decav/1.00/woa13_decav_t15mn01v2.csv.gz",
    destfile = "data-raw/woa13_decav_t15mn01v2.csv.gz"
  )
}
if (!file.exists("data-raw/woa13_decav_t16mn01v2.csv.gz")) {
  download.file(
    url = "http://data.nodc.noaa.gov/woa/WOA13/DATAv2/temperature/csv/decav/1.00/woa13_decav_t16mn01v2.csv.gz",
    destfile = "data-raw/woa13_decav_t16mn01v2.csv.gz"
  )
}




woa13 <- read_csv("data-raw/woa13_decav_t13mn01v2.csv.gz", skip = 1) %>%
  rename(LATITUDE = `#COMMA SEPARATED LATITUDE`) %>%
  select(LATITUDE, LONGITUDE, `50`) %>%
  filter(between(LONGITUDE, -120, 50), between(LATITUDE, -60, 90))
woa14 <- read_csv("data-raw/woa13_decav_t14mn01v2.csv.gz", skip = 1) %>%
  rename(LATITUDE = `#COMMA SEPARATED LATITUDE`) %>%
  select(LATITUDE, LONGITUDE, `50`) %>%
  filter(between(LONGITUDE, -120, 50), between(LATITUDE, -60, 90))
woa15 <- read_csv("data-raw/woa13_decav_t15mn01v2.csv.gz", skip = 1) %>%
  rename(LATITUDE = `#COMMA SEPARATED LATITUDE`) %>%
  select(LATITUDE, LONGITUDE, `50`) %>%
  filter(between(LONGITUDE, -120, 50), between(LATITUDE, -60, 90))

woa16 <- read_csv("data-raw/woa13_decav_t16mn01v2.csv.gz", skip = 1) %>%
  rename(LATITUDE = `#COMMA SEPARATED LATITUDE`) %>%
  select(LATITUDE, LONGITUDE, `50`) %>%
  filter(between(LONGITUDE, -120, 50), between(LATITUDE, -60, 90))

# get thermal summer temperature
woa <- bind_rows(
  jfm = woa13,
  amj = woa14,
  jas = woa15,
  ond = woa16,
  .id = "season"
) %>%
  filter(!is.na(`50`)) %>%
  group_by(LONGITUDE, LATITUDE) %>%
  summarise(summ50 = max(`50`))

ggplot(woa, aes(x = LONGITUDE, y = LATITUDE, fill = summ50)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = 15)

woa <- na.omit(woa)

# interpolate temperatures to site

Atlantic$summ50 <- interp::interp(
  x = woa$LONGITUDE, y = woa$LATITUDE, z = woa$summ50,
  output = "points", xo = Atlantic$Longitude, yo = Atlantic$Latitude
)$z

Atlantic %>% select(Core, Latitude, Longitude, summ50, everything())

mp <- map_data("world")
ggplot(Atlantic, aes(x = Longitude, y = Latitude, colour = summ50)) +
  geom_map(
    map = mp, data = mp, aes(map_id = region),
    inherit.aes = FALSE,
    fill = "grey50"
  ) +
  geom_point()

# save data
usethis::use_data(Atlantic, overwrite = TRUE)


#### STOR ####

STOR <- read_csv(
  file =
    "depthup,  depthdo, cageup,  cagedo
69.5,    70.5,    -56,    -46
119.5,   120.0,   3465,   3565
156.6,   157.0,   5775,   5925
204.5,   205.0,   6760,   6900
244.5,   245.0,   8410,   8625
274.5,   275.0,   9005,   9085
296.5,   297.0,   7880,   8010
314.5,   315.0,  10235,  10400
334.5,   335.0,  10605,  11060
349.5,   350.0,  11255,  11555
389.5,   390.5,  13450,  13820"
)

usethis::use_data(STOR, overwrite = TRUE)


#### arctic pollen ####
arctic <- read_excel("data-raw/pollen-climate data.xls", sheet = "POLLEN+CLIMATE+AVHRR", skip = 14)

# remove dupicate sites (pers comm)
# ID2 3893 is identical to ID2 3930
# ID2 330 is identical to ID2 3984
# ID2 213 is identical to ID2 4124
arctic <- arctic %>% filter(!ID2 %in% c(3930, 3984, 4124))

arctic.pollen <- arctic %>%
  select(starts_with("F-")) %>%
  as.data.frame()
arctic.env <- arctic %>%
  select(-matches("^[FC]-")) %>%
  as.data.frame()

usethis::use_data(arctic.pollen, overwrite = TRUE)
usethis::use_data(arctic.env, overwrite = TRUE)
