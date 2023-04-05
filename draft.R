nc <- nc_open('data/IPSL-CM6A-LR/tas/tas_IPSL-CM6A-LR_historical_r1i1p1f1_19500101-20141230_merged_regridded_v20190614.nc')
var <- 'tas'
# Temporal ranges
year_present <- 1970:1990
year_future <- 1999:2019
yyyy <- substr(as.character(nc.get.time.series(nc)), 1, 4)
iyyyy <- which(yyyy %in% year_present)
tmp <- ncvar_get(nc, var, start = c(1, 1, min(iyyyy)), count = c(-1, -1, length(iyyyy)))
