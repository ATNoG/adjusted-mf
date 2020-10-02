library("rnoaa")
library("RcppCNPy")

#------------------------------------------------------------------------
# NOTE
# In order to collect the data, a token is required 
# (it should be requested from https://www.ncdc.noaa.gov/cdo-web/token)
# the token should be stored in variable "token" (token='...')
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# IBERIAN_TEMP DATASET
#------------------------------------------------------------------------
# get lists of stations in the zone
stations_details = ncdc_stations(extent=c(37,-9,42,-5),token=token)
stationsids=stations_details$data$id
idxs=c(4,5,14,15,16,17,18,19,25)
stationsids=stationsids[idxs] 

# get data
data = matrix(0,length(stationsids),365)
discard_idxs=c()
for (i in 1:length(stationsids)){
  station_data=ncdc(datasetid = 'GHCND', datatypeid='TAVG', stationid = stationsids[i], startdate = '2010-01-01',enddate = '2010-12-31',limit=400,token=token)
  if (dim(station_data$data)[1]==365){
    data[i,]=station_data$data$value
  }else{discard_idxs[length(discard_idxs)+1]=i}
}

# store data matrix
npySave('iberian_temp.npy',data/10)


#------------------------------------------------------------------------
# WORLD_TEMP DATASET
#------------------------------------------------------------------------
# set list of stations to consider
stationsids=c('GHCND:USW00026616','GHCND:USW00014734','GHCND:USW00012919', 'GHCND:ARM00087532', 'GHCND:BR00D6-0010',  'GHCND:PEM00084628','GHCND:SFM00068538','GHCND:MOM00060150','GHCND:PO000008562', 'GHCND:SP000008280','GHCND:FRE00106209' , 'GHCND:UKM00003917')

#get data
data = matrix(0,length(stationsids),365)
discard_idxs=c()
for (i in 1:length(stationsids)){
  station_data=ncdc(datasetid = 'GHCND', datatypeid='TAVG', stationid = stationsids[i], startdate = '2019-01-01',enddate = '2019-12-31',limit=400,token=token)
  if (dim(station_data$data)[1]==365){
    data[i,]=station_data$data$value
  }else{discard_idxs[length(discard_idxs)+1]=i}
}

# store data matrix
npySave('world_temp.npy',data)

