rm(list=ls(all=T))
require(raster)
require(fields)
require(RSQLite)
require(spatstat)
require(maps)
require(ncdf4)

# functions ####
readmap.nc <- function(filename,var="",lo="lon",la="lat",starttime=1){
  nc <- nc_open(filename)
  if(var=="") var <- names(nc$var)[1]
  lon <- ncvar_get(nc,lo)
  if(min(lon)>=0){
    cat("WARNING! Longitude does not contain negative values, shifting >180 by 360\n")
    lon[lon>180] <- lon[lon>180]-360
  }
  lat <- ncvar_get(nc,la)
  if(lat[1] > lat[length(lat)]){
    cat("WARNING, inverting latitudes\n")
    
  }
  if(starttime==1)
    buf <- ncvar_get(nc,var)
  else
    buf <-ncvar_get(nc,var,start=c(1,1,starttime))
  nc_close(nc)
  if(length(dim(buf))==2)
    buf <- buf[order(lon),order(lat)]
  else if(length(dim(buf))==3)
    buf <- buf[order(lon),order(lat),]
  else if(length(dim(buf))>3)
    cat("WARNING, cannot adjust lon/lat setting for 4-dim array\n")
  buf
}

find.dominant <- function(data1,data2,data3,data4,data5,data6,data7,
                          data8,data9,data10,data11,data12,data13,
                          mu1,mu2,mu3,mu4,mu5,mu6,mu7,
                          mu8,mu9,mu10,mu11,mu12,mu13,
                          ar,nsoil_pix1,nsoil_pix2,index)
{
  nsoil_pix <- nsoil_pix1*nsoil_pix2
  #cat("noilpix",nsoil_pix,nsoil_pix1,nsoil_pix2,"\n")
  dominant_mu <- NA # initialize to NA in case there is no valid soil texture class
  weight <- array(0,dim=dim(data1))
  # if landuse is zero, use dominant soil type of all land
  if(sum(ar,na.rm=T)==0)
  {
    # find most common element, use first of most common ones
    if(all(!is.finite(data1)) & all(!is.finite(data2)) & all(!is.finite(data3)) & all(!is.finite(data4)) & all(!is.finite(data5)) &
       all(!is.finite(data6)) & all(!is.finite(data7)) & all(!is.finite(data8)) & all(!is.finite(data9)) & all(!is.finite(data10)) & 
       all(!is.finite(data11)) & all(!is.finite(data12)) & all(!is.finite(data13))){
      dominant <- NA
      dominant.mu <- NA
      weight[] <- 1
    } else {
      area <- integer(13)
      for(so in 1:13){ # loop over soil texture classes
        data <- get(paste0("data",so))
        #area[so] <- length(which(is.finite(data)))
        area[so] <- sum(data,na.rm=T) 
      }
      #dominant <- as.integer(names(sort(table(data),decreasing = T)[1]))
      dominant <- which.max(area) # always takes first if several are of same max value
      data <- get(paste0("data",dominant))
      weight[is.finite(data)] <- 1
      #cat("dominant",dominant,"\n")
      #cat(weight,"\n")
    }
  } else {
    #soils <- unique(as.vector(data))
    #soils <- soils[is.finite(soils)]
    # sieve out all texture classes that have no pixels
    soils <- 1:13
    remove <- NULL
    for(i in soils){
      data <- get(paste0("data",i))
      if(all(!is.finite(data)))
        remove <- c(remove,i)
    }
    soils <- soils[-remove]
    if(length(soils)<1){
      dominant <- dominant.mu <- NA
      weight[] <- 1
    } else {
      # suffle randomly in case the first element is selected
      if(length(soils)>1) soils <- sample(soils)
      areas <- numeric(length(soils))
      pix_with_lu <- which(ar>0,arr.ind=T)
      for(i in 1:length(soils))
      {
        weight[] <- 0
        areas[i] <- 0
        # sequence of soils has been suffled, so we need to take right data array
        data <- get(paste0("data",soils[i]))
        #cat("data",dim(data),"nsoil_pix1/2",nsoil_pix1,nsoil_pix2,"\n")
        for(p in 1:dim(pix_with_lu)[1]){
          #cat("p",p,"\n")
          dat <- data[c(1:nsoil_pix1)+(pix_with_lu[p,1]-1)*nsoil_pix1,c(1:nsoil_pix2)+(pix_with_lu[p,2]-1)*nsoil_pix2]
          #ndat <- length(which(dat==soils[i]))
          ndat <- sum(dat,na.rm=T)
          if(ndat>0){
            # assumes land use is distributed equally to all finer soil data grids
            areas[i] <- areas[i] + ar[pix_with_lu[p,1],pix_with_lu[p,2]]*ndat/nsoil_pix # total area of that soil texture class
            # fill weight only for pixels that have land-use (initialized to zero above), assuming even distribution across all nsoil_pix pixels.
            weight[c(1:nsoil_pix1)+(pix_with_lu[p,1]-1)*nsoil_pix1,c(1:nsoil_pix2)+(pix_with_lu[p,2]-1)*nsoil_pix2] <- ar[pix_with_lu[p,1],pix_with_lu[p,2]]/nsoil_pix
          }
        }
        assign(paste0("weight",soils[i]),weight)
      }
      #take (random) first soil type if several are dominant
      dominant <- soils[which.max(areas)]
      weight <- get(paste0("weight",dominant))
      data <- get(paste0("data",dominant))
      weight[!is.finite(data)] <- 0 # set all weights in pixels with non-dominant soil type to zero again
      if(all(!is.finite(weight*data)) | all(weight==0)){
        weight[] <- 1
      }   
    }
  }
  # weight has values only for pixels with land use and with dominant soil type
  if(!is.na(dominant))
  {
    mu_pix <- get(paste0("mu",dominant))
    dominant_mu <- as.integer(names(sort(table(mu_pix[weight>0]),decreasing = T)[1]))
    cat("dominant",dominant,"dominant_mu",dominant_mu,"\n")
  }
  list(type=dominant,weight=weight,mu=dominant_mu) 
}

find.dominant.neighbor <- function(data,index,maxdist=200)
{
  neighbor <- integer(dim(index)[1])
  for(i in 1:dim(index)[1]){
    #cat("index",range(index),"\n")
    #str(index)
    for(dist in 1:maxdist){
      #cat(i,"\n")
      xis <- c(-dist:dist)+index[i,1]
      #cat(i,":",xis,"\n")
      remove <- which(xis>dim(data)[1])
      if(length(remove>0))
         xis <- xis[-remove]
      yis <- c(-dist:dist)+index[i,2]
      #cat(yis,"\n")
      remove <- which(yis>dim(data)[2])
      if(length(remove)>0)
        yis <- yis[-remove]
      #cat("neighbor:",dim(data),dim(index),"\n")
      #cat(xis,"\n")
      #cat(yis,"\n")
      find <- which(data[xis[xis>0],yis[yis>0]]>0,arr.ind=T)
      if(dim(find)[1]>0){
        #cat(i,": find for dist",dist,"\n")
        #str(find)
        neighbor[i] <- as.integer(names(sort(table(data[xis[xis>0],yis[yis>0]]),decreasing = T)[1]))
        #cat(neighbor[i],"\n")
        break
      }
    }
  }
  neighbor
}

weighted.median2 <- function(x,w){
  if(all(!is.finite(x*w))){
    return(NA)
  } else {
    weighted.median(x,w)
  }
}

# settings ####
setwd("/p/projects/macmit/data/GGCMI/phase3/hwsd/")

read.df <- T
get.data <- F
process.hwsd <- F
do.aggregation <- T
get.landuse <- F
process.maps <- F
do.cropland.weighting <- T
do.allland <- T

# processing HWSD data
if(get.data){
  # get data: HWSD 1.21  and prepare ####
  
  system("wget http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HWSD_Data/HWSD_RASTER.zip")
  system("unzip HWSD_RASTER.zip")
  system("wget https://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HWSD_Data/HWSD.mdb")
  # convert mdb to sqlite
  # use tool from 
  # https://sourceforge.net/projects/nsbase/files/latest/download
  
  # version 2.0 or later
  # convert mdb to csv
  # mdbtools-0.8.2/bin/mdb-export HWSD.mdb HWSD_DATA > hwsd_data.csv
  # mdbtools-0.8.2/bin/mdb-export HWSD.mdb HWSD_DATA > hwsd_data.csv
}

if(process.hwsd){
  hwsd <- raster("hwsd.bil")
  proj4string(hwsd) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

  
  # used before version 2.0
  
  # m <- dbDriver("SQLite")
  # con <- dbConnect(m,dbname="HWSD.sqlite")
  # # some testing
  # #dbListTables(con)
  # #dbGetQuery(con, "pragma table_info(HWSD_DATA)")$name
  # #dbGetQuery(con, "pragma table_info(HWSD_DATA)")$type
  # #dbGetQuery(con, "select count(*) as grid_total from HWSD_DATA")
  # 
  # dbWriteTable(con, name="global_HWSD", value=data.frame(smu_id=unique(hwsd)), overwrite=T)
  # 
  # records <- dbGetQuery(con, "select T.* from HWSD_DATA as T join global_HWSD as U on T.mu_global=u.smu_id order by su_sym90")
  # #coverage <- dbGetQuery(con, "select * from D_COVERAGE")
  # # convert strings to numeric for pH values (also other variables have this issue, but we're not processing them here)
  # records$T_PH_H2O <- as.numeric(gsub(",", ".", records$T_PH_H2O))
  # records$T_BULK_DENSITY <- as.numeric(gsub(",", ".", records$T_BULK_DENSITY))
  # records$T_OC <- as.numeric(gsub(",", ".", records$T_OC))
  # records$T_ECE <- as.numeric(gsub(",", ".", records$T_ECE))
  
  # new in version 2.0 and larger 
  records <- read.csv("/p/projects/macmit/data/GGCMI/phase3/hwsd/hwsd_data.csv")
  
  
  
  # gap fill all parameters with 
  
  # select of each mapping unit only the one with the largest SHARE
  muunique <- unique(records$MU_GLOBAL)

  miss.ph <- miss.caco3 <- miss.bs <- miss.ece <- miss.cec <- miss.oc <- 
    miss.gravel <- miss.silt <- miss.sand <- miss.clay <- miss.bulkdens <- 
    miss.root <- miss.il <- miss.awc <- miss.issoil <- 0
  keep.all <- NULL
  for(i in 1:length(muunique))
  {
    if(i %% 100 ==0) cat("processing",i,"of",length(muunique),"\n")
    set1 <- which(records$MU_GLOBAL==muunique[i])
    tex_classes <- unique(records$T_USDA_TEX_CLASS[set1])
    if(length(tex_classes)<1){
      cat("no valid USDA texture class in MU_GLOBAL",muunique[i],"\n")
    } else {
      # new in version 2.0 and larger: gap filling first with profiles of same USDA TEX CLASS and same MU_GLOBAL, if not possible with same MU_GLOBAL
      # loop over all soils in same USDA TEX CLASS, if gap filling is not possible with that, use all soils in MU_GLOBAL
      for(tc in 1:length(tex_classes)){
        set <- which(records$T_USDA_TEX_CLASS==tex_classes[tc] & records$MU_GLOBAL==muunique[i])
        if(length(set)>0){
          if(length(set)>1) cat("found",length(set),"unique TEX_CLASS in",muunique[i],"\n")
          # select first of those records with largest share, drop rest
          keep <- set[which.max(records$SHARE[set])[1]]
          keep.all <- c(keep.all,keep)
          drop <- set[-which(set==keep)] # all records of same MU_GLOBAL and same USDA TEX CLASS
            drop1 <- set1[-which(set1==keep)] # all records of same MU_GLOBAL
          # test if one of parameters of interest is not included in record with max SHARE, but in others
          if(!is.finite(records$T_PH_H2O[keep]))
          {
            if(any(is.finite(records$T_PH_H2O[drop])))
            {
              find <- drop[which(is.finite(records$T_PH_H2O[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              cat(i,"found in same USDA TEX CLASS",set,"keep",keep,"drop",drop,"replace ",records$T_PH_H2O[keep],"with",records$T_PH_H2O[take])
              records$T_PH_H2O[keep] <- records$T_PH_H2O[take]
              cat(" to",records$T_PH_H2O[keep],"\n")
            } else if(any(is.finite(records$T_PH_H2O[drop1])))
            {
              find <- drop[which(is.finite(records$T_PH_H2O[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              cat(i,"found in same MU_GLOBAL",set,"keep",keep,"drop",drop1,"replace ",records$T_PH_H2O[keep],"with",records$T_PH_H2O[take])
              records$T_PH_H2O[keep] <- records$T_PH_H2O[take]
              cat(" to",records$T_PH_H2O[keep],"\n")
            } else {
              miss.ph <- miss.ph + 1
            }
          }
          if(!is.finite(records$T_CACO3[keep]))
          {
            if(any(is.finite(records$T_CACO3[drop])))
            {
              find <- drop[which(is.finite(records$T_CACO3[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_CACO3[keep] <- records$T_CACO3[take]
            } else if(any(is.finite(records$T_CACO3[drop1])))
            {
              find <- drop[which(is.finite(records$T_CACO3[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_CACO3[keep] <- records$T_CACO3[take]
            } else{
              miss.caco3 <- miss.caco3 + 1
            }
          }
          if(!is.finite(records$T_BS[keep]))
          {
            if(any(is.finite(records$T_BS[drop])))
            {
              find <- drop[which(is.finite(records$T_BS[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_BS[keep] <- records$T_BS[take]
            } else if(any(is.finite(records$T_BS[drop1])))
            {
              find <- drop[which(is.finite(records$T_BS[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_BS[keep] <- records$T_BS[take]
            } else {
              miss.bs <- miss.bs + 1
            }
          }
          if(!is.finite(records$T_ECE[keep]))
          {
            if(any(is.finite(records$T_ECE[drop])))
            {
              find <- drop[which(is.finite(records$T_ECE[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_ECE[keep] <- records$T_ECE[take]
            } else if(any(is.finite(records$T_ECE[drop1])))
            {
              find <- drop[which(is.finite(records$T_ECE[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_ECE[keep] <- records$T_ECE[take]
            } else {
              miss.ece <- miss.ece + 1
            }
          }
          if(!is.finite(records$T_CEC_SOIL[keep]))
          {
            if(any(is.finite(records$T_CEC_SOIL[drop])))
            {
              find <- drop[which(is.finite(records$T_CEC_SOIL[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_CEC_SOIL[keep] <- records$T_CEC_SOIL[take]
            } else if(any(is.finite(records$T_CEC_SOIL[drop1])))
            {
              find <- drop[which(is.finite(records$T_CEC_SOIL[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_CEC_SOIL[keep] <- records$T_CEC_SOIL[take]
            } else {
              miss.cec <- miss.cec + 1
            }
          }
          if(!is.finite(records$T_OC[keep]))
          {
            if(any(is.finite(records$T_OC[drop])))
            {
              find <- drop[which(is.finite(records$T_OC[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_OC[keep] <- records$T_OC[take]
            } else if(any(is.finite(records$T_OC[drop1])))
            {
              find <- drop[which(is.finite(records$T_OC[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_OC[keep] <- records$T_OC[take]
            } else {
              miss.oc <- miss.oc + 1
            }
          }
          if(!is.finite(records$T_GRAVEL[keep]))
          {
            if(any(is.finite(records$T_GRAVEL[drop])))
            {
              find <- drop[which(is.finite(records$T_GRAVEL[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_GRAVEL[keep] <- records$T_GRAVEL[take]
            } else if(any(is.finite(records$T_GRAVEL[drop1])))
            {
              find <- drop[which(is.finite(records$T_GRAVEL[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_GRAVEL[keep] <- records$T_GRAVEL[take]
            } else {
              miss.gravel <- miss.gravel + 1
            }
          }
          if(!is.finite(records$T_SAND[keep]))
          {
            if(any(is.finite(records$T_SAND[drop])))
            {
              find <- drop[which(is.finite(records$T_SAND[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_SAND[keep] <- records$T_SAND[take]
            } else if(any(is.finite(records$T_SAND[drop1])))
            {
              find <- drop[which(is.finite(records$T_SAND[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_SAND[keep] <- records$T_SAND[take]
            } else {
              miss.sand <- miss.sand + 1
            }
          }
          if(!is.finite(records$T_SILT[keep]))
          {
            if(any(is.finite(records$T_SILT[drop])))
            {
              find <- drop[which(is.finite(records$T_SILT[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_SILT[keep] <- records$T_SILT[take]
            } else if(any(is.finite(records$T_SILT[drop1])))
            {
              find <- drop[which(is.finite(records$T_SILT[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_SILT[keep] <- records$T_SILT[take]
            } else {
              miss.silt <- miss.silt + 1
            }
          }
          if(!is.finite(records$T_CLAY[keep]))
          {
            if(any(is.finite(records$T_CLAY[drop])))
            {
              find <- drop[which(is.finite(records$T_CLAY[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_CLAY[keep] <- records$T_CLAY[take]
            } else if(any(is.finite(records$T_CLAY[drop1])))
            {
              find <- drop[which(is.finite(records$T_CLAY[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_CLAY[keep] <- records$T_CLAY[take]
            } else {
              miss.clay <- miss.clay + 1
            }
          }
          if(!is.finite(records$T_BULK_DENSITY[keep]))
          {
            if(any(is.finite(records$T_BULK_DENSITY[drop])))
            {
              find <- drop[which(is.finite(records$T_BULK_DENSITY[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_BULK_DENSITY[keep] <- records$T_BULK_DENSITY[take]
            } else if(any(is.finite(records$T_BULK_DENSITY[drop1])))
            {
              find <- drop[which(is.finite(records$T_BULK_DENSITY[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$T_BULK_DENSITY[keep] <- records$T_BULK_DENSITY[take]
            } else {
              miss.bulkdens <- miss.bulkdens + 1
            }
          }
          if(!is.finite(records$T_CACO3[keep]))
          {
            if(any(is.finite(records$ROOTS[drop])))
            {
              find <- drop[which(is.finite(records$ROOTS[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$ROOTS[keep] <- records$ROOTS[take]
            } else if(any(is.finite(records$ROOTS[drop1])))
            {
              find <- drop[which(is.finite(records$ROOTS[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$ROOTS[keep] <- records$ROOTS[take]
            } else {
              miss.root <- miss.root + 1
            }
          }
          if(!is.finite(records$IL[keep]))
          {
            if(any(is.finite(records$IL[drop])))
            {
              find <- drop[which(is.finite(records$IL[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$IL[keep] <- records$IL[take]
            } else if(any(is.finite(records$IL[drop1])))
            {
              find <- drop[which(is.finite(records$IL[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$IL[keep] <- records$IL[take]
            } else {
              miss.il <- miss.il + 1
            }
          }
          if(!is.finite(records$AWC_CLASS[keep]))
          {
            if(any(is.finite(records$AWC_CLASS[drop])))
            {
              find <- drop[which(is.finite(records$AWC_CLASS[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$AWC_CLASS[keep] <- records$AWC_CLASS[take]
            } else if(any(is.finite(records$AWC_CLASS[drop1])))
            {
              find <- drop[which(is.finite(records$AWC_CLASS[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$AWC_CLASS[keep] <- records$AWC_CLASS[take]
            } else {
              miss.awc <- miss.awc + 1
            }
          }
          if(!is.finite(records$ISSOIL[keep]))
          {
            if(any(is.finite(records$ISSOIL[drop])))
            {
              find <- drop[which(is.finite(records$ISSOIL[drop]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$ISSOIL[keep] <- records$ISSOIL[take]
            } else if(any(is.finite(records$ISSOIL[drop1])))
            {
              find <- drop[which(is.finite(records$ISSOIL[drop1]))]
              take <- find[which.max(records$SHARE[find])[1]]
              records$ISSOIL[keep] <- records$ISSOIL[take]
            } else {
              miss.issoil <- miss.issoil + 1
            }
          }
          #  cat("removing",length(drop),"from MU_GLOBAL",muunique[i],"\n")
          #records <- records[-drop,]
        }
      }
    }
  }
  records <- records[keep.all,]
  

    # convert to spatial data frame
  if(read.df){
    load("hwsd2.0.df.Rdata")
  } else {
    hwsd.df <- as(hwsd,"SpatialGridDataFrame")
    
    ma <- match(hwsd.df@data$hwsd, records$MU_GLOBAL)
    
    save(ma,hwsd.df,file="hwsd2.0.df.Rdata")
    save(records,file="hwsd.records2.0.Rdata")
  }
  
  # reduce to the important columns and convert to arrays
  #pars <- c(2,17,20:21,24:28,30:32,34,35,37,40)
  #pars <- c(2,28)
  #for(i in pars){
  # for MU_GLOBAL of max share soil type
  for(i in 1:13)
    assign(paste0("rec",i),records[which(records$T_USDA_TEX_CLASS==i),])
  
  # keep copy
  hwsd.df2 <- hwsd.df

  
  # for plotting only
  largest <- records
  for(mu in unique(largest$MU_GLOBAL)){
    set <- which(largest$MU_GLOBAL==mu)
    if(length(set)>1){
      keep <- set[which.max(largest$SHARE[set])[1]]
      drop <- set[-which(set==keep)]
      largest <- largest[-drop,]
    }
  }
  if(F){
    ma <- match(hwsd.df@data$hwsd,largest$MU_GLOBAL)
    hwsd.df@data <- data.frame(largest[ma,2])
    mapi <- as.array(hwsd.df)
    png(paste0(names(records)[2],"_global_30sec2.0.png"),width=8*300,height=4*300,res=300,pointsize=10)
    image.plot(mapi[,21600:1],x=seq(-180,180,length.out=43200),y=seq(-90,90,length.out=21600))
    map(add=T)
    dev.off()
    save(ma,mapi,file=paste0(names(records)[2],"_map_2.0.Rdata"))
  }
  if(process.maps){
    for(i in 1:13){
  #for(i in 10:13){# loop over all soil texture classes
      rec <- get(paste0("rec",i))
      #cat("i:",i,"\n")
      #str(rec)
      hwsd.df <- hwsd.df2
      ma <- match(hwsd.df@data$hwsd,rec$MU_GLOBAL)
      #assign(paste0("ma",i),ma)
      #cat("ma\n")
      #str(ma)
      hwsd.df@data <- data.frame(rec[ma,6]) # share
      mapi <- as.array(hwsd.df)
      #cat("mapi\n")
      #str(mapi)
      hwsd.df <- hwsd.df2
      hwsd.df@data <- data.frame(rec[ma,2]) # MU_GLOBAL
      mapi_mu <- as.array(hwsd.df)
      #png(paste0(names(records)[28],"_",i,"shares_global_30sec2.0.png"),width=8*300,height=4*300,res=300,pointsize=10)
      #image.plot(mapi[,21600:1],x=seq(-180,180,length.out=43200),y=seq(-90,90,length.out=21600))
      #map(add=T)
      #dev.off()
      #png(paste0(names(records)[28],"_",i,"MU_GLOBAL_global_30sec2.0.png"),width=8*300,height=4*300,res=300,pointsize=10)
      #image.plot(mapi_mu[,21600:1],x=seq(-180,180,length.out=43200),y=seq(-90,90,length.out=21600))
      #map(add=T)
      #dev.off()
      # use 'do.call' in order to make 'save' and 'get' work in combination
      #do.call(save, list(paste0("ma",i),mapi,mapi_mu,file=paste0("texture_class_",i,"_2.0.Rdata")))
      save(mapi,mapi_mu,file=paste0("texture_class_",i,"_2.0.Rdata"))
    }
    rm(ma,mapi,mapi_mu)
  }  
  #dbDisconnect(con)
}

# aggregate to 0.5 degree ####
if(do.aggregation){
  if(get.landuse){
    # read sealandmask
    landmask <- readmap.nc("/p/projects/isimip/isimip/ISIMIP3b/InputData/geo_conditions/landseamask/landseamask_water-global.nc")
    land <- which(landmask==1,arr.ind=T)
    # read landuse (from MIRCA as we want higher spatial resolution than 0.5 degrees)
    #landuse <- readmap.nc("/p/projects/isimip/isimip/ISIMIP3b/InputData/socioeconomic/landuse/2015soc/landuse-totals_2015soc_annual_2015_2100.nc")[,,1]
    landuse <- array(0,dim=c(2160,4320))
    for(i in 1:26){
      landuse <- landuse + as.matrix(raster(paste0("/p/projects/lpjml/raw_data/MIRCA2000/HA_05min/ANNUAL_AREA_HARVESTED_RFC_CROP",i,"_HA.ASC")))
      landuse <- landuse + as.matrix(raster(paste0("/p/projects/lpjml/raw_data/MIRCA2000/HA_05min/ANNUAL_AREA_HARVESTED_IRC_CROP",i,"_HA.ASC")))
    }
    landuse <- t(landuse)[,2160:1]
    save(landmask,land,landuse,file="landinformation2.0.Rdata")
  } else {
    load("landinformation2.0.Rdata")
  }
  landuse2 <- landuse
  landuse2[] <- 0
  
  # set seed for randomizing order of soils considered per 0.5 grid cell so that 
  # in case of equal area shares there is no bias in which soil type to use
  set.seed(123)
  
  for(i in 1:13){
    load(paste0("texture_class_",i,"_2.0.Rdata"))
    assign(paste0("share",i),mapi[,21600:1])
    assign(paste0("mu",i),mapi_mu[,21600:1])
    #cat(i,"dims",dim(mapi),"mu",dim(mapi_mu),"\n")
    rm(mapi,mapi_mu)
  }
  # load("T_USDA_TEX_CLASS2.0.Rdata")
  # usda_tex_class <- mapi[,21600:1]
  # rm(mapi)
  # load("MU_GLOBAL2.0.Rdata")
  # mu <- mapi[,21600:1]
  # rm(mapi)
  data_tex_class <- data_mu <- 
    data_tex_class_nolu <- data_mu_nolu <- 
    array(NA,dim=c(720,360))
  
  # soil data come in 30 arc-sec resolution, land-use areas in 5 arc-minute resolution
  # i.e. there are 100 soil pixels (10x10) per unit of land-use area
  nsoil_pix1 <- dim(share1)[1]/dim(landuse)[1]
  nsoil_pix2 <- dim(share1)[2]/dim(landuse)[2]
  #cat("nsoils 1/2",nsoil_pix1,nsoil_pix2,dim(share1),dim(landuse),"\n")
  
  for(i in 1:dim(land)[1])
  {
    if(i%%100==0) 
      cat("processing",i,"\n")
    #cat("dims",dim(share1),"mu",dim(mu1),"lu",dim(landuse),"\n")
    # cat(share1[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share2[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share3[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share4[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share5[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share6[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share7[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share8[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share9[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share10[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share11[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share12[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      share13[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],"\n")
    # cat( mu1[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu2[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu3[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu4[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu5[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu6[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu7[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu8[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu9[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu10[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu11[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu12[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #      mu13[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],"\n")
    # cat(landuse[c(1:6)+(land[i,1]-1)*6,c(1:6)+(land[i,2]-1)*6],"\n")
    # dominant <- find.dominant(usda_tex_class[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
    #                           landuse[c(1:6)+(land[i,1]-1)*6,c(1:6)+(land[i,2]-1)*6],
    #                           mu[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],nsoil_pix1,nsoil_pix2,i)
    dominant <- find.dominant(share1[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share2[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share3[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share4[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share5[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share6[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share7[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share8[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share9[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share10[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share11[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share12[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share13[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu1[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu2[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu3[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu4[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu5[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu6[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu7[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu8[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu9[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu10[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu11[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu12[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu13[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              landuse[c(1:6)+(land[i,1]-1)*6,c(1:6)+(land[i,2]-1)*6],
                              nsoil_pix1,nsoil_pix2,i)
    dominant2 <- find.dominant(share1[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share2[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share3[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share4[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share5[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share6[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share7[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share8[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share9[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share10[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share11[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share12[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              share13[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu1[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu2[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu3[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu4[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu5[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu6[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu7[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu8[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu9[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu10[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu11[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu12[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              mu13[c(1:60)+(land[i,1]-1)*60,c(1:60)+(land[i,2]-1)*60],
                              landuse2[c(1:6)+(land[i,1]-1)*6,c(1:6)+(land[i,2]-1)*6],
                              nsoil_pix1,nsoil_pix2,i)
    data_tex_class[land[i,1],land[i,2]] <- dominant$type
    data_mu[land[i,1],land[i,2]] <- dominant$mu
    data_tex_class_nolu[land[i,1],land[i,2]] <- dominant2$type
    data_mu_nolu[land[i,1],land[i,2]] <- dominant2$mu
  }
  
  save(data_tex_class,data_mu,data_tex_class_nolu,data_mu_nolu,
       file="dominant_soil_type_data2.0.Rdata")
  
  png(paste0("data_tex_class_global_30min2.0.png"),width=8*300,height=4*300,res=300,pointsize=10)
  image.plot(data_tex_class,x=seq(-180,180,length.out=720),y=seq(-90,90,length.out=360))
  map(add=T)
  dev.off()
  png(paste0("MU_global_30min2.0.png"),width=8*300,height=4*300,res=300,pointsize=10)
  image.plot(data_mu,x=seq(-180,180,length.out=720),y=seq(-90,90,length.out=360))
  map(add=T)
  dev.off()
  png(paste0("data_tex_class_nolu_global_30min2.0.png"),width=8*300,height=4*300,res=300,pointsize=10)
  image.plot(data_tex_class_nolu,x=seq(-180,180,length.out=720),y=seq(-90,90,length.out=360))
  map(add=T)
  dev.off()
  png(paste0("MU_GLOBAL_nolu_global_30min2.0.png"),width=8*300,height=4*300,res=300,pointsize=10)
  image.plot(data_mu_nolu,x=seq(-180,180,length.out=720),y=seq(-90,90,length.out=360))
  map(add=T)
  dev.off()
  
  

  load("hwsd.records2.0.Rdata")
  
  # version 2.2 and up: gap fill missing values for specific properties, Sam Rabin analyzed what records are the most similar
  # bulk-density
  records$T_BULK_DENSITY[which(records$MU_GLOBAL==26333 & records$T_USDA_TEX_CLASS==9)] <- records$T_BULK_DENSITY[which(records$MU_GLOBAL==27915 & records$T_USDA_TEX_CLASS==9)]
  cat("filling bulk density of 26333",records$T_BULK_DENSITY[which(records$MU_GLOBAL==26333 & records$T_USDA_TEX_CLASS==9)],"with 27915",records$T_BULK_DENSITY[which(records$MU_GLOBAL==27915 & records$T_USDA_TEX_CLASS==9)],"\n")
  records$T_BULK_DENSITY[which(records$MU_GLOBAL==26308 & records$T_USDA_TEX_CLASS==11)] <- records$T_BULK_DENSITY[which(records$MU_GLOBAL==26314 & records$T_USDA_TEX_CLASS==11)]
  cat("filling bulk density of 26308",records$T_BULK_DENSITY[which(records$MU_GLOBAL==26308 & records$T_USDA_TEX_CLASS==11)],"with 27915",records$T_BULK_DENSITY[which(records$MU_GLOBAL==26314 & records$T_USDA_TEX_CLASS==11)],"\n")
  # AWC
  records$AWC_CLASS[which(records$MU_GLOBAL==11385 & records$T_USDA_TEX_CLASS==9)] <- records$AWC_CLASS[which(records$MU_GLOBAL==11755 & records$T_USDA_TEX_CLASS==9)]
  cat("filling awc of 11385",records$AWC_CLASS[which(records$MU_GLOBAL==11385 & records$T_USDA_TEX_CLASS==9)],"with 11755",records$AWC_CLASS[which(records$MU_GLOBAL==11755 & records$T_USDA_TEX_CLASS==9)],"\n")
  records$AWC_CLASS[which(records$MU_GLOBAL==11271 & records$T_USDA_TEX_CLASS==9)] <- records$AWC_CLASS[which(records$MU_GLOBAL==11275 & records$T_USDA_TEX_CLASS==9)]
  cat("filling awc of 11271",records$AWC_CLASS[which(records$MU_GLOBAL==11271 & records$T_USDA_TEX_CLASS==9)],"with 11275",records$AWC_CLASS[which(records$MU_GLOBAL==11275 & records$T_USDA_TEX_CLASS==9)],"\n")
  records$AWC_CLASS[which(records$MU_GLOBAL==11268 & records$T_USDA_TEX_CLASS==9)] <- records$AWC_CLASS[which(records$MU_GLOBAL==11254 & records$T_USDA_TEX_CLASS==9)]
  cat("filling awc of 11268",records$AWC_CLASS[which(records$MU_GLOBAL==11268 & records$T_USDA_TEX_CLASS==9)],"with 11254",records$AWC_CLASS[which(records$MU_GLOBAL==11254 & records$T_USDA_TEX_CLASS==9)],"\n")
  # replace AWD codes with values
  #m <- dbDriver("SQLite")
  #con <- dbConnect(m,dbname="HWSD.sqlite")
  #awdcode <- dbGetQuery(con, "select * from D_AWC")
  #dbDisconnect(con)
  # version 2.0 or higher, use CSV extration of D_AWC table
  awdcode <- read.csv("/p/projects/macmit/data/GGCMI/phase3/hwsd/hwsd_awc.csv")
  for(i in awdcode$CODE) records$AWC_CLASS[records$AWC_CLASS==i] <- awdcode$VALUE[i]
  
  
  # fill missing values in ISIMIP landmask with nearest neighbor HWSD mu_global
  # do this before assinging all other parameters based on MU_GLOBAL
  
  # cropland weighted selection
  if(do.cropland.weighting){
    missing <- landmask
    missing[is.finite(data_mu)] <- NA
    miss.index <- which(is.finite(missing),arr.ind=T)
    nearest.neighbor <- find.dominant.neighbor(data_mu,miss.index)
    data_issoil <- data_mu
    data_issoil[is.finite(data_mu)] <- 1
    for(i in 1:dim(miss.index)[1])
    {
      # assign 0 to issoil variable for soils that have been extrapolated
      data_issoil[miss.index[i,1],miss.index[i,2]] <- 0
      data_mu[miss.index[i,1],miss.index[i,2]] <- nearest.neighbor[i]
      if(!is.finite(data_tex_class[miss.index[i,1],miss.index[i,2]])){
        tc <- records$T_USDA_TEX_CLASS[which(records$MU_GLOBAL==nearest.neighbor[i])]
        data_tex_class[miss.index[i,1],miss.index[i,2]] <- tc[which(tc>0)][1]
      } else {
        cat("missing MU but exisitig TEX_CLASS",data_tex_class[miss.index[i,1],miss.index[i,2]],"in pixel",i,"\n")
      }
    }
    data_ph <- data_co <- data_bd <- data_cec <- data_oc <- 
      data_ro <- data_il <- data_awc <- data_sand <- data_silt <- 
      data_clay <- data_gravel <- data_ece <- data_bs <- data_mu
    muunique <- unique(data_mu[is.finite(data_mu)])
    #cat("data_mu\n",range(data_mu),range(data_mu,na.rm=T),"\n")
    #str(data_mu)
    #str(muunique)
    #for(i in 1:length(records$MU_GLOBAL)){
    for(i in muunique){
      index <- which(records$MU_GLOBAL==i)
      if(i%%100==0) cat("processing",i,"for",index,"\n")
      if(length(index>1)){
        for(j in 1:length(index)){
          index2 <- index[j]
          data_ph[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_PH_H2O[index2]
          data_co[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_CACO3[index2]
          data_bd[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_BULK_DENSITY[index2]
          data_cec[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_CEC_SOIL[index2]
          data_oc[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_OC[index2]
          data_ro[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$ROOTS[index2]
          data_il[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$IL[index2]
          data_awc[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$AWC_CLASS[index2]
          data_sand[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_SAND[index2]
          data_silt[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_SILT[index2]
          data_clay[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_CLAY[index2]
          data_gravel[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_GRAVEL[index2]
          data_ece[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_ECE[index2]
          data_bs[data_mu==i & data_tex_class==records$T_USDA_TEX_CLASS[index2]] <- records$T_BS[index2]
        }
        
      } else {
        data_ph[data_mu==i] <- records$T_PH_H2O[index]
        data_co[data_mu==i] <- records$T_CACO3[index]
        data_bd[data_mu==i] <- records$T_BULK_DENSITY[index]
        data_cec[data_mu==i] <- records$T_CEC_SOIL[index]
        data_oc[data_mu==i] <- records$T_OC[index]
        data_ro[data_mu==i] <- records$ROOTS[index]
        data_il[data_mu==i] <- records$IL[index]
        data_awc[data_mu==i] <- records$AWC_CLASS[index]
        data_sand[data_mu==i] <- records$T_SAND[index]
        data_silt[data_mu==i] <- records$T_SILT[index]
        data_clay[data_mu==i] <- records$T_CLAY[index]
        data_gravel[data_mu==i] <- records$T_GRAVEL[index]
        data_ece[data_mu==i] <- records$T_ECE[index]
        data_bs[data_mu==i] <- records$T_BS[index]
      }
    }
    
    
    # write netcdf file ####  
    dim_lon <- ncdim_def("lon","degrees_east",seq(-179.75,179.75,len=360/0.5))
    #change order of latitudes
    dim_lat <- ncdim_def("lat","degrees_north",seq(89.75,-89.75,len=180/0.5))
    mv <- 1e20
    ncv_st <-ncvar_def("texture_class","-",list(dim_lon,dim_lat),mv,
                       longname="USDA soil texture class dominant HWSD on cropland",
                       compression=6)
    ncv_mu <-ncvar_def("mu_global","-",list(dim_lon,dim_lat),mv,
                       longname="dominant HWSD soil mapping unit within dominant USDA soil texture class on cropland",
                       compression=6)
    ncv_ph <-ncvar_def("soil_ph","-",list(dim_lon,dim_lat),mv,
                       longname="Topsoil pH(H2O)",
                       compression=6)
    ncv_co <-ncvar_def("soil_caco3","percent weight",list(dim_lon,dim_lat),mv,
                       longname="Topsoil Calcium Carbonate",
                       compression=6)
    ncv_bd <-ncvar_def("bulk_density","kg dm-1",list(dim_lon,dim_lat),mv,
                       longname="Topsoil Bulk Density",
                       compression=6)
    ncv_cec <-ncvar_def("cec_soil","cmol kg-1",list(dim_lon,dim_lat),mv,
                        longname="Topsoil Cation Exchange Capacity (soil)",
                        compression=6)
    ncv_oc <-ncvar_def("oc","percent weight",list(dim_lon,dim_lat),mv,
                       longname="Topsoil Organic  Carbon",
                       compression=6)
    ncv_ro <-ncvar_def("root_obstacles","- (class)",list(dim_lon,dim_lat),mv,
                       longname="depth of Obstacles to Roots (ESDB) (class6: 0cm - 20cm, class5: 0cm - 80cm, class4: 20cm - 40cm, class3: 40cm - 60cm, class2: 60cm - 80cm, class1: >80cm",
                       compression=6)
    ncv_il <-ncvar_def("impermeable_layer","- (class)",list(dim_lon,dim_lat),mv,
                       longname="depth of Impermeable Layer (ESDB), (class4: <40cm, class3: 40cm - 80cm, class2: 80cm - 150cm, class1: >150cm)",
                       compression=6)
    ncv_awc <-ncvar_def("awc","mm",list(dim_lon,dim_lat),mv,
                        longname="Available Water Content",
                        compression=6)
    ncv_sand <-ncvar_def("sand","percent",list(dim_lon,dim_lat),mv,
                         longname="Topsoil Sand Fraction",
                         compression=6)
    ncv_clay <-ncvar_def("clay","percent",list(dim_lon,dim_lat),mv,
                         longname="Topsoil Clay Fraction",
                         compression=6)
    ncv_silt <-ncvar_def("silt","percent",list(dim_lon,dim_lat),mv,
                         longname="Topsoil Silt Fraction",
                         compression=6)
    ncv_gravel <-ncvar_def("gravel","percent",list(dim_lon,dim_lat),mv,
                           longname="Topsoil Gravel Content",
                           compression=6)
    ncv_ece <-ncvar_def("ece","dS m-1",list(dim_lon,dim_lat),mv,
                        longname="Topsoil Salinity (ECe)",
                        compression=6)
    ncv_bs <-ncvar_def("bs_soil","percent weight",list(dim_lon,dim_lat),mv,
                       longname="Topsoil Base Saturation",
                       compression=6)
    ncv_issoil <-ncvar_def("issoil","",list(dim_lon,dim_lat),mv,
                           longname="flag for valid soils",
                           compression=6)
    
    
    # versions < 1.0 are for internal discussions only
    # version 1.0 was released to the group to work with these data
    # version 1.01 includes also the ISSOIL variable, otherwise unchanged
    # version 2.0 uses mdbtools to convert the HWSD_DATA table to a csv table for propper treatment of missing values
    # version 2.1 with issoil=0 for gap-filled pixels
    # version 2.2 with hand-filled missing values for bulk density (2 MU_GLOBAL), and AWC ()
    # version 2.3 corrected units for impermeable_layer and root_obstacles
    
    ncf <- nc_create("HWSD_soil_data_on_cropland_v2.3.nc",
                     list(ncv_st,ncv_mu,ncv_ph,ncv_co,ncv_bd,ncv_cec,ncv_oc,
                          ncv_ro,ncv_il,ncv_awc,ncv_sand,ncv_silt,ncv_clay,ncv_gravel,
                          ncv_ece,ncv_bs,ncv_issoil))
    ncatt_put(ncf,varid=0,"author","cmueller@pik-potsdam.de")
    ncatt_put(ncf,varid=0,"version","Version 2.3, proper treating for missing values in processing HWSD, 2.1 issoi=0 for gap filled pixles, 2.2 hand-filled missing values for BD and AWC, 2.3 corrected units for impermeable_layer and root_obstacles")
    ncatt_put(ncf,varid=0,"comment","GGCMI Phase 3 soil input data set for usage in ISIMIP/GGCMI Phase 3 simulations, data aggregated by dominant soil profile (MU_GLOBAL) within dominant soil texture class from HWSD on current cropland (MIRCA2000 at 5 arc-minutes)")
    
    ncvar_put(ncf,ncv_ph,data_ph[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_st,data_tex_class[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_mu,data_mu[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_co,data_co[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_bd,data_bd[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_cec,data_cec[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_oc,data_oc[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_ro,data_ro[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_il,data_il[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_awc,data_awc[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_sand,data_sand[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_silt,data_silt[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_clay,data_clay[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_gravel,data_gravel[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_ece,data_ece[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_bs,data_bs[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_issoil,data_issoil[,360:1],start=c(1,1),count=c(-1,-1))
    
    nc_close(ncf)
  }

  # all land selection
  if(do.allland){
    missing <- landmask
    missing[is.finite(data_mu_nolu)] <- NA
    miss.index <- which(is.finite(missing),arr.ind=T)
    nearest.neighbor <- find.dominant.neighbor(data_mu_nolu,miss.index)
    data_issoil <- data_mu_nolu
    data_issoil[is.finite(data_mu_nolu)] <- 1
    for(i in 1:dim(miss.index)[1])
    {
      # assign 0 to issoil variable for soils that have been extrapolated
      data_issoil[miss.index[i,1],miss.index[i,2]] <- 0
      data_mu_nolu[miss.index[i,1],miss.index[i,2]] <- nearest.neighbor[i]
      if(!is.finite(data_tex_class_nolu[miss.index[i,1],miss.index[i,2]])){
        tc <- records$T_USDA_TEX_CLASS[which(records$MU_GLOBAL==nearest.neighbor[i])]
        data_tex_class_nolu[miss.index[i,1],miss.index[i,2]] <- tc[which(tc>0)][1]
      } else {
        cat("missing MU but exisitig TEX_CLASS",data_tex_class_nolu[miss.index[i,1],miss.index[i,2]],"in pixel",i,"\n")
      }
    }
    data_ph <- data_co <- data_bd <- data_cec <- data_oc <- 
      data_ro <- data_il <- data_awc <- data_sand <- data_silt <- 
      data_clay <- data_gravel <- data_ece <- data_bs <- data_mu_nolu
    muunique <- unique(data_mu_nolu[is.finite(data_mu_nolu)])
    for(i in muunique){
      index <- which(records$MU_GLOBAL==i)
      if(i%%100==0) cat("processing",i,"for",index,"\n")
      if(length(index>1)){
        for(j in 1:length(index)){
          index2 <- index[j]
          data_ph[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_PH_H2O[index2]
          data_co[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_CACO3[index2]
          data_bd[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_BULK_DENSITY[index2]
          data_cec[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_CEC_SOIL[index2]
          data_oc[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_OC[index2]
          data_ro[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$ROOTS[index2]
          data_il[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$IL[index2]
          data_awc[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$AWC_CLASS[index2]
          data_sand[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_SAND[index2]
          data_silt[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_SILT[index2]
          data_clay[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_CLAY[index2]
          data_gravel[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_GRAVEL[index2]
          data_ece[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_ECE[index2]
          data_bs[data_mu_nolu==i & data_tex_class_nolu==records$T_USDA_TEX_CLASS[index2]] <- records$T_BS[index2]
        }
        
      } else {
        data_ph[data_mu_nolu==i] <- records$T_PH_H2O[index]
        data_co[data_mu_nolu==i] <- records$T_CACO3[index]
        data_bd[data_mu_nolu==i] <- records$T_BULK_DENSITY[index]
        data_cec[data_mu_nolu==i] <- records$T_CEC_SOIL[index]
        data_oc[data_mu_nolu==i] <- records$T_OC[index]
        data_ro[data_mu_nolu==i] <- records$ROOTS[index]
        data_il[data_mu_nolu==i] <- records$IL[index]
        data_awc[data_mu_nolu==i] <- records$AWC_CLASS[index]
        data_sand[data_mu_nolu==i] <- records$T_SAND[index]
        data_silt[data_mu_nolu==i] <- records$T_SILT[index]
        data_clay[data_mu_nolu==i] <- records$T_CLAY[index]
        data_gravel[data_mu_nolu==i] <- records$T_GRAVEL[index]
        data_ece[data_mu_nolu==i] <- records$T_ECE[index]
        data_bs[data_mu_nolu==i] <- records$T_BS[index]
      }
    }
    
    
    # write netcdf file ####  
    dim_lon <- ncdim_def("lon","degrees_east",seq(-179.75,179.75,len=360/0.5))
    #change order of latitudes
    dim_lat <- ncdim_def("lat","degrees_north",seq(89.75,-89.75,len=180/0.5))
    mv <- 1e20
    ncv_st <-ncvar_def("texture_class","-",list(dim_lon,dim_lat),mv,
                       longname="USDA soil texture class dominant HWSD on all land",
                       compression=6)
    ncv_mu <-ncvar_def("mu_global","-",list(dim_lon,dim_lat),mv,
                       longname="dominant HWSD soil mapping unit within dominant USDA soil texture class on all land",
                       compression=6)
    ncv_ph <-ncvar_def("soil_ph","-",list(dim_lon,dim_lat),mv,
                       longname="Topsoil pH(H2O)",
                       compression=6)
    ncv_co <-ncvar_def("soil_caco3","percent weight",list(dim_lon,dim_lat),mv,
                       longname="Topsoil Calcium Carbonate",
                       compression=6)
    ncv_bd <-ncvar_def("bulk_density","kg dm-1",list(dim_lon,dim_lat),mv,
                       longname="Topsoil Bulk Density",
                       compression=6)
    ncv_cec <-ncvar_def("cec_soil","cmol kg-1",list(dim_lon,dim_lat),mv,
                        longname="Topsoil Cation Exchange Capacity (soil)",
                        compression=6)
    ncv_oc <-ncvar_def("oc","percent weight",list(dim_lon,dim_lat),mv,
                       longname="Topsoil Organic  Carbon",
                       compression=6)
    ncv_ro <-ncvar_def("root_obstacles","- (class)",list(dim_lon,dim_lat),mv,
                       longname="depth of Obstacles to Roots (ESDB) (class6: 0cm - 20cm, class5: 0cm - 80cm, class4: 20cm - 40cm, class3: 40cm - 60cm, class2: 60cm - 80cm, class1: >80cm",
                       compression=6)
    ncv_il <-ncvar_def("impermeable_layer","- (class)",list(dim_lon,dim_lat),mv,
                       longname="depth of Impermeable Layer (ESDB), (class4: <40cm, class3: 40cm - 80cm, class2: 80cm - 150cm, class1: >150cm)",
                       compression=6)
    ncv_awc <-ncvar_def("awc","mm",list(dim_lon,dim_lat),mv,
                        longname="Available Water Content",
                        compression=6)
    ncv_sand <-ncvar_def("sand","percent",list(dim_lon,dim_lat),mv,
                         longname="Topsoil Sand Fraction",
                         compression=6)
    ncv_clay <-ncvar_def("clay","percent",list(dim_lon,dim_lat),mv,
                         longname="Topsoil Clay Fraction",
                         compression=6)
    ncv_silt <-ncvar_def("silt","percent",list(dim_lon,dim_lat),mv,
                         longname="Topsoil Silt Fraction",
                         compression=6)
    ncv_gravel <-ncvar_def("gravel","percent",list(dim_lon,dim_lat),mv,
                           longname="Topsoil Gravel Content",
                           compression=6)
    ncv_ece <-ncvar_def("ece","dS m-1",list(dim_lon,dim_lat),mv,
                        longname="Topsoil Salinity (ECe)",
                        compression=6)
    ncv_bs <-ncvar_def("bs_soil","percent weight",list(dim_lon,dim_lat),mv,
                       longname="Topsoil Base Saturation",
                       compression=6)
    ncv_issoil <-ncvar_def("issoil","",list(dim_lon,dim_lat),mv,
                           longname="flag for valid soils",
                           compression=6)
    
    
    ncf <- nc_create("HWSD_soil_data_all_land_v2.3.nc",
                     list(ncv_st,ncv_mu,ncv_ph,ncv_co,ncv_bd,ncv_cec,ncv_oc,
                          ncv_ro,ncv_il,ncv_awc,ncv_sand,ncv_silt,ncv_clay,ncv_gravel,
                          ncv_ece,ncv_bs,ncv_issoil))
    ncatt_put(ncf,varid=0,"author","cmueller@pik-potsdam.de")
    ncatt_put(ncf,varid=0,"version","Version 2.3, proper treating for missing values in processing HWSD, 2.1 issoi=0 for gap filled pixles, 2.2 hand-filled missing values for BD and AWC, 2.3 corrected units for impermeable_layer and root_obstacles")
    ncatt_put(ncf,varid=0,"comment","GGCMI Phase 3 soil input data set for usage in ISIMIP/GGCMI Phase 3 simulations, data aggregated by dominant soil profile (MU_GLOBAL) within dominant soil texture class from HWSD on all land")
    
    ncvar_put(ncf,ncv_ph,data_ph[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_st,data_tex_class_nolu[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_mu,data_mu_nolu[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_co,data_co[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_bd,data_bd[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_cec,data_cec[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_oc,data_oc[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_ro,data_ro[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_il,data_il[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_awc,data_awc[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_sand,data_sand[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_silt,data_silt[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_clay,data_clay[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_gravel,data_gravel[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_ece,data_ece[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_bs,data_bs[,360:1],start=c(1,1),count=c(-1,-1))
    ncvar_put(ncf,ncv_issoil,data_issoil[,360:1],start=c(1,1),count=c(-1,-1))
    
    nc_close(ncf)
  }
  
}

