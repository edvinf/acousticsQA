library(RstoxData)
library(ggplot2)


#' Horisontal Profile
#'
#' @description
#'  sa in horisontal bands of the water column, along with their distance from seadbed (bottom)
#'
#' @details
#'  data table with columns:
#'  \describe{
#'   \item{height}{Factor for each horisontal band. Names indicate interval in meters.}
#'   \item{sa}{Aggregated sa in horsintal band}
#'  }
#' @name horisontalProfile
NULL


#' Extract Horisontal profile
#' @description
#'  Makes horisontal profile (sa by approximate distance from bottom) from pelagic channels
#' @details
#'  Consult reference table for acoustic categories.
#'  https://referenceeditor.hi.no/apps/referenceeditor/v2/tables/acousticCategory
#'
#'  some common options:
#'  1 (other)
#'  2
#'  6
#'  12 (herring)
#'  22 (saithe)
#'  30 (haddock)
#'  31 (cod)
#'
#'  Note that acoustic categories are not uniquely specified by species,
#'  consult survey manual and inspect survey data to make sure you cover them all.
#'
#' @param echosounder luf20 as read by RstoxData::readXmlFile
#' @param targetAcocat acoustic category to extract profile for
#' @param bins number of bins for the profile
#' @param freq frequncy of transreceiver to extract data from.
#' @return \code{\link[acousticsQA]{horisontalProfile}}
#' @export
horisontalProfileLUF20 <- function(echosounder, targetAcocat, bins=20, freq=38000){
  target <- echosounder$sa_by_acocat[echosounder$sa_by_acocat$acocat %in% targetAcocat,]
  targetPC <- target[target$type=="P",]
  targetPC38 <- targetPC[targetPC$freq==freq,]
  targetPC38A <- merge(targetPC38, echosounder$frequency, all.x=T)
  targetSA <- merge(targetPC38A, echosounder$sa, all.x=T, by=c("log_start", "start_time", "freq", "transceiver", "type", "acocat"))
  targetSA <- merge(targetSA, echosounder$distance, all.x=T, by=c("log_start", "start_time"))
  stopifnot(length(unique(targetSA$transceiver))==1)
  targetSA$height <- targetSA$max_bot_depth - targetSA$ch*targetSA$pel_ch_thickness
  targetSA$heightbins <- cut(targetSA$height, breaks=bins, ordered_result = T)
  profile <- aggregate(list(sa=targetSA$sa), by=list(height=targetSA$heightbins), FUN=sum)

  return(profile)
}

#' Vertical Profile
#'
#' @description
#'  sa in vertical bands of the water column, along with log (sailed distance),
#'  position parameters and time parameters
#'
#' @details
#'  data table with columns:
#'  \describe{
#'   \item{sa}{Aggregated sa in vertical band}
#'   \item{log}{log, sailed distance in nmi}
#'   \item{time}{time as POSIXct}
#'   \item{latitude}{latitude in decimal degrees, WGS84}
#'   \item{longitude}{longitude in decimal degrees, WGS84}
#'   \item{distance}{width of vertical band in nmi}
#'  }
#' @name horisontalProfile
NULL

#' Vertical profile
#' @description
#'  Makes vertical profile (sa by sailed distance)
#' @details
#'  Consult reference table for acoustic categories.
#'  https://referenceeditor.hi.no/apps/referenceeditor/v2/tables/acousticCategory
#'
#'  some common options:
#'  1 (other)
#'  2
#'  6
#'  12 (herring)
#'  22 (saithe)
#'  30 (haddock)
#'  31 (cod)
#'
#' @param echosounder luf20 as read by RstoxData::readXmlFile
#' @param targetAcocat acoustic category to extract profile for
#' @param channelType the channel type to extract profile from (e.g. Pelagic (P) or Bottom (B))
#' @param freq frequncy of transreceiver to extract data from.
#' @return \code{\link[acousticsQA]{verticalProfile}}
#' @export
verticalProfileLUF20 <- function(echosounder, targetAcocat, channelType="P", freq=38000){
  target <- echosounder$sa_by_acocat[echosounder$sa_by_acocat$acocat %in% targetAcocat,]
  targetPC <- target[target$type==channelType,]
  targetPC38 <- targetPC[targetPC$freq==freq,]
  targetSA <- merge(targetPC38, echosounder$sa, all.x=T, by=c("log_start", "start_time", "freq", "transceiver", "type", "acocat"))
  stopifnot(length(unique(targetSA$freq))==1)
  stopifnot(length(unique(targetSA$transceiver))==1)
  stopifnot(length(unique(targetSA$type))==1)
  targetSAInt <- targetSA[, list(sa=sum(get("sa"))), by=list(log_start=get("log_start"), start_time=get("start_time"))]
  dist <- echosounder$distance[,c("log_start", "start_time", "lat_start", "lon_start", "integrator_dist"), with=F]
  targetSAInt <- merge(targetSAInt, dist, all.y=T, by=c("log_start", "start_time"))
  targetSAInt$start_time <- as.POSIXct(targetSAInt$start_time)
  targetSAInt <- targetSAInt[order(targetSAInt$start_time),]
  targetSAInt <- targetSAInt[,c("sa", "log_start", "start_time", "lat_start", "lon_start", "integrator_dist")]
  names(targetSAInt) <- c("sa", "log", "time", "latitude", "longitude", "distance")


  return(targetSAInt)
}

#' Trawl locations
#' @description
#'  Trawl locations
#' @details
#'  \code{\link[data.table]{data.table}} with columns
#'  \describe{
#'   \item{log}{log, sailed distance in nmi}
#'   \item{time}{time as POSIXct}
#'   \item{latitude}{latitude in decimal degrees, WGS84}
#'   \item{longitude}{longitude in decimal degrees, WGS84}
#'  }
#'
#' @name trawlLocation
#'
NULL

#' Extract Trawl
#' @description
#'  Extract trawl locations from Biotic data
#' @param biotic NMDbiotic data (v3.x) as parsed by RstoxData::readXmlFile
#' @param serialnumbers serialnumbers to include. If NULL all are included. Compared to column (fishstation/serialnumber)
#' @param sampleQuality sample quality codes to include. If NULL all are included. Compared to column (fishstation/samplequality)
#' @param serialnumbers gear codes to include. If NULL all are included. Compared to column (fishstation/gear)
#' @return \code{\link[acousticsQA]{trawlLocation}}
#' @export
extractTrawlsBiotic <- function(biotic, serialnumbers=NULL, sampleQuality=NULL, gearCodes=NULL){

  st <- biotic$fishstation
  if (!is.null(serialnumbers)){
    st <- st[st$serialnumber %in% serialnumbers]
  }
  if (!is.null(sampleQuality)){
    st <- st[st$samplequality %in% sampleQuality]
  }
  if (!is.null(gearCodes)){
    st <- st[st$gear %in% gearCodes]
  }

  st$time <- as.POSIXct(paste(st$stationstartdate, st$stationstarttime))
  st <- st[,c("logstart", "time", "latitudestart", "longitudestart")]
  names(st) <- c("log", "time", "latitude", "longitude")

  return(st)
}

#' Plots horisontal profile
#' @description
#'  Plots aggregated sa binned by distance from bottom
#' @details
#'  Plots aggregated sa binned by distance from bottom
#'  Cumulative fractions are plotted on an alternate y-axis
#' @param profile \code{\link[acousticsQA]{horisontalProfile}}
#' @param header header for plot
#' @export
plotHorisontalProfile <- function(profile, header=""){
  profile$sa2 <- max(profile$sa)*profile$sa/sum(profile$sa)
  coeff <- 1/sum(profile$sa2)
  ggplot2::ggplot(profile, ggplot2::aes(x=height)) +
    ggplot2::geom_col(ggplot2::aes(y=sa)) +
    ggplot2::geom_point(ggplot2::aes(y=cumsum(sa2)), group=1) +
    ggplot2::scale_y_continuous(

      # Features of the first axis
      name = "total sA",

      # Add a second axis and specify its features
      sec.axis = ggplot2::sec_axis(~.*coeff, name="cum Frac sA")
    ) +
    ggplot2::ggtitle(header) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))

}

#' Plot vertical profile
#' @description
#'  Plots aggregated sa for each vertical band
#' @details
#'  Plots aggregated sa for each vertical band
#'  Cumulative fractions are plotted on an alternate y-axis
#'  Trawl locations are plotted as horisontal lines if provided.
#' @param profile \code{\link[acousticsQA]{verticalProfile}}
#' @param header header for plot
#' @param trawls trawl locations, \code{\link[acousticsQA]{trawlLocation}}. If NULL, trawl locations are not plotted.
#' @export
plotStretch <- function(profile, header="", trawls=NULL){

  profile$sa[is.na(profile$sa)] <- 0

  profile$sa2 <- max(profile$sa)*profile$sa/sum(profile$sa)
  coeff <- 1/sum(profile$sa2)
  pl <- ggplot2::ggplot(profile) +
    ggplot2::geom_col(ggplot2::aes_string(x="log", y="sa"), width = .99) +
    ggplot2::geom_line(ggplot2::aes(x=log, y=cumsum(sa2))) +
    ggplot2::scale_y_continuous(

      # Features of the first axis
      name = "sA",

      # Add a second axis and specify its features
      sec.axis = ggplot2::sec_axis(~.*coeff, name="cum Frac sA")
    ) +
    ggplot2::xlab("log (nmi)") +
    ggplot2::ggtitle(header) +
    ggplot2::theme_bw()

  if (!is.null(trawls)){
    trawls$cumSa <- as.numeric(NA)
    for (i in 1:nrow(trawls)){
      trawls$cumSa[i] <- sum(profile$sa[profile$log<=trawls$log[i]]) * max(profile$sa) / sum(profile$sa)
    }
    pl <- pl + ggplot2::geom_point(data=trawls, ggplot2::aes_string(x="log", y="cumSa"))
    #pl <- pl + ggplot2::geom_hline(data=trawls, ggplot2::aes_string(yintercept="cumSa"))
  }

  pl
}

#' Plot map
#' @description
#'  Plots map of acoustic registration, and optionally trawl positions.
#'  Vertical bands in acoustic registrations are represented as points with area proportional to registration (sa).
#' @param profile \code{\link[acousticsQA]{verticalProfile}} from acoustic registrations
#' @param trawls \code{\link[acousticsQA]{trawlLocation}}
#' @param lonLim limits for plots (longitudes). Vector of size 2. Will be calculated from 'profile' if NULL.
#' @param latLim limits for plots (latitudes). Vector of size 2. Will be calculated from 'profile' if NULL.
#' @param projection proj4string or EPSG code specifying the desired projection, see \code{\link[sf]{st_crs}}. Defaults to mercator projection.
#' @param header header for plot
#' @param maxSaSize controls the maximum size of points representing vertical bands in 'profile'
#' @param saColor color of points representing vertical bands in 'profile'
#' @param trawlColor color of points representing trawl locations.
#' @param trawlPointShape code for point shape for trawl locations (ggplot2-codes)
#' @export
plotMap <- function(profile, trawls=NULL, lonLim=NULL, latLim=NULL, projection=NULL, header="", maxSaSize=5, saColor="black", trawlColor="red", trawlPointShape=10){

  profile$sa[is.na(profile$sa)] <- 0

  if (is.null(lonLim)){
    lonLim <- c(min(profile$longitude), max(profile$longitude))
  }
  if (is.null(latLim)){
    latLim <- c(min(profile$latitude), max(profile$latitude))
  }
  if (is.null(projection)) {
    projection <- "+proj=merc +datum=WGS84"
  }
  newcrs <- sf::st_crs(projection)

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  limboxP <- sf::st_bbox(sf::st_transform(sf::st_sfc(sf::st_point(x = c(min(lonLim), min(latLim))),
                                                     sf::st_point(x = c(max(lonLim), min(latLim))),
                                                     sf::st_point(x = c(min(lonLim), max(latLim))),
                                                     sf::st_point(x = c(max(lonLim), max(latLim))),
                                                     sf::st_point(x = c(min(lonLim), mean(latLim))),
                                                     sf::st_point(x = c(max(lonLim), mean(latLim))),
                                                     sf::st_point(x = c(mean(lonLim), min(latLim))),
                                                     sf::st_point(x = c(mean(lonLim), max(latLim))),crs = sf::st_crs(4326)), newcrs))

  points <- sf::st_as_sf(profile, coords=c("longitude","latitude"), crs=sf::st_crs(4326))

  pl <- ggplot2::ggplot(points) +
    ggplot2::geom_sf(data=world) +
    ggplot2::geom_sf(data=points, mapping=ggplot2::aes_string(size="sa"), shape=21, alpha = 0.7, colour = "black",fill=saColor,stroke = .2) +
    ggplot2::scale_size_area(max_size=maxSaSize) +
    ggplot2::ggtitle(header) +
    ggplot2::theme_bw()

  if (!is.null(trawls)){
    trawlPoints <- sf::st_as_sf(trawls, coords=c("longitude","latitude"), crs=sf::st_crs(4326))
    pl <- pl + ggplot2::geom_sf(data=trawlPoints, shape=trawlPointShape, alpha = 0.7, colour = trawlColor, fill=trawlColor)
  }

  pl <- pl + ggplot2::coord_sf(crs = newcrs, xlim = c(limboxP$xmin, limboxP$xmax), ylim = c(limboxP$ymin, limboxP$ymax), expand = T)
  pl
}

#' Plot concentration
#' @description
#'  Plots cumulative sa against cumulative distance,
#'  in order of decreasing sa
#' @param profile \code{\link[acousticsQA]{verticalProfile}}
#' @param header header for plot
#' @export
plotConcentration <- function(profile, header=""){
  profile$sa[is.na(profile$sa)] <- 0
  profile <- profile[order(profile$sa, decreasing = T),]
  profile$cumSum <- cumsum(profile$sa)
  profile$cumFrac <- profile$cumSum / sum(profile$sa)
  profile$distance <- cumsum(profile$distance)
  ggplot2::ggplot(profile) +
    ggplot2::geom_path(ggplot2::aes_string(x="distance", y="cumSum")) +
    ggplot2::ylab("Cum. sa") +
    ggplot2::xlab("Cum. distance (nmi)") +
    ggplot2::ggtitle(header) +
    ggplot2:: theme_bw()

}

#' Plot purity
#' @description
#'  Plots contributions to total echo by the purity of target acoustic category.
#' @details
#'  Each interpreted segment is assigned a purity, which is the fraction of the total
#'  echo assigned to the acoustic category of interest (target acoustic category).
#'  The cumulative contribution to the total echo of the target acoustic category is plotted against
#'  different thresholds for purity. In the example provided below, one can for example infer that between
#'  80\% and 90\% of the total echo assigned to cod came from segments where at least 80 % of the echo was
#'  assigned to COD.
#'
#'  Purity of assignment in a segment may serve as a proxy for confidence in the assignment,
#'  it typically reflects that categories are well separable acoustically, and that trawl-samples supports
#'  a relatively pure assignment.
#'
#' @param profileTarget \code{\link[acousticsQA]{verticalProfile}} from acoustic registrations of target acoustic category
#' @param profileAll \code{\link[acousticsQA]{verticalProfile}} from acoustic registrations of all acoustic categories
#' @param header header for plot
#' @param breaks vector of purities ([0,1]) to plot contributions for
#' @param ylim limits for y-axis.
#' @examples
#'  prof <- verticalProfileLUF20(acousticsQA::echosounderSkrei2019, 31)
#'  profAll <- verticalProfileLUF20(acousticsQA::echosounderSkrei2019,
#'      unique(acousticsQA::echosounderSkrei2019$acocat$acocat))
#'  plotPurity(prof, profAll, header="COD")
#' @export
plotPurity <- function(profileTarget, profileAll, header="", breaks=seq(0.5,1,.05), ylim=c(0,1)){
  profileAll$sa[is.na(profileAll$sa)] <- 0
  profileTarget$sa[is.na(profileTarget$sa)] <- 0

  tot <- profileAll[,c("log", "time", "sa")]
  names(tot) <- c("log", "time", "saTotal")
  profile <- merge(profileTarget, tot)
  profile$purity <- profile$sa / profile$saTotal

  purity <- data.table::data.table(purity=breaks)
  purity$contribution <- as.numeric(NA)

  for (i in 1:nrow(purity)){
    purity$contribution[i] <- sum(profile$sa[profile$saTotal > 0 & profile$purity >= purity$purity[i]]) / sum(profile$sa)
  }

  ggplot2::ggplot(purity) +
    ggplot2::geom_col(ggplot2::aes_string(x="purity", y="contribution")) +
    ggplot2::xlab("purity threshold") +
    ggplot2::ylab("contribution") +
    ggplot2::ggtitle(header) +
    ggplot2::ylim(ylim) +
    ggplot2::theme_bw()
}

#' Plot lognormal QQ plot
#' @description
#'  Plot quantiles of log-transformed sa against quantiles of normal distribution.
#' @param profile \code{\link[acousticsQA]{verticalProfile}}, NAs are ignored.
#' @param header header for plot
#' @export
plotLogNormalQQ <- function(profile, header=""){
  profile <- profile[!is.na(profile$sa),]

  profile$logsa <- log(profile$sa)
  ggplot2::ggplot(profile, ggplot2::aes_string(sample="logsa")) +
    ggplot2::stat_qq() +
    ggplot2::geom_qq_line() +
    ggplot2::ylab("ln sa quantiles") +
    ggplot2::xlab("std normal quantiles") +
    ggplot2::ggtitle(header) +
    ggplot2::theme_bw()
}

#' Plot histogram of ln sa
#' @description
#'  Plot distribution of ln sa over vertical bands as histogram
#' @param profile \code{\link[acousticsQA]{verticalProfile}}, NAs are ignored.
#' @param header header for plot
#' @export
plotLogSaHistogram <- function(profile, header=""){
  profile <- profile[!is.na(profile$sa),]
  profile$logsa <- log(profile$sa)

  ggplot2::ggplot(profile, ggplot2::aes_string(sample="logsa")) +
    ggplot2::geom_histogram(ggplot2::aes_string(x="logsa"), bins = 30) +
    ggplot2::xlab("ln sa") +
    ggplot2::ggtitle(header) +
    ggplot2::theme_bw()
}

#' @noRd
plotDistanceTrawl <- function(profile, trawls, header=""){
  profile$trawlDist <- as.numeric(NA)
  for (i in 1:nrow(profile)){
    minDist <- NA
    for (j in 1:nrow(trawls)){
      dist <- geosphere::distHaversine(c(profile$longitude[i], profile$latitude[i]), c(trawls$longitude[j], trawls$latitude[j])) / 1852
      minDist <- min(minDist, dist, na.rm=T)
    }
    profile$trawlDist[i] <- minDist
  }
  profile$sa[is.na(profile$sa)] <- 0
  ggplot2::ggplot(profile) +
    ggplot2::geom_point(ggplot2::aes_string(y="sa", x="trawlDist")) +
    ggplot2::ggtitle(header) +
    ggplot2::ylab("allocated sa") +
    ggplot2::xlab("nearest trawl (nmi)") +
    ggplot2::theme_bw()

}

# plot sa vs distance to trawl

# add fraction of total fish in 'plotConcentration' as color ?

# make variant with 2-d grid mot mercator proj coordinates and survey trace, support adding trawled positions
