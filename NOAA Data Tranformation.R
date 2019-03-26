setwd("D:/University of Wollongong/Data Mining project")

# Load read netCDF file function
library(ncdf4)

# Import file into environment
filePath <- "D:/University of Wollongong/Data Mining project/"
fileName <- "air.mon.anom"
file <- paste(filePath, fileName, ".nc", sep = "")
raw_data <- nc_open(file)

# Get data structure
print(raw_data)

# Get location points
lat <- ncvar_get(raw_data, "lat")
lon <- ncvar_get(raw_data, "lon")
no_lat <- dim(lat)
no_lon <- dim(lon)
coordinates <- as.matrix(expand.grid(lat, lon))

# Get time and its units
time <- ncvar_get(raw_data, "time")
no_time <- dim(time)
time_units <- ncatt_get(raw_data, "time", "units")

# Get air temperature
air_temp <- ncvar_get(raw_data, "air")
data_long_name <- ncatt_get(raw_data, "air", "long_name")
temp_units <- ncatt_get(raw_data, "air", "units")
fillValue <- ncatt_get(raw_data, "air", "valid_range")

# Get global attributes
title <- ncatt_get(raw_data, 0, "title")
source <- ncatt_get(raw_data, 0, "Source")
history <- ncatt_get(raw_data, 0, "history")
references <- ncatt_get(raw_data, 0, "References")

# Close the file to keep original dataset safe
nc_close(raw_data)

# Load some more necessary libraries
library(chron)

# Convert time into readable format
time_units_string <- strsplit(time_units$value, " ")
time_date_split <- strsplit(unlist(time_units_string)[3], "-")
day <- as.integer(unlist(time_date_split)[3])
month <- as.integer(unlist(time_date_split)[2])
year <- as.integer(unlist(time_date_split)[1])
time <- chron(time, format = c(dates = "dd/mm/yyyyy"), origin = c(day, month, year))

# Get vectors of month indexes
month <- list(vector(), vector(), vector(), vector(), vector(), vector(),
              vector(), vector(), vector(), vector(), vector(), vector())
for (i in 1:12) {
  start_month <- i
  if (start_month <= 3) {
    j <- 138
    while (j >= 0) {
      month[[i]][139-j] <- start_month + (138-j)*12
      j <- j - 1
    }
  }
  else {
    j <- 137
    while (j >= 0) {
      month[[i]][138-j] <- start_month + (137-j)*12
      j <- j - 1
    }
  }
}

# Get a slice of data on selected month
month_slice <- air_temp[,,month[[1]][50]] ## Number on the left side can be changed (1 is for January, 2 is for February...)
                                          ## Number on the right side is the year (1 means 1880, 2 means 1881...)

# Load necessary library
library(maps)

# Set to Plots folder to create plots
setwd("D:/University of Wollongong/Data Mining project/Plots")

# Create list of month names
month_name <- c("January", "February", "March", "April", "May", "June",
                "July", "August", "September", "October", "November", "December")

# Dynamically create folders
for(i in 1880:2018) {
  dir.create(as.character(i))
}

# Dynamically produce plots
for (i in 1:12) {
  if (i <= 3) {
    for (j in 1880:2018) {
      png(paste(j, "/", i, ". ", month_name[i], "-", j, " Anomaly.png", sep = ""), width=1500, height=750)
      month_slice <- air_temp[,,month[[i]][(j+1)-1880]]
      map_matrix <- pmax(pmin(month_slice,6),-6)
      par(mar = c(4,5,3,0))
      int <- seq(-6,6,length.out = 81)
      rgb.palette <- colorRampPalette(c("black","blue","darkgreen","green",
                                        "yellow","pink","red","maroon"),
                                      interpolate="spline")
      map_matrix <- map_matrix[,seq(length(map_matrix[1,]),1)]
      filled.contour(lon, lat, map_matrix, color.palette=rgb.palette, levels=int,
                     plot.title=title(main=paste("NOAAGlobalTemp Anomalies of ", month_name[i], "-", j, " [deg C]"),
                                      xlab="Latitude",ylab="Longitude", cex.lab=1.5),
                     plot.axes={axis(1, cex.axis=1.5);
                       axis(2, cex.axis=1.5);map("world2", add=TRUE);grid()},
                     key.title=title(main="[oC]"),
                     key.axes={axis(4, cex.axis=1.5)})
      dev.off()
    }
  }
  else {
    for (j in 1880:2017) {
      png(paste(j, "/", i, ". ", month_name[i], "-", j, " Anomaly.png", sep = ""), width=1500, height=750)
      month_slice <- air_temp[,,month[[i]][(j+1)-1880]]
      map_matrix <- pmax(pmin(month_slice,6),-6)
      par(mar = c(4,5,3,0))
      int <- seq(-6,6,length.out = 81)
      rgb.palette <- colorRampPalette(c("black","blue","darkgreen","green",
                                        "yellow","pink","red","maroon"),
                                      interpolate="spline")
      map_matrix <- map_matrix[,seq(length(map_matrix[1,]),1)]
      filled.contour(lon, lat, map_matrix, color.palette=rgb.palette, levels=int,
                     plot.title=title(main=paste("NOAAGlobalTemp Anomalies of ", month_name[i], "-", j, " [deg C]", sep = ""),
                                      xlab="Latitude",ylab="Longitude", cex.lab=1.5),
                     plot.axes={axis(1, cex.axis=1.5);
                       axis(2, cex.axis=1.5);map("world2", add=TRUE);grid()},
                     key.title=title(main="[oC]"),
                     key.axes={axis(4, cex.axis=1.5)})
      dev.off()
    }
  }
}

# Return back to original folder
setwd("D:/University of Wollongong/Data Mining project")

# Calculate the percentage of missing values for each year
temperature_vec <- as.vector(air_temp)
temperature_matrix <- matrix(temperature_vec, nrow = no_lat * no_lon, ncol = no_time)
dataFrame <- data.frame(cbind(coordinates, temperature_matrix))

# Change the name of columns
column_names <- vector()
column_names[[1]] <- "Latitude"
column_names[[2]] <- "Longitude"
for (i in 1:1659) {
  column_names[[i+2]] <- substr(format(time[i], format = c(dates = "dd/mm/yyyy")),4,12)
}
names(dataFrame) <- column_names

temp <- dataFrame
area_weight <- matrix(0, nrow = 2592, ncol = 1661)
area_weight[,1] <- temp[,1]
area_weight[,2] <- temp[,2]
veca <- cos(temp[,1] * pi / 180)
for(j in 3:1661) {
  for (i in 1:2592) {
    if(!is.na(temp[i,j])) area_weight[i,j]=veca[i]
  }
}
temp[is.na(temp)] <- 0
spatial_average <- area_weight * temp
spatial_average[,1:2] <- temp[,1:2]
average_vector <- colSums(spatial_average[,3:1661]) / colSums(area_weight[,3:1661])
time_month <- seq(1880,2018,length=1659)
plot(time_month, average_vector, type="l",
     cex.lab=1.4, xlab="Year", ylab="Temperature anomaly [deg C]",
     main="Area-weighted global average of monthly SST anomalies: Jan 1880 - Jan 2018")
abline(h = 0, col = "red", lty = 2)
abline(model <- lm(average_vector ~ time_month), col="blue", lwd = 2)
text(1930, 0.7, paste("Linear trend:", round(100*coef(model)[2],3), "[deg C] per century"),
     cex=1.4, col="blue")
summary(model)

# Calculate the percentage covered
rcover <- 100*colSums(area_weight[,3:1661])/sum(veca)
month_time <- seq(1880, 2018, length=1659)
plot(month_time, rcover, type="l", ylim=c(0,100),
     main="NOAAGlobalTemp Data Coverage: Jan 1880-Mar 2018",
     xlab="Year",ylab="Percent area covered [%]")

average_month <- matrix(average_vector[445:1656], ncol=12, byrow=TRUE)
year_vector <- seq(0,length = 101)
for(i in 1:101) year_vector[i] <- mean(average_month[i,])
average_monthly <- cbind(average_month, year_vector)
year_series <- seq(1917, 2017, len = 101)
png("Monthly_Anomaly.png", 800, 600)
par(mfrow = c(4, 3))
par(mgp = c(2,1,0))
for (i in 1:12) {
  plot(year_series, average_monthly[,i], type="l", ylim=c(-1.0,1.0),
       xlab="Year",ylab="Temp [deg C]",
       main = month_name[[i]])
  abline(lm(average_monthly[,i] ~ year_series), col="red")
  text(1945,0.7, paste("Trend oC/century=",
                       round(digits=3,
                             (100*coefficients(lm(average_monthly[,i] ~ year_series))[2]))),
       col="red")
}
dev.off()

# Fit a polynomial regression line into the dataset
poly9deg<-lm(average_vector ~ poly(time_month, 9, raw= FALSE))
summary(poly9deg)
poly25deg<-lm(average_vector ~ poly(time_month, 25, raw= FALSE))
summary(poly25deg)
plot(time_month, average_vector, type="s",
     cex.lab=1.4, lwd=2,
     xlab="Year", ylab="Temperature anomaly [deg C]",
     main="Annual SAT time series and its orthogonal polynomial 9 degree fits: 1880-2018")
lines(time_month,predict(poly9deg),col="blue", lwd=3)
plot(time_month, average_vector, type="s",
     cex.lab=1.4, lwd=2,
     xlab="Year", ylab="Temperature anomaly [deg C]",
     main="Annual SAT time series and its orthogonal polynomial 25 degree fits: 1880-2018")
lines(time_month,predict(poly25deg),col="blue", lwd=3)

# Test the residual standard error of the model
summary(model)$r.squared
summary(poly9deg)$r.squared
summary(poly20deg)$r.squared
plot(model$fitted.values, model$residuals, main = "Fitted vs. Residuals Linear Regression", xlab = "Fitted Values", ylab = "Residuals", type = "p")
axis(2, at = seq(-0.5,0.5,length=11))
plot(poly9deg$fitted.values, poly9deg$residuals, main = "Fitted vs. Residuals Polynomial 9 Degree", xlab = "Fitted Values", ylab = "Residuals", type = "p")
axis(2, at = seq(-0.5,0.5,length=11))
abline(h = c(-0.3, 0.3), col = "blue", lwd = 2, lty = 2)
plot(poly25deg$fitted.values, poly25deg$residuals, main = "Fitted vs. Residuals Polynomial 25 Degree", xlab = "Fitted Values", ylab = "Residuals", type = "p")
axis(2, at = seq(-0.5,0.5,length=11))
abline(h = c(-0.3, 0.3), col = "blue", lwd = 2, lty = 2)
