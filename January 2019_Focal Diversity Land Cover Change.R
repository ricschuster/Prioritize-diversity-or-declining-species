## Paper 2 figs
library(raster)
library(readxl)
library(dplyr)
library(plyr)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(here)
library(leaflet)
library(png)
library(fields)
rm(list=ls())

setwd("C:/Users/scott/Documents/R/working/Scott_trends")
owd<-getwd()

memory.size()
memory.size(max=FALSE)

# set to 1 if first time this code is run (HF processing, )
first <- 0
write <- 1

load("raster_setup.RData")
here("Rup") %>% 
  list.files(full.names = TRUE) %>% 
  walk(source)

# natural earth political boundaries
ne_land <- read_sf("data/ne-land.gpkg") %>% st_geometry()
ne_country_lines <- read_sf("data/ne-country-lines.gpkg") %>% st_geometry()
ne_state_lines <- read_sf("data/ne-state-lines.gpkg") %>% st_geometry()

idx.r <- nb.stack[[1]]
idx.r[] <- 1:length(nb.stack[[1]][])
figs <- read_xlsx("diversity.land.xlsx")
figs2 <- mutate(figs,
	Fig5a = figs$y2000,
	Fig5b = figs$ssp1,
	Fig5c = figs$ssp2,
	Fig5d = figs$ssp3)

# 1=forest, 2=mosaic forest, 3=peri-urban, 4=urban, 5=grassland/bare, 6=cropland/mosaic cropland-grass  
idx.df <- data.frame(idx = idx.r[])
idx.df <- join(idx.df, figs2, by = "idx")

fig.r <- idx.r
fig.r[] <- idx.df$Fig5a
plot(fig.r)

fig.lst <- list()
for(ii in 4:7){
  fig.r[] <- idx.df[,ii]
  fig.lst[[ii - 3]] <- fig.r
}
fig.st <- stack(fig.lst)
names(fig.st) <- c("Fig5a", "Fig5b", "Fig5c", "Fig5d")

add_legend <- function(title, palette, bump = 0, low_high = FALSE, 
                       text_col = "black") {
  if (low_high) {
    labs <- list(at = (seq(0:5)-1)/5, labels = c("one","two","three","four","five","six"), line = -1,
                 cex.axis = 0.75, fg = NA, col.axis = text_col)
    
  } else {
    labs <- list(at = c(0, 1), labels = NA, line = 0, fg = NA)
  }
  fields::image.plot(zlim = c(0,1), legend.only = TRUE, col = palette(1:6),
                     legend.width = 1, horizontal = TRUE,
                     smallplot = c(0.05, 0.35, 0.05 + bump, 0.075 + bump),
                     axis.args = labs,
                     legend.args = list(text = proper(title), side = 1,
                                        col = text_col)) 
}

proper <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
}

stem_crop <- function(x) {
  stopifnot(inherits(x, "Raster"))
  
  # aggregate for faster processing
  x_agg <- raster::aggregate(x, fact = 3)
  
  # extent of non-NA
  x_agg <- stem_to_na(x_agg)
  #x_agg <- raster::trim(x_agg, values = NA)
  #x_ext <- raster::extent(x_agg)
  x_ext <- extent_na(x_agg)
  raster::crop(x, x_ext)
}

extent_na <- function(x) {
  pts <- raster::rasterToPoints(x)
  x_rng <- range(pts[, "x"])
  y_rng <- range(pts[, "y"])
  raster::extent(x_rng[1] - res(x)[1] / 2, x_rng[2] + res(x)[1] / 2,
                 y_rng[1] - res(x)[2] / 2, y_rng[2] + res(x)[2] / 2)
}

stem_to_na <- function(x, value = 0) {
  stopifnot(inherits(x, "Raster"))
  stopifnot(is.numeric(value), length(value) == 1)
  
  if (inherits(x, "RasterLayer")) {
    x[x[] == value] <- NA_real_
  } else {
    for (i in seq.int(raster::nlayers(x))) {
      x[[i]][x[[i]][] == value] <- NA_real_
    }
  }
  return(x)
}

# process for visualization
crs <- "+proj=moll +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
abd_plot <- fig.st %>% 
  stem_to_na() %>% 
  projectRaster(crs = crs, method = "ngb") %>% 
  #sqrt() %>% 
  stem_crop()

# plot
e <- extent(abd_plot)
text_col <- "black"

fig <- names(abd_plot)
#palette <- c("Greens", "Greys")
palette <- c("Greens", "Blues", "BuPu", "Greys")
legend_offsets <- c(0.01, 0.06, 0.11, 0.16)
# prepare vector layers
land <- st_transform(ne_land, crs = crs)
country <- st_transform(ne_country_lines, crs = crs)
state <- st_transform(ne_state_lines, crs = crs)

# plot fig
for (ii in seq_along(fig)) {
title <- paste(fig[ii])
  # print map
  # here("out_all_spp/fig", paste0(fig[ii], ".png")) %>% 
  png(filename = paste0(fig[ii], ".png"), width = 3000, height = 3000, res = 300)
  par(mar = c(0, 0, 0, 0), oma = c(0,0,0,0), bg = "white")
  plot(land, col = "gray87", border = NA, xlim = e[1:2], ylim = e[3:4])

##*********************
#pal<-brewer.pal(6,"Greens")
#pal <- colorRampPalette(pal)
pal<-colorFactor(palette=c('forestgreen','chartreuse','royalblue1','magenta2',
'goldenrod3','yellow'),domain=NULL) # specify individual colors
nms <- c("one","two","three","four","five","six") # try this with the actual names and a rotated legend, try 
#pal <- colorQuantile(pal, values(abd_plot[fig]), n = 8, #probs = seq(0, 1, length.out = n + 1),
#              na.color = "#808080", alpha = FALSE, reverse = FALSE)
plot(abd_plot[[gsub(" ", ".",fig[ii])]], add = TRUE, col = pal(1:6), legend = FALSE, 
     maxpixels = ncell(abd_plot))
add_legend("", pal, legend_offsets[3], low_high = TRUE,
           text_col = text_col) # check legend rotation

# boundaries
plot(state, col = "black", lwd = 0.5, lty = 1, add = TRUE)
plot(country, col = "black", lwd = 1, add = TRUE)

# title
# plot bounds
usr <- par("usr")
xwidth <- usr[2] - usr[1]
yheight <- usr[4] - usr[3]
# labels
text(x = usr[1] + 0.05 * xwidth, y = usr[3] + 0.21 * yheight,
     labels = "something", pos = 4, font = 1, cex = 1.5, col = text_col)

#rasterImage(logo,usr[1] + 0.01 * xwidth, usr[3] + 0.03 * yheight,
#            usr[1] + 0.38 * xwidth, usr[3] + 0.09 * yheight)
dev.off()
}


