#####################################
####   Script for plotting a map ####
#####################################

##### goals #####

### Restart R and clear the current workspace

.rs.restartR()
rm(list=ls())

### because tidyverse and ggplot is totally f***ed we have to read in the map data before doing anything else (&%$&U^*%$)
map2 <- ggplot2::map_data(map = "world", region = c("UK","Ireland"))


### install if necessary and then load the libraries you need

j <- c("rstudioapi","ggplot2","RColorBrewer","ggmap","maps","grid")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded


### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## finally, let's plot a map figure

### read in the data - this is the .txt file you produced in the PreparingData.R script. 
coords<-read.csv("Data/Coordinates.csv", header=TRUE)

coords$Site <- as.factor(coords$Site)
summary(coords)

# check the points plot ok with a ggplot

gtest <- ggplot(coords, aes(x=LongDec,y=LatDec,group=Site))+
  geom_point(aes(colour=Site))

gtest



# now we want to apply them to a map
register_google(key = "XXXX", account_type = "premium", day_limit = 100000)


map1 <- get_map(location = c(lon=-1.045,lat=53.945), maptype="terrain", zoom=14,color="color")


ggmap(map1)


# set things up for a scale bar
# get match attributes
bb <- attr(map1,"bb")


# figure out points to define scale bar
sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                   lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                   lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                   lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))

sbar$distance <- geosphere::distVincentyEllipsoid(c(sbar$lon.start,sbar$lat.start),
                                                  c(sbar$lon.end,sbar$lat.end))

ptspermm <- 2.83464567  # apparently, we "need this because geom_text uses mm, and themes use pts. Urgh."

# set the length of the scale bar - here 0.5 km
scalebar.length <- 0.5
sbar$lon.end <- sbar$lon.start +
  ((sbar$lon.end-sbar$lon.start)/sbar$distance)*scalebar.length*1000


g1 <- ggmap(map1, extent = "device", crop = T) +
  geom_point(data = coords, size=5,colour="black",fill="white",stroke=2,alpha=0.8,
             aes(x=LongDec,y=LatDec,shape=Site))+
  scale_shape_manual(values=c(0,1,2,5),
                     labels=c("1","2","3","RIS trap"))+
  theme(legend.justification=c(0,0), legend.position=c(0.07,0.77))+ 
  theme(legend.text = element_text(size=20))+
  theme(legend.title = element_text(size=20))+
  theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="black"))+
  geom_segment(data = sbar, size=1.5,
               aes(x = lon.start+0.005,
                   xend = lon.end+0.005,
                   y = lat.start+0.001,
                   yend = lat.end+0.001),
               arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                           ends = "both", type = "open")) +
  geom_text(data = sbar,
            aes(x = ((lon.start + lon.end)/2)+0.005,
                y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat) + 0.001,
                label = paste(format(scalebar.length),
                              'km')),
            hjust = 0.5,
            vjust = 0,
            size = 20/ptspermm)  +
  coord_map(projection = "mercator",
            xlim=c(bb$ll.lon, bb$ur.lon),
            ylim=c(bb$ll.lat, bb$ur.lat))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.key = element_blank(),
        legend.key.height=unit(2,"line"),
        legend.key.width=unit(2,"line"))


g1



ggsave("FigS1close.svg", plot = g1, device = "svg", path = "Plots", width = 30, height = 30, units = "cm")



# we also need a big map to show where this map is
# we read in the map at the top

g2 <- ggplot(map2, aes(x = long, y = lat, group = group))+
  geom_polygon(fill = "white", colour = "black")+
  geom_point(data = coords[1,], aes(x = LongDec, y = LatDec, group = Site), size = 2)+
  geom_text(data = coords[1,], aes(x = LongDec, y = 0.5+LatDec, group = Site, label = "York"), size = 6)+
  theme_classic()+
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"))

g2

ggsave("FigS1far.svg", plot = g2, device = "svg", path = "Plots", width = 12, height = 15, units = "cm")


# combine the plots

plot2 <- ggplotGrob(g2)

# Calculate position in plot1 coordinates
# Extract x and y values from plot1
xleft   = 0.70
xright  = 0.96
ybottom = 0.05
ytop    = 0.40 


l1 = ggplot_build(g1)
x1 = l1$layout$panel_params[[1]]$x.range[1]
x2 = l1$layout$panel_params[[1]]$x.range[2]
y1 = l1$layout$panel_params[[1]]$y.range[1]
y2 = l1$layout$panel_params[[1]]$y.range[2]
xdif = x2-x1
ydif = y2-y1
xmin  = x1 + (xleft*xdif)
xmax  = x1 + (xright*xdif)
ymin  = y1 + (ybottom*ydif)
ymax  = y1 + (ytop*ydif) 


# inset g2 into g1

g3 <- g1 + inset(plot2,
                 xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

g3

ggsave("FigS1.svg", plot = g3, device = "svg", path = "Plots", width = 30, height = 30, units = "cm")

