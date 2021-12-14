library(ggplot2)
library(ggmap)
library(grid)
library(maps)
library(tigris)
library(sf)

f=map_data('state')
try <- subset(f, region =  "massachusetts" )

xx=range(try$long)
yy=range(try$lat)
g  <- ggplot() + coord_cartesian(xlim=xx, ylim = yy) +
  geom_polygon( data=try, 
                aes(x=long, y=lat, group=group), fill='gray80') +
  theme(
    panel.border = element_rect(colour = "gray80", fill='NA', size=1),
        axis.text =element_blank(),
        axis.ticks=element_blank(),
        axis.title =element_blank(),
        panel.grid=element_blank(),
    panel.background = element_blank()) +
    geom_point(aes(y=41.551115,x=-70.619438),fill='red',
           pch=22,,stroke=1,size=7)
g


f2=map_data('county')
NewStates2 <- filter(f2,subregion ==  "barnstable" | subregion=="nantucket")
xx=range(NewStates2$long)
yy=range(NewStates2$lat)
g2  <- ggplot() + coord_cartesian(xlim=xx, ylim = yy) +
  geom_polygon( data=NewStates, 
                aes(x=long, y=lat, group=group), fill='gray80') +
  #coord_fixed(xlim=c(-70.618,-70.62),ylim=c(41.54,41.56),clip='on') +
  theme(
    panel.border = element_rect(colour = "gray80", fill='NA', size=1),
    axis.text =element_blank(),
    axis.ticks=element_blank(),
    axis.title =element_blank(),
    panel.grid=element_blank(),
    panel.background = element_blank()) +
  geom_point(aes(y=41.551115,x=-70.619438),fill='red',
             pch=21,,stroke=1,size=7)
g2

counties_MA <- counties(state="MA",cb=T,class="sf",)
barn=county_subdivisions(state="MA", county = "Barnstable", cb=T)
barn1=subset(barn, NAME=="Falmouth")



counties_MA <- st_crop(counties_MA,
                       xmin=min_lon,xmax=max_lon,
                       ymin=min_lat,ymax=max_lat)

ggplot() + 
  geom_sf(data=barn1,fill="gray90",lwd=1, colour='black')+ 
  geom_point(aes(y=41.551115,x=-70.619438),fill='red',
             pch=21,,stroke=1,size=7) +
  theme_void()


### here 
a=area_water(state="MA", county="Barnstable", class='sf')
water <- st_crop(a,
                 xmin=-70.6,xmax=-70.7,
                 ymin=41.52,ymax=41.57)
land2 <- st_crop(barn1,
                 xmin=-70.6,xmax=-70.7,
                 ymin=41.52,ymax=41.57)


ggplot() + 
  geom_sf(data=water,fill="gray90",lwd=1, colour='black')+ 
  geom_point(aes(y=41.551115,x=-70.619438),fill='red',
             pch=21,,stroke=1,size=7) +
  theme_void()

ggplot() + 
  geom_sf(data=land2, fill='forestgreen')+
  geom_sf(data=a,fill='blue',alpha=0.3,
          inherit.aes = F,
          col="black")+
  coord_sf(xlim = c(-70.6, -70.69), 
           ylim = c(41.51, 41.57))+
  theme(legend.position = F) + theme_void() +
  geom_point(aes(y=41.551115,x=-70.619438),fill='red',
             pch=21,,stroke=1,size=7)

st_erase <- function(x, y) {
  st_difference(x, st_union(y))
}
counties_MA <- st_erase(land2,a)


ggplot() + 
  geom_sf(data=land2, fill='forestgreen')+
  geom_sf(data=a,fill='blue',alpha=0.3,
          inherit.aes = F,
          col="black")+
  coord_sf(xlim = c(-70.6, -70.69), 
           ylim = c(41.51, 41.57))+
  theme(legend.position = F) + theme_void() +
  geom_point(aes(y=41.551115,x=-70.619438),fill='red',
             pch=21,,stroke=1,size=7)





ggplot() + 
  geom_sf(data=counties_MA,fill='gray70',lwd=0)+
  theme(legend.position = F) + theme_void() +
  geom_point(aes(y=41.551115,x=-70.619438),fill='red',
             pch=21,,stroke=1,size=7) +
  theme(panel.background=element_rect(fill = 'white'))



land3 <- st_crop(a,
                 xmin=-70.6187,xmax=-70.6258,
                 ymin=41.545,ymax=41.556)





### good plot 
barn=county_subdivisions(state="MA", county = "Barnstable", cb=T)
barn1=subset(barn, NAME=="Falmouth")
a=area_water(state="MA", county="Barnstable", class='sf')
water <- st_crop(a,
                 xmin=-70.6,xmax=-70.7,
                 ymin=41.52,ymax=41.57)
land2 <- st_crop(barn1,
                 xmin=-70.6,xmax=-70.7,
                 ymin=41.52,ymax=41.57)

st_erase <- function(x, y) {
  st_difference(x, st_union(y))
}
counties_MA <- st_erase(land2,a)


saveRDS(counties_MA, "MAP_Barnstable.RDS")


siderspond <- st_crop(a,
                 xmin=-70.6187,xmax=-70.6256,
                 ymin=41.545,ymax=41.556)
try2= ggplot() + 
  geom_sf(data=counties_MA,fill='forestgreen',lwd=0.5, alpha=0.5)+
  geom_sf(data=siderspond,fill='navy',lwd=0) +
  theme(legend.position = F) + theme_void() +
  theme(panel.background=element_rect(fill = 'white'))
try2

NewStates <- county_subdivisions(state="MA")
#NewStates_ag=area_water(state="MA", county = "Barnstable")
# look up coastline 
#NewStates_ag=coastline()
NewStates_ag=linear_water(state="MA", county=NA)
f=map_data('state')
try <- subset(f, region ==  "massachusetts" , county=NA)
dim(try)
xx=range(try$long)
yy=range(try$lat)



f=map_data('state')
try <- subset(f, region =  "massachusetts" )

xx=range(try$long)
yy=range(try$lat)

g  <- ggplot() + coord_cartesian(xlim=xx, ylim = yy) +
  geom_polygon( data=try, 
                aes(x=long, y=lat, group=group), fill='forestgreen', alpha=0.5, lwd=0.5, colour='black') +
  theme(
    panel.border = element_rect(colour = "gray80", fill='NA', size=1),
    axis.text =element_blank(),
    axis.ticks=element_blank(),
    axis.title =element_blank(),
    panel.grid=element_blank(),
    panel.background = element_blank()) + theme_void() +
  #geom_point(aes(y=41.551115,x=-70.619438),fill='red',
   #          pch=22,,stroke=1,size=7)
 # geom_rect(aes(xmin=-70.6,xmax=-70.7,
   #             ymin=41.52,ymax=41.57), fill='black')
   geom_rect(aes(xmin=-70.68778,xmax=-70.498,
               ymin=41.51478,ymax=41.65841), fill='gray50')        
g

library(ggtern)


grid.draw(g)
grid.draw(try2)
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) # plot area for the main map
v2<-viewport(width = 0.5, height = 0.6, x = 0.3, y = 0.27) # plot area for the inset map
dev.off()
print(g,vp=v1) 
print(try2,vp=v2)
dev.off()




