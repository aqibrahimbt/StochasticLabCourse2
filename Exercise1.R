#specify the packages of interest
options(warn=-1)
#options(repr.plot.width=10, repr.plot.height=8)

packages = c("haven", "tidyverse", "plyr", "ggplot2", "dplyr", "maptools", "rgdal", "raster")

## Check to see if package is available and load else install the package and its dependencies
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

## data processing
##Read in data
children <- read_dta("data/childrenfinal.dta")

#remove all variables that start with “s“, “v” and “m”, followed by a number
children1<- children %>% 
  dplyr::select(-matches("^[svm][0-9]"))

#convert labelled double into double variables
children2<- children1 %>%
  mutate_if(is.double, as.double)

children3<- children2 %>%
  dplyr::select(c(hypage, ruralfacto, female, zstunt, zweight, zwast, adm2))

#Make a scatter plot of zstunt against hypage, Add a smooth line to the plot
ggplot(children3, aes(x = hypage, y = zstunt)) +
  geom_point() +
  geom_smooth(se = F, color = "green") + theme_classic()

#smooth plots of zstunt against age for gender
ggplot(children3, aes(x = hypage, y = zstunt,  colour = factor(female))) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = F) +
  scale_colour_manual(labels = c("male", "female"), values = c("red", "blue")) +
  guides(colour = guide_legend(title="Gender")) + theme_classic()

#plot zstunt against age for urban and rural children
ggplot(children3, aes(x = hypage, y = zstunt,  colour = factor(ruralfacto))) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = F) +
  scale_colour_manual(labels = c("urban", "rural"), values = c("orange", "purple")) +
  guides(colour = guide_legend(title="Area")) + theme_classic()


Kenya1<-getData("GADM", country="KE", level=1) #Download Kenya shapefile data

Kenya1_UTM<-spTransform(Kenya1, CRS("+init=epsg:32537"))
colnames(children3)[7] <- "NAME_1" 

children3<- children3[order(children3$NAME_1),]
Kenya1_UTM@data<- Kenya1_UTM@data[order(Kenya1_UTM@data$NAME_1),]

#summarising children3 data by the mean of zstunt in the corresponding county. Not sure why summarize doesnt work
#3 without specifying dplyr. BUG???
children4 <- children3 %>%
  dplyr::group_by(NAME_1) %>%
  dplyr::summarise(mean = mean(zstunt), n = n())

#Adding the missing county Isiolo
children4[nrow(children4) + 1,] <- NA
children4$NAME_1[47] <- "Isiolo"
children4<- children4[order(children4$NAME_1),]


#Prepare the dataframe for ggplot
Kenya1_UTM@data$id <- rownames(Kenya1_UTM@data)
Kenya1_UTM@data <- mutate(Kenya1_UTM@data, zstunt.mean= children4$mean)
Kenya1_df <- fortify(Kenya1_UTM)
Kenya1_df <- full_join(Kenya1_df,Kenya1_UTM@data, by="id")


## Creating a dataframe for the names on the map
centroids_df <- as.data.frame(coordinates(Kenya1_UTM))
names(centroids_df) <- c("long", "lat")
children4<- children4[order(children4$NAME_1),]
centroids_df$NAME_1 <- Kenya1_UTM@data$NAME_1
centroids_df$zstunt.mean <- children4$mean

#Generating the map
ggplot(data = Kenya1_df, aes(x = long, y = lat, group = group, fill = zstunt.mean)) + 
  geom_polygon(color = "black", size = 0.25) +
  geom_text(data = centroids_df, aes(x = long, y = lat, label = NAME_1, group = NULL), size = 3, position=position_jitter(width=1,height=1)) +
  scale_fill_distiller(name="Mean zstunt by county", palette = "Set1") +
  theme(aspect.ratio = 1) + theme_classic()


