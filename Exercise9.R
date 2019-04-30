#specify the packages of interest
options(warn=-1)
options(repr.plot.width=6, repr.plot.height=4)

packages = c("tidyverse", "JoSAE", "nlme","mapdata", "maps", "maptools", "rgdal", "ggrepel", "Matrix")

## Check to see if package is available and load else install the package and its dependencies
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Variables of interest:
#  SegmentsInCounty: number of of segments of county.
#  SegementID: identificator for segment.
#  HACorn: hectares of corn for given segment.
#  HASoybeans: hectares of soybeans for given segment.
#  PixelsCorn: pixels for corn for given segment.
#  PixelsSoybeans: pixels for soybeans for given segment.
#  MeanPixelsCorn: mean of pixels for corn over all segments in given county.
#  MeanPixelsSoybeans: mean of pixels for soybeans over all segments in given county.
#  CountyName: county identificator of the segment.


######
### a)

# Load and prepare data
data(landsat)
crops_data <- landsat[, -9]

# Create a groupedData object
corn_data_grouped <- groupedData(HACorn ~ PixelsCorn | CountyName,
                                 data = crops_data)
soy_data_grouped <- groupedData(HASoybeans ~ PixelsSoybeans | CountyName,
                                data = crops_data)

# Make individual lm fits according to the grouped data
lm_corn_grouped <- lmList(corn_data_grouped)
lm_corn_grouped
lm_soy_grouped <- lmList(soy_data_grouped)
lm_soy_grouped

######
### b)

# Fit a linear mixed model for both crops such that segments share the same
# countywide random effect
lme_corn <- lme(HACorn ~ PixelsCorn, data = corn_data_grouped, random = ~ 1)
lme_corn
beta_corn <- lme_corn$coefficients$fixed

lme_soy <- lme(HASoybeans ~ PixelsSoybeans, data = soy_data_grouped, random = ~ 1)
lme_soy
beta_soy <- lme_soy$coefficients$fixed


######
### c)

# Population mean of the explanatory variables
pop_mean_corn <- unique(crops_data$MeanPixelsCorn)
pop_mean_soy <- unique(crops_data$MeanPixelsSoybeans)

# Mean over the observed segments
seg_mean <- aggregate(crops_data[3:6], by = list(crops_data$CountyName), mean)

# Number of observations in each county
num_obs_df <- plyr::count(crops_data, "CountyName")
num_obs <- num_obs_df[, 2]
county_names <- num_obs_df[, 1]

## Model for corn and soy
var_est_corn <- VarCorr(lme_corn)
var_est_soy <- VarCorr(lme_soy)
sigma_eps_corn <- as.numeric(var_est_corn[2])
sigma_eps_soy <- as.numeric(var_est_soy[2])
sigma_rand_corn <- as.numeric(var_est_corn[1])
sigma_rand_soy <- as.numeric(var_est_soy[1])

# Calculate covariance matrix V_hat of beta_hat
M_corn <- list()
M_soy <- list()
for (i in 1:12){
  M_corn[[i]] = matrix(sigma_rand_corn, num_obs[i], num_obs[i])
  M_soy[[i]] = matrix(sigma_rand_soy, num_obs[i], num_obs[i])
}

I_n_corn <- diag(x = sigma_eps_corn, sum(num_obs), sum(num_obs))
I_n_soy <- diag(x = sigma_eps_soy, sum(num_obs), sum(num_obs))
block_M_corn <- bdiag(M_corn)
block_M_soy <- bdiag(M_soy)
V_corn <- I_n_corn + block_M_corn
V_soy <- I_n_soy + block_M_soy

aux_corn = cbind(1, crops_data$PixelsCorn)
aux_soy = cbind(1, crops_data$PixelsSoybeans)
V_hat_corn <- solve(t(aux_corn) %*% solve(V_corn) %*% aux_corn)
V_hat_soy <- solve(t(aux_soy) %*% solve(V_corn) %*% aux_soy)

gamma_corn <- sigma_rand_corn/(sigma_rand_corn + sigma_eps_corn/num_obs)
gamma_soy <- sigma_rand_soy/(sigma_rand_soy + sigma_eps_soy/num_obs)

## Calculate the different predictors

# Regression predictor
reg_pred_corn <- cbind(1, pop_mean_corn) %*% beta_corn
reg_pred_soy <- cbind(1, pop_mean_soy) %*% beta_soy

# Adjusted survey predictor
adj_surv_pred_corn <- cbind(1, pop_mean_corn) %*% beta_corn +
  (seg_mean$HACorn - (cbind(1, seg_mean$PixelsCorn) %*% beta_corn))
adj_surv_pred_soy <- cbind(1, pop_mean_soy) %*% beta_soy +
  (seg_mean$HASoybeans - (cbind(1, seg_mean$PixelsSoybeans) %*% beta_soy))

# (Empirical) BLUP
BLUP_pred_corn <- cbind(1, pop_mean_corn) %*% beta_corn +
  gamma_corn*(seg_mean$HACorn - (cbind(1, seg_mean$PixelsCorn) %*% beta_corn))
BLUP_pred_soy <- cbind(1, pop_mean_soy) %*% beta_soy +
  gamma_soy*(seg_mean$HASoybeans - (cbind(1, seg_mean$PixelsSoybeans) %*% beta_soy))

# Survey predictor
surv_pred_corn <- seg_mean$HACorn
surv_pred_soy <- seg_mean$HASoybeans

# Save resulting predictors in a data frame
df_pred <- data.frame(County = county_names, 
                      Reg_pred_corn = reg_pred_corn,
                      Reg_pred_soy = reg_pred_soy,
                      Adj_surv_pred_corn = adj_surv_pred_corn,
                      Adj_surv_pred_soy = adj_surv_pred_soy,
                      BLUP_pred_corn = BLUP_pred_corn,
                      BLUP_pred_soy = BLUP_pred_soy,
                      Surv_pred_corn = surv_pred_corn, 
                      Surv_pred_soy = surv_pred_soy)
df_pred

# Estimate MSE
MSE_pred <- function(d, crop){
  # corn: whether we have the model for "corn" or "soy"
  # d: characteristic number for predictor
  # d = 0:      Regression predictor
  # d = 1:      Adjusted survey predictor
  # d = gamma:  (Empirical) BLUP
  # d = 2:      Survey predictor
  if (length(d) == 1){
    d = rep(d, 12)
  }
  if (crop == "corn"){
    sigma_eps = sigma_eps_corn
    sigma_rand = sigma_rand_corn
    gamma = gamma_corn
    pop_mean = pop_mean_corn
    seg_mean = seg_mean$PixelsCorn
    V_hat = V_hat_corn
  } else {
    sigma_eps = sigma_eps_soy
    sigma_rand = sigma_rand_soy
    gamma = gamma_soy
    pop_mean = pop_mean_soy
    seg_mean = seg_mean$PixelsSoybeans
    V_hat = V_hat_soy
  }
  
  res = rep(0, 12)
  for (i in 1:12){
    if (d[1] == 2){
      aux = (cbind(1, pop_mean[i]) - cbind(1, seg_mean[i]))
      
      res[i] = sigma_eps/num_obs[i] + aux %*% V_hat %*% t(aux)
    } else {
      aux1 = (cbind(1, pop_mean[i]) - d[i]*cbind(1, seg_mean[i]))
      aux2 = cbind(1, seg_mean[i])
      
      term1 = (1 - d[i])**2*sigma_rand + d[i]**2*sigma_eps/num_obs[i]
      term2 = 2*(d[i] - gamma[i])*aux1 %*% V_hat %*% t(aux2)
      term3 = aux1 %*% V_hat %*% t(aux1)
      
      res[i] = term1 + term2 + term3
    }
  }
  return(res)
}

MSE_reg_pred_corn <- MSE_pred(0, crop = "corn")
MSE_reg_pred_soy <- MSE_pred(0, crop = "soy")
MSE_adj_surv_pred_corn <- MSE_pred(1, crop = "corn")
MSE_adj_surv_pred_soy <- MSE_pred(1, crop = "soy")
MSE_BLUP_pred_corn <- MSE_pred(gamma_corn, crop = "corn")
MSE_BLUP_pred_soy <- MSE_pred(gamma_soy, crop = "soy")
MSE_surv_pred_corn <- MSE_pred(2, crop = "corn")
MSE_surv_pred_soy <- MSE_pred(2, crop = "soy")

df_MSE <- data.frame(County = county_names,
                     Reg_pred_corn = MSE_reg_pred_corn,
                     Reg_pred_soy = MSE_reg_pred_soy,
                     Adj_surv_pred_corn = MSE_adj_surv_pred_corn,
                     Adj_surv_pred_soy = MSE_adj_surv_pred_soy,
                     BLUP_pred_corn = MSE_BLUP_pred_corn,
                     BLUP_pred_soy = MSE_BLUP_pred_soy,
                     Surv_pred_corn = MSE_surv_pred_corn,
                     Surv_pred_soy = MSE_surv_pred_soy)
df_MSE


######
### d)

## Estimate the total county field size for both crops

county_num_segments <- unique(crops_data$MeanPixelsCorn)

# Estimates for total county field size
total_reg_pred_corn <- reg_pred_corn * county_num_segments
total_reg_pred_soy <- reg_pred_soy * county_num_segments
total_adj_surv_pred_corn <- adj_surv_pred_corn * county_num_segments
total_adj_surv_pred_soy <- adj_surv_pred_soy * county_num_segments
total_BLUP_pred_corn <- BLUP_pred_corn * county_num_segments
total_BLUP_pred_soy <- BLUP_pred_soy * county_num_segments
total_surv_pred_corn <- surv_pred_corn * county_num_segments
total_surv_pred_soy <- surv_pred_soy * county_num_segments

df_total_HA <- data.frame(county = county_names,
                          Reg_pred_corn = total_reg_pred_corn,
                          Reg_pred_soy = total_reg_pred_soy,
                          Adj_surv_pred_corn = total_adj_surv_pred_corn,
                          Adj_surv_pred_soy = total_adj_surv_pred_soy,
                          BLUP_pred_corn = total_BLUP_pred_corn,
                          BLUP_pred_soy = total_BLUP_pred_soy,
                          Surv_pred_corn = total_surv_pred_corn,
                          Surv_pred_soy = total_surv_pred_soy)
df_total_HA



## Plot the results by the BLUP from part (c) as well as the predictor
## only relying on the survey data in a table

df_total <- data.frame(County = county_names, 
                       BLUP_corn = total_BLUP_pred_corn,
                       BLUP_soy = total_BLUP_pred_soy,
                       Survey_corn = total_surv_pred_corn,
                       Survey_soy = total_surv_pred_soy)
df_total


## Plot the results onto a map of Iowa
states <- map_data("state")
iowa <- subset(states, region == "iowa")

counties <- map_data("county")
iowa_counties <- subset(counties, region == "iowa")

tl <- tolower(as.character(county_names))

iowa_counties_polygon <- select(iowa_counties, long, lat, subregion)
centroids <- aggregate(iowa_counties_polygon[,1:2], by=list(iowa_counties_polygon$subregion), FUN = mean)
centroids$County <- unique(iowa_counties_polygon[,3])
centroids <- filter(centroids, County %in% tl)

label_corn <- data.frame(total_BLUP = "BLUP:", BLUP = round(df_total$BLUP_corn, 0),
                         total_Survey = "Survey:", Survey = round(df_total$Survey_corn , 0))
label_soy <- data.frame(total_BLUP = "BLUP:", BLUP = round(df_total$BLUP_soy, 0),
                        total_Survey = "Survey:", Survey = round(df_total$Survey_soy, 0))
label_corn <- data.frame(BLUP = paste( "", label_corn$total_BLUP, "", label_corn$BLUP),
                         Survey = paste(label_corn$total_Survey, label_corn$Survey))
label_soy <- data.frame(BLUP = paste( "", label_soy$total_BLUP, "", label_soy$BLUP),
                        Survey = paste(label_soy$total_Survey, label_soy$Survey))
label_corn <- data.frame(Total_HA = paste("", label_corn$BLUP, "\n", label_corn$Survey))
label_soy <- data.frame(Total_HA = paste("", label_soy$BLUP, "\n", label_soy$Survey))

iowa_counties$fill_value <- 0
iowa_counties$fill_value[iowa_counties$subregion %in% tl] <- 1

# Map for corn
ggplot(data = iowa_counties, aes(x = long, y = lat), fill = ) + 
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "palegreen", color = "black") +
  geom_polygon(data = iowa_counties, aes(x = long, y = lat, group = group, fill = factor(fill_value)), 
               color = "white", show.legend = FALSE) +
  geom_polygon(data = iowa, aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  geom_label_repel(data = centroids, aes(long, lat, label = label_corn$Total_HA),
                   size = 3.5, alpha = 0.7, point.padding = 1.5,
                   min.segment.length = 0, segment.size = 0.6) +
  scale_fill_manual(values = c("grey80", "turquoise")) +
  theme_void()

# Map for soybeans
ggplot(data = iowa_counties, aes(x = long, y = lat), fill = ) + 
  geom_polygon(aes(x = long, y = lat, group = group), 
               fill = "palegreen", color = "black") +
  geom_polygon(data = iowa_counties, aes(x = long, y = lat, group = group, fill = factor(fill_value)), 
               color = "white", show.legend = FALSE) +
  geom_polygon(data = iowa, aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  geom_label_repel(data = centroids, aes(long, lat, label = label_soy$Total_HA),
                   size = 3.5, alpha = 0.7, point.padding = 1.5,
                   min.segment.length = 0, segment.size = 0.6) +
  scale_fill_manual(values = c("grey80", "khaki")) +
  theme_void()




