#####################################################################
# Exploring Spatiotemporal Trends in Cetacean Stranding in Scotland #
#####################################################################
# Set up ----
## Packages
library(dplyr)
library(mgcv) 
library(lmtest)
library(lme4)
library(rphylopic)
library(ggeffects)
library(ggplot2)
library(scales)
library(DHARMa)
library(gridExtra)
library(factoextra)
library(spatstat)
library(sp)
library(sf)
library(spdep)
library(cowplot)

## set working directory 

## get dataframe
df <- read.csv("TDDTT_df.csv")

## set correct data types
df$year <- as.numeric(df$year)
df$month <- as.numeric(df$month)

## make knots for continuation of month
knots <- list(month = c(0.5, seq(1, 12, length = 10), 12.5)) #So month is continuous and distance between Dec - Jan is the same

## GAMM Over dispersion function
overdispersion_gamm <- function(model, data) {
  # Calculate the residual deviance
  rdev <- sum(residuals(model$lme)^2)
  
  # Calculate the number of fixed effect degrees of freedom
  mdf <- length(fixef(model$lme))
  
  # Calculate the residual degrees of freedom
  rdf <- nrow(data) - mdf
  
  # Calculate and return the rdev/rdf ratio
  rdev_rdf_ratio <- rdev / rdf
  return(rdev_rdf_ratio)
}

## GAM Overdispersion function
overdispersion_gam <- function(gam_model) {
  # Calculate the residual deviance
  residual_deviance <- gam_model$null.deviance - gam_model$deviance
  
  # Calculate the residual degrees of freedom
  residual_degrees_of_freedom <- gam_model$df.residual
  
  # Calculate the residual scaled deviance to residual degrees of freedom ratio
  ratio <- residual_deviance / residual_degrees_of_freedom
  
  return(ratio)
}

#####################################################################
# Exploration ----
## DATASETS FOR EACH SPECIES ----

## get the unique levels of spec.r
unique_spec <- unique(df$spec.r)

## initialise empty list to store results of each group
spec_list <- list()

## loop through each species to generate a count with zeroes added in for month and year where relevant
for(spec in unique_spec) {
  # Filter the dataframe for the current species
  spec_df <- filter(df, df$spec.r == spec)
  
  # Add zero data
  t1 <- data.frame(year = c(1992:2022), month = rep(c(1:12), each = 31), n = c(0))
  
  # Count the occurrences for the current species
  spec_m <- count(spec_df, year, month)
  
  # Merge with zero data
  spec_m <- merge(t1, spec_m, all = TRUE)
  
  # Handle duplicates
  spec_m$duplicate <- duplicated(spec_m[c("year", "month")])
  spec_m$duplicate <- ifelse(duplicated(spec_m[,c("year", "month")], fromLast=T), yes=T, no = spec_m$duplicate)
  spec_result <- spec_m[!(spec_m$duplicate == TRUE & spec_m$n == 0),] 
  
  # Store the result in the list with the species name
  spec_list[[spec]] <- spec_result
}

## SUMMARY STATS ----
## Make list to store summary information for each species 
summary_spec <- list()

## Run loop to generate summary information for each species 
for(spec in names(spec_list)) {
  # Extract df for each species
  spec_df <- spec_list[[spec]]
  
  # Run summary through the spec_list dataframes
  summary_result <- summary(spec_df)
  
  # Store the result in the list with the species name
  summary_spec[[spec]] <- summary_result
  
}

## VISUALISE SPECIES TRENDS TO LOOK FOR GROUPINGS ----
#loess fit used here as flexible and will illustrate the type of model that would be most appropriate to
#use for the grouping, as well as illustrating similar trends across species which will aid in grouping validation. 

## Annual trends ##
# Make a list to store the figures in annual
annual_spec <- list()

# Loop through each dataframe in the group_list
for(spec in names(spec_list)) {
  
  # Extract the current dataframe
  spec_df <- spec_list[[spec]]
  
  # Annual visualisation
  spec_p <- ggplot(spec_df, aes(year, n)) +
    geom_smooth(method = "loess") +
    ggtitle(spec) +
    scale_x_continuous(breaks=pretty_breaks()) +
    scale_y_continuous(breaks=pretty_breaks()) +
    theme_classic()
  
  # Produce list of plots
  annual_spec[[spec]] <- spec_p
}

# Arrange the plots into a grid with 4 columns
do.call(grid.arrange, c(annual_spec, ncol=4))

## Seasonal trends ##
# Make a list to store the figures in annual
seasonal_spec <- list()

# Loop through each dataframe in the group_list
for(spec in names(spec_list)) {
  
  # Extract the current dataframe
  spec_df <- spec_list[[spec]]
  
  # Seasonal visualisation
  spec_p <- ggplot(spec_df, aes(month, n)) +
    geom_smooth(method = "loess") +
    ggtitle(spec) +
    scale_x_continuous(breaks=pretty_breaks()) +
    scale_y_continuous(breaks=pretty_breaks()) +
    theme_classic()
  
  # Produce list of plots
  seasonal_spec[[spec]] <- spec_p
}

# Arrange the plots into a grid with 4 columns

tiff('species_seasonal.tiff', units="in", width=15, height=17, res=1000)

do.call(grid.arrange, c(seasonal_spec, ncol=4))

dev.off()

## common dolphins distinct to other pelagics - group seperately. 

## DATASETS FOR EACH GROUP ----
# Generate unique count datasets for each group

## get the unique levels of group.r
unique_groups <- unique(df$group.r)

## initialise empty list to store results of each group
group_list <- list()

## loop through each group to generate a count with zeroes added in for month and year where relevant
for(group in unique_groups) {
  # Filter the dataframe for the current group
  group_df <- filter(df, df$group.r == group)
  
  # Add zero data
  t1 <- data.frame(year = c(1992:2022), month = rep(c(1:12), each = 31), n = c(0))
  
  # Count the occurrences for the current group
  group_m <- count(group_df, year, month)
  
  # Merge with zero data
  group_m <- merge(t1, group_m, all = TRUE)
  
  # Handle duplicates
  group_m$duplicate <- duplicated(group_m[c("year", "month")])
  group_m$duplicate <- ifelse(duplicated(group_m[,c("year", "month")], fromLast=T), yes=T, no = group_m$duplicate)
  group_result <- group_m[!(group_m$duplicate == TRUE & group_m$n == 0),] 
  
  # Store the result in the list with the group name
  group_list[[group]] <- group_result
}

## COLINEARITY TEST ##  ----

## Make list to store summary information for each species 
ct_group <- list()

## Run loop to generate colinearity for each species time series data 
for(group in names(group_list)) {
  # Extract df for each species
  group_df <- group_list[[group]]
  
  # Run model on all species 
  model_group <- lm(n ~ year, data = group_df)
  
  # Run Durbin-Watson colinearity test
  ct_test_group <- dwtest(model_group)
  
  # Print results in list 
  ct_group[[group]] <- ct_test_group
  
}
## Look at ct_group lists for collinearity
## All significant, colinearity in timeseries data as expected - include ARMA in models 

## DISTRIBUTION TEST (ZERO INFLATION) ----

# Make a list to store the histograms 
hist_group <- list()

# Loop through each dataframe in the group_list
for(group in names(group_list)) {
  
  # Extract the current dataframe
  group_df <- group_list[[group]]
  
  # Create histogram using ggplot
  group_hist <- ggplot(group_df, aes(x = n)) +
    geom_histogram(binwidth = 1, fill = "lightblue", color = "black") + 
    ggtitle(group) +
    theme_classic()
  
  # Produce list of plots
  hist_group[[group]] <- group_hist
}

# Arrange the plots into a grid with 4 columns
do.call(grid.arrange, c(hist_group, ncol=2))
## common dolphins, baleen whales, deep divers all have lots of zeroes. Check ZIp distribution options in models 


#####################################################################
# Temporal Trends ----
## Model Comparisons ----
# Run models on each group

# Initialize a list to store the AIC values for each model
aic_results <- list()

# Loop through each dataframe in the group_list
for(group in names(group_list)) {
  
  # Extract the current dataframe
  group_df <- group_list[[group]]
  
  # Fit each model
  
  ## GLM
  model_glm <- glm(n ~ month + year, data = group_df, family = poisson(link = "log"))
  
  ## GAM with no autocorrelation
  model_gam_no_ac <- gam(n ~ s(month, bs = "cc", k = 12) + s(year), 
                         knots = knots, method = "REML", data = group_df, 
                         family = poisson(link = "log"))
  
  ## GAMM with autocorrelation (AC 1)
  model_gamm_ac1 <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), 
                         knots = knots, method = "REML", data = group_df, 
                         family = poisson(link = "log"), 
                         correlation = corARMA(form = ~ 1|year, p = 1))
  
  ## GAMM with random effect of year
  model_gamm_re <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year, bs = "re"), 
                        knots = knots, method = "REML", data = group_df, 
                        family = poisson(link = "log"), 
                        correlation = corARMA(form = ~ time, p = 1))
  
  ## GAMM without year
  model_gamm_no_year <- gamm(n ~ s(month, bs = "cc", k = 12), 
                             knots = knots, method = "REML", data = group_df, 
                             family = poisson(link = "log"), 
                             correlation = corARMA(form = ~ time, p = 1))
  
  ## GAMM with interaction of year
  model_gamm_no_year <- gamm(n ~ s(month, bs = "cc", k = 12, by = as.factor(year)) + as.factor(year), 
                             knots = knots, method = "REML", data = group_df, 
                             family = poisson(link = "log"), 
                             correlation = corARMA(form = ~ time, p = 1))
  
  ## GAMM with negative binomial (NB) with autocorrelation
  model_gamm_nb_ac <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), 
                           knots = knots, method = "REML", data = group_df, 
                           family = negative.binomial(1), 
                           correlation = corARMA(form = ~ time, p = 1))
  
  ## Model with ZIP
  model_gam_ziP_no_ac <- gam(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                        data = group_df, family = ziP())
  
  
  # Store the AIC values for each model in a list
  aic_values <- AIC(model_glm, model_gam_no_ac, 
                    model_gamm_ac1$lme, model_gamm_re$lme, 
                    model_gamm_no_year$lme, model_gamm_nb_ac$lme, 
                    model_gamm_nb_no_ac$lme, model_gam_ziP_no_ac)
  
  # Add the AIC values to the results list with the group name
  aic_results[[group]] <- aic_values
}

# Convert the results list to a data frame for easy viewing
aic_table <- do.call(rbind, aic_results)

# AIC RESULTS
#Baleen whales best model = GAM with no AC, but colinearity in the dataset so test this *
#Common dolphin best model = GAM with no AC, but colinearity in the dataset so test this *
#Harbour porpoise best model = GAMM with AC1 & negative binomial distribution, no interaction terms or random effects *
#Pelagic dolphin best model = GAMM with AC1 & Poisson distribution, no interaction terms or random effects
#Deep diver best model = GAM with ziP distribution 

## BALEEN WHALES ----

## extract df from results list
baleen_df <- group_list[["Baleen Whales"]]

  ## run best fit model according to AIC table
  b_model <- gam(n ~ s(month, bs = "cc", k = 12) + s(year), 
                knots = knots, method = "REML", data = baleen_df, 
                family = poisson(link = "log"))

  ## run model with AIC and compare fit
  b_model2 <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), 
                knots = knots, method = "REML", data = baleen_df, 
                family = poisson(link = "log"), 
                correlation = corARMA(form = ~ 1|year, p = 1))
#compare
gam.check(b_model)
gam.check(b_model2$gam) # better fit

## conduct diagnostics 
summary(b_model2$lme) # phi = -0.06
summary(b_model2$gam) # r2 = .37 (its an estimate for a gamm)
overdispersion_gamm(b_model, baleen_df) # 1.16, not amazing but ok as not above 1.5

## test for residual autocorrelation
## change layout so you can see both plots at once
layout(matrix(1:2, ncol = 2))

## extract the normalised residuals
res <- resid(b_model2$lme, type = "normalized")

## run acf and partial acf plots 
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors") # looks like it is captured well

## return layout to normal
layout(1) 

## look at results
plot(b_model2$gam, resid = TRUE) ##looks good
summary(b_model2$gam) #both significant

## visualise
## make a df count for year
b_year <- count(df, group.r, year)
b_year <- filter(b_year, group.r == "Baleen Whales")

## make a simplistic year model for plotting
byear <- gamm(n ~ s(year), data = b_year, family = poisson(link = "log"), method = "REML", 
              correlation = corARMA(form = ~ year, p = 1))

## get minke phlopic for plot
minke_uuid <- get_uuid(name = "Balaenoptera acutorostrata")
minke <- get_phylopic(uuid = minke_uuid, format = "vector")

## make minke stranded
dead_minke <- rotate_phylopic(minke, angle = 180)
dead_minke <- flip_phylopic(minke, vertical = TRUE, horizontal = FALSE)

## generate predictions from model to plot
b_month_predict <- ggpredict(b_model2$gam, terms = "month", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 
b_year_predict <- ggpredict(byear$gam, terms = "year", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 

## year plot
b_year_predict$year <- c(1992:2022)
b_year_predict$year <- as.numeric(b_year_predict$year) 
b_year$year <- c(1992:2022)

bal_year <- ggplot (b_year_predict, aes(year, predicted)) +
  #geom_point(mapping = aes(year, n), data = b_year, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab (" ") +
  xlab (" ") + 
  ggtitle("Baleen Whales") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_minke, x = 1995, y = 40, ysize = 12, alpha = 1)

## seasonal plot
bal_season <- ggplot (b_month_predict, aes(x, predicted, colour = NULL)) +
  #geom_point(mapping = aes(month, n), data = baleen_df, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab ("") +
  xlab ("") + 
  ggtitle("Baleen Whales") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_minke, x = 2, y = 2, ysize = 0.5, alpha = 1)

## SENSITIVITY CHECKS ##

# PER SPECIES
# species plots on line 179 plots shows trends are similar in baleen whales, but check for differences 

## Make dataset
baleen.group <- count(df, spec.r, group.r, year, month) %>%
  filter(group.r == "Baleen Whales")

baleen.group$spec.r <- as.factor(baleen.group$spec.r)
baleen.group <- filter(baleen.group, spec.r != "Cetacean (indeterminate species)")

## Run basic model to observe differences between species
group.mod <- gamm(n ~ s(month, bs = "cc", k = 12, by = spec.r) + spec.r*month + s(year), knots = knots, method = "REML", 
                  data = baleen.group, family = poisson(link="log"), 
                  correlation = corARMA(form = ~ 1|year, p = 1))

## Look at diagnostics 
gam.check(group.mod$gam)
summary(group.mod$lme) #0.42
overdispersion_gamm(group.mod, baleen.group)

## Results
summary(group.mod$gam)
plot(group.mod$gam) 
## No difference between species

## Minkes make up most of baleens, check trend without minke. 
## Make df 
baleen.df <- filter(df, df$group.r == "Baleen Whales")
no.mink <- filter(baleen.df, baleen.df$spec.r != "Minke whale")
no.mink <- count(no.mink, year, month)

## run model 
nominke.mod <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                    data = no.mink, family = poisson(link="log"), 
                    correlation = corARMA(form = ~ 1|year, p = 1))

## diagnostics
gam.check(nominke.mod$gam)
summary(nominke.mod$lme) #0.07

## results
summary(nominke.mod$gam) 
plot(nominke.mod$gam) 
## Trends hold, but no longer significant without minke's. Likely driven by lack of fin & humpback datapoints. 
# Can confirm group is robust. 

## TIME INFLUENCE
# Add a pseudo year with strandings to see if trend holds

## Make df
random <- data.frame(year=2023, month= c(1:12), n=rep(4, each = 12))
yeartest.baleen <- merge(random, baleen.group, all = TRUE)

## Run model 
year.mod.baleen <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                 data = yeartest.baleen, family = poisson(link="log"), 
                 correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.mod.baleen$gam)
summary(year.mod.baleen$lme) #phi = -0.038

## Results
plot(year.mod.baleen$gam, residuals = TRUE) 
summary(year.mod.baleen$gam)
## Trend holds. 


#Remove a year and see if trend holds

## Remove year
baleen.less <- filter(baleen.group, year != "2018")

## Run model
year.less.baleen <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                  data = baleen.less, family = poisson(link="log"), 
                  correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.less.baleen$gam)
summary(year.less.baleen$lme) #phi = -0.5
overdispersion_gamm(year.less.baleen, baleen.less) #1.19

## Results 
summary(year.less.baleen$gam)
plot(year.less.baleen$gam, residuals = TRUE) 
# Trend holds. 

## COMMON DOLPHINS ----
## extract df from results list
common_df <- group_list[["Common Dolphin"]]

## RUN MODELS
  ## run best fit model according to AIC table + includes AC structure
  c_model <- gam(n ~ s(month, bs = "cc", k = 12) + s(year), 
                knots = knots, method = "REML", data = common_df, 
                family = negative.binomial(1))

  ## run with AC
  c_model2 <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), 
                  knots = knots, method = "REML", data = common_df, 
                  family = negative.binomial(1), 
                  correlation = corARMA(form = ~ 1|year, p = 1))
#compare
gam.check(c_model)
gam.check(c_model2$gam) #fit better

## conduct diagnostics 
summary(c_model2$lme) # phi = 0.07
summary(c_model2$gam) # r2 = .52 (its an estimate for a gamm)
overdispersion_gamm(c_model2, common_df) # 1.53 but a nb

## test for residual autocorrelation
## change layout so you can see both plots at once
layout(matrix(1:2, ncol = 2))

## extract the normalised residuals
res <- resid(c_model2$lme, type = "normalized")

## run acf and partial acf plots 
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors") # looks like it is captured well

## return layout to normal
layout(1) 

## look at results
plot(c_model2$gam, resid = TRUE) ##looks good
summary(c_model2$gam) #both significant

## visualise
## make a df count for year
c_year <- count(df, group.r, year)
c_year <- filter(c_year, group.r == "Common Dolphin")

## make a simplistic year model for plotting
cyear <- gamm(n ~ s(year), data = c_year, family = poisson(link = "log"), method = "REML", 
              correlation = corARMA(form = ~ year, p = 1))

## get common dolphin phylopic for plot
cd_uuid <- get_uuid(name = "Delphinus delphis")
cd <- get_phylopic(uuid = cd_uuid, format = "vector")

## make common dolphin stranded
dead_cd <- rotate_phylopic(cd, angle = 180)

## generate predictions from model to plot
c_month_predict <- ggpredict(c_model2$gam, terms = "month", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 
c_year_predict <- ggpredict(cyear$gam, terms = "year", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 

## year plot
c_year_predict$year <- c(1992:2022)
c_year_predict$year <- as.numeric(c_year_predict$year) 
c_year$year <- c(1992:2022)

cd_years <- ggplot (c_year_predict, aes(year, predicted)) +
  #geom_point(mapping = aes(year, n), data = cd_year, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab (" ") +
  xlab (" ") + 
  ggtitle("Common Dolphins") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_cd, x = 1994, y = 75, ysize = 19, alpha = 1)

## seasonal plot
cd_season <- ggplot (c_month_predict, aes(x, predicted, colour = NULL)) +
  #geom_point(mapping = aes(month, n), data = common, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab ("") +
  xlab ("") + 
  ggtitle("Common Dolphins") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_cd, x = 2, y = 1.75, ysize = 0.4,alpha = 1)



## SENSITIVTY CHECKS ##

## Make dataframe
common.group <- count(df, spec.r, group.r, year, month) %>%
  filter(group.r == "Common Dolphin")

common.group$spec.r <- as.factor(common.group$spec.r)
common.group <- filter(common.group, spec.r != "Cetacean (indeterminate species)")


## TIME INFLUENCE
# Add a pseudo year with strandings to see if trend holds

## Make df
random <- data.frame(year=2023, month= c(1:12), n=rep(4, each = 12))
yeartest.common <- merge(random, common.group, all = TRUE)

## Run model 
year.mod.common <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                        data = yeartest.common, family = poisson(link="log"), 
                        correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.mod.common$gam)
summary(year.mod.common$lme) #phi = -0.038

## Results
plot(year.mod.common$gam, residuals = TRUE) 
summary(year.mod.common$gam)
## Trend holds. 

# Remove a year and see if trend holds
## Remove year
common.less <- filter(common.group, year != "2018")

## Run model
year.less.common <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                         data = common.less, family = poisson(link="log"), 
                         correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.less.common$gam)
summary(year.less.common$lme) #phi = -0.5
overdispersion_gamm(year.less.common, common.less) #1.19

## Results 
summary(year.less.common$gam)
plot(year.less.common$gam, residuals = TRUE) 
# Trend holds. 

## HARBOUR PORPOISE ----
## get best model, run diagnostics and visuals

## extract df from results list
harbour_df <- group_list[["Harbour porpoise"]]

## RUN MODELS ##
  ## run best fit model according to AIC table + includes AC structure
  hp_model <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML",
                 data = harbour_df, family = negative.binomial(1),
                 correlation = corARMA(form = ~ 1|year, p = 1))

  ## Test with poisson, distribution looked ok in exploration
  hp_model1 <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML",
                 data = harbour_df, family = poisson(link = "log"),
                 correlation = corARMA(form = ~ 1|year, p = 1))
  
# Look at rdev/rdf ratio
overdispersion_gamm(hp_model1, harbour_df) # 0.26
  
## conduct diagnostics 
gam.check(hp_model1$gam) #ok
summary(hp_model1$lme) # phi = 0.21
summary(hp_model1$gam) # r2 = .40 (its an estimate for a gamm)

## test for residual autocorrelation
## change layout so you can see both plots at once
layout(matrix(1:2, ncol = 2))

## extract the normalised residuals
res <- resid(hp_model1$lme, type = "normalized")

## run acf and partial acf plots 
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors") # looks like it is captured well

## return layout to normal
layout(1) 

## look at results
plot(hp_model1$gam, resid = TRUE) ##looks good
summary(hp_model1$gam) #both significant

## visualise
## make a df count for year
h_year <- count(df, group.r, year)
h_year <- filter(h_year, group.r == "Harbour porpoise")

## make a simplistic year model for plotting
hyear <- gamm(n ~ s(year), data = h_year, family = poisson(link = "log"), method = "REML", 
              correlation = corARMA(form = ~ year, p = 1))

## get porpoise phylopic for plot
hp_uuid <- get_uuid(name = "Phocoena phocoena")
porpoise <- get_phylopic(uuid = hp_uuid, format = "vector")

## make porpoise stranded
dead_hp <- rotate_phylopic(porpoise, angle = 180)

## generate predictions from model to plot
hp_month_predict <- ggpredict(hp_model1$gam, terms = "month", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 
hp_year_predict <- ggpredict(hyear, terms = "year", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 

## year plot
hp_year_predict$year <- c(1992:2022)
hp_year_predict$year <- as.numeric(hp_year_predict$year) 
h_year$year <- c(1992:2022)

hp_years <- ggplot (hp_year_predict, aes(year, predicted, colour = NULL)) +
  #geom_point(mapping = aes(year, n), data = hp_year, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab (" ") +
  xlab (" ") + 
  ggtitle("Harbour Porpoises")+
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_hp, x = 1994, y = 170, ysize = 23, alpha = 1)


## seasonal plot
hp_season <- ggplot (hp_month_predict, aes(x, predicted, colour = NULL)) +
  #geom_point(mapping = aes(month, n), data = hp_count, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab ("") +
  xlab ("") + 
  ggtitle("Harbour Porpoises") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_hp, x = 1.75, y = 11, ysize = 1.5, alpha = 1)


## SENSITIVTY CHECKS ##

## Make dataframe
hp.group <- count(df, spec.r, group.r, year, month) %>%
  filter(group.r == "Harbour porpoise")

hp.group$spec.r <- as.factor(hp.group$spec.r)
hp.group <- filter(hp.group, spec.r != "Cetacean (indeterminate species)")


## TIME INFLUENCE
# Add a pseudo year with strandings to see if trend holds

## Make df
random <- data.frame(year=2023, month= c(1:12), n=rep(4, each = 12))
yeartest.hp <- merge(random, hp.group, all = TRUE)

## Run model 
year.mod.hp <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                        data = yeartest.hp, family = poisson(link="log"), 
                        correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.mod.hp$gam)
summary(year.mod.hp$lme) #phi = -0.038

## Results
plot(year.mod.hp$gam, residuals = TRUE) 
summary(year.mod.hp$gam)
## Trend holds. 

# Remove a year and see if trend holds
## Remove year
hp.less <- filter(hp.group, year != "2018")

## Run model
year.less.hp <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                         data = hp.less, family = poisson(link="log"), 
                         correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.less.hp$gam)
summary(year.less.hp$lme) #phi = -0.5
overdispersion_gamm(year.less.hp, hp.less) #1.19

## Results 
summary(year.less.hp$gam)
plot(year.less.hp$gam, residuals = TRUE) 
# Trend holds, oscillation slightly less obvious but still remains

## Pelagic Dolphins######################################################################################################
# Pelagic Dolphins
## get best model, run diagnostics and visuals

## extract df from results list
pelagic_df <- group_list[["Pelagic Dolphins"]]

## run best fit model according to AIC table + includes AC structure
p_model <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML",
                 data = pelagic_df, family = poisson(link = "log"),
                 correlation = corARMA(form = ~ 1|year, p = 1))

## conduct diagnostics 
gam.check(p_model$gam) #ok
summary(p_model$lme) # phi = 0.08
summary(p_model$gam) # r2 = .11 (its an estimate for a gamm)
overdispersion_gamm(p_model, pelagic_df) # 0.45

## test for residual autocorrelation
## change layout so you can see both plots at once
layout(matrix(1:2, ncol = 2))

## extract the normalised residuals
res <- resid(p_model$lme, type = "normalized")

## run acf and partial acf plots 
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors") # looks like it is captured well

## return layout to normal
layout(1) 

## look at results
plot(p_model$gam, resid = TRUE) ##looks good
summary(p_model$gam) #both significant

## visualise
## make a df count for year
p_year <- count(df, group.r, year)
p_year <- filter(p_year, group.r == "Pelagic Dolphins")

## make a simplistic year model for plotting
pyear <- gamm(n ~ s(year), data = p_year, family = poisson(link = "log"), method = "REML", 
              correlation = corARMA(form = ~ year, p = 1))

## get pilot whale phylopic for plot
p_uuid <- get_uuid(name = "Globicephala melas")
p <- get_phylopic(uuid = p_uuid, format = "vector")

## make porpoise stranded
dead_p <- rotate_phylopic(p, angle = 180)

## generate predictions from model to plot
p_month_predict <- ggpredict(p_model$gam, terms = "month", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 
p_year_predict <- ggpredict(pyear$gam, terms = "year", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 

## year plot
p_year_predict$year <- c(1992:2022)
p_year_predict$year <- as.numeric(p_year_predict$year) 
p_year$year <- c(1992:2022)

pel_year <- ggplot (p_year_predict, aes(year, predicted)) +
  #geom_point(mapping = aes(year, n), data = pelagic_year, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab (" ") +
  xlab (" ") + 
  ggtitle("Pelagic Dolphins") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_p, x = 1995, y = 53, ysize = 7, alpha = 1)

## seasonal plot
pel_season <- ggplot (p_month_predict, aes(x, predicted, colour = NULL)) +
  #geom_point(mapping = aes(month, n), data = pelagic, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab ("") +
  xlab ("") + 
  ggtitle("Pelagic Dolphins") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_p, x = 2, y = 5, ysize = 0.75, alpha = 1)


## SENSITIVITY CHECKS ##

# PER SPECIES
# species plots on line 179 plots shows trends are similar in pelagic dolphins
# though strength of trend varies a bit so check for differences 

## Make dataset
pelagic.group <- count(df, spec.r, group.r, year, month) %>%
  filter(group.r == "Pelagic Dolphins")

pelagic.group <- filter(pelagic.group, spec.r != "Cetacean (indeterminate species)")
pelagic.group$spec.r <- as.factor(pelagic.group$spec.r)

## Run basic model to observe differences between species
group.mod <- gamm(n ~ s(month, bs = "cc", k = 12, by = spec.r) + spec.r*month + s(year), knots = knots, method = "REML", 
                  data = pelagic.group, family = poisson(link="log"), 
                  correlation = corARMA(form = ~ 1|year, p = 1))

## Look at diagnostics 
gam.check(group.mod$gam)
summary(group.mod$lme) #0.42
overdispersion_gamm(group.mod, pelagic.group)

## Results
summary(group.mod$gam)
plot(group.mod$gam) 
## No difference between species.

## Looking at trends shows striped dolphins are slightly different, but still retain 

## TIME INFLUENCE
pelagic <- count(df, group.r, year, month) %>%
  filter(group.r == "Pelagic Dolphins")

# Add a pseudo year with strandings to see if trend holds

## Make df
set.seed(123) 
random <- data.frame(year=2023, month= c(1:12), n = sample(1:5, 12, replace = TRUE))
yeartest.pelagic <- merge(random, pelagic, all = TRUE)

## Run model 
year.mod.pelagic <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                        data = yeartest.pelagic, family = poisson(link="log"), 
                        correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.mod.pelagic$gam)
summary(year.mod.pelagic$lme) #phi = -0.038

## Results
plot(year.mod.pelagic$gam) 
summary(year.mod.pelagic$gam)
## Trend holds


#Remove a year and see if trend holds
## Remove year
pelagic.less <- filter(pelagic, year != "2018")

## Run model
year.less.pelagic <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                         data = pelagic.less, family = poisson(link="log"), 
                         correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.less.pelagic$gam)
summary(year.less.pelagic$lme) #phi = -0.5
overdispersion_gamm(year.less.pelagic, pelagic.less) #1.19

## Results 
summary(year.less.pelagic$gam)
plot(year.less.pelagic$gam) 
# Trend holds.


## Deep Divers###########################################################################################################
# Deep Divers
## get best model, run diagnostics and visuals
hist(deep_df$n)

## extract df from results list
deep_df <- group_list[["Deep Divers"]]
deep_df <- filter(deep_df, year != "2018")

## run best fit model according to AIC table
d_model <- gam(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
               data = deep_df, family = ziP())

## conduct diagnostics 
gam.check(d_model) #ok
overdispersion_gam(d_model) # 0.04

## test for residual autocorrelation
## change layout so you can see both plots at once
layout(matrix(1:2, ncol = 2))

## extract Pearson's residuals
res <- resid(d_model)

## run acf and partial acf plots 
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors") # looks like it is captured well

## return layout to normal
layout(1) 

## look at results
plot(d_model, resid = TRUE) ##looks good
summary(d_model$gam) #both significant

## visualise
## make a df count for year
d_year <- count(df, group.r, year)
d_year <- filter(d_year, group.r == "Deep Divers")
d_year <- filter(d_year, year != "2018")

## make a simplistic year model for plotting
dyear <- gam(n ~ s(year), data = d_year, family = poisson(link = "log"), method = "REML")

## get sperm whale phylopic for plot
d_uuid <- get_uuid(name = "Physeter macrocephalus")
d <- get_phylopic(uuid = d_uuid, format = "vector")

## make sperm whale stranded
dead_d <- rotate_phylopic(d, angle = 180)

## generate predictions from model to plot
d_month_predict <- ggpredict(d_model, terms = "month", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 
d_year_predict <- ggpredict(dyear, terms = "year", ci.lvl = 0.95, type = "fixed", back.transform = TRUE) 

## year plot
d_year_predict$year <- d_year$year
d_year_predict$year <- as.numeric(d_year_predict$year) 

d_years <- ggplot (d_year_predict, aes(year, predicted)) +
  #geom_point(mapping = aes(year, n), data = d_year, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab (" ") +
  xlab (" ") + 
  ggtitle("Deep Divers") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_d, x = 1996, y = 15, ysize = 2.5, alpha = 1)

## seasonal plot
d_season <- ggplot (d_month_predict, aes(x, predicted, colour = NULL)) +
  #geom_point(mapping = aes(month, n), data = deep.less, col = "grey90") +
  geom_line(aes(linetype = NULL), linewidth = 0.75)  + 
  geom_ribbon(aes(ymin =conf.low, ymax = conf.high), alpha = 0.2) +
  ylab ("") +
  xlab ("") + 
  ggtitle("Deep Divers") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() +
  add_phylopic(dead_d, x = 2, y = 0.75, ysize = 0.3, alpha = 1)

## SENSITIVITY CHECKS ##

# PER SPECIES
# species plots on line 179 plots shows trends are similar in deep divers
# though strength of trend varies a bit so check for differences 

## Make dataset
deep.group <- count(df, spec.r, group.r, year, month) %>%
  filter(group.r == "Deep Divers")

deep.group <- filter(deep.group, spec.r != "Cetacean (indeterminate species)")
deep.group$spec.r <- as.factor(deep.group$spec.r)

## Run basic model to observe differences between species
group.mod <- gamm(n ~ s(month, bs = "cc", k = 12, by = spec.r) + spec.r*month + s(year), knots = knots, method = "REML", 
                  data = deep.group, family = poisson(link="log"), 
                  correlation = corARMA(form = ~ 1|year, p = 1))

## Look at diagnostics 
gam.check(group.mod$gam)
summary(group.mod$lme) #0.42
overdispersion_gamm(group.mod, deep.group)

## Results
summary(group.mod$gam)
plot(group.mod$gam) 
## No difference between species.
## Looking at trends shows cuvier's  hold the trend the most. 

## TIME INFLUENCE
deep <- count(df, group.r, year, month) %>%
  filter(group.r == "Deep Divers")

# Add a pseudo year with strandings to see if trend holds

## Make df
set.seed(123) 
random <- data.frame(year=2023, month= c(1:12), n = sample(1:8, 12, replace = TRUE))
yeartest.deep <- merge(random, deep, all = TRUE)

## Run model 
year.mod.deep <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                         data = yeartest.deep, family = poisson(link="log"), 
                         correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.mod.deep$gam)
summary(year.mod.deep$lme) #phi = -0.038

## Results
plot(year.mod.deep$gam) 
summary(year.mod.deep$gam)
## Trend hard to see in visuals


#Remove a year and see if trend holds
## Remove year
deep.less <- filter(deep, year != "2018")

## Run model
year.less.deep <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year), knots = knots, method = "REML", 
                          data = deep.less, family = poisson(link="log"), 
                          correlation = corARMA(form = ~ 1|year, p = 1))

## Check diagnostics 
gam.check(year.less.deep$gam)
summary(year.less.deep$lme) #phi = -0.5
overdispersion_gamm(year.less.deep, deep.less) #1.19

## Results 
summary(year.less.deep$gam)
plot(year.less.deep$gam) 
# Trend holds.

# Figures################################################################################################################
# Figures


## all seasonal plots
tiff('all_seasonal.tiff', units="in", width=7, height=5, res=1000)

grid.arrange(bal_season,cd_season, d_season, pel_season, hp_season, ncol=2, left = "Seasonal Count of Strandings", bottom = "Month")

dev.off()

## all annual plots
tiff('all_annual.tiff', units="in", width=7, height=5, res=1000)

grid.arrange(bal_year,cd_years, d_years, pel_year, hp_years, ncol=2, left = "Annual Count of Strandings", bottom = "Year")

dev.off()

#####################################################################
# Demographic Trends ----
## SEX DATASETS ----
## initialise empty list to store results of each group
sex_list <- list()

## loop through each group to generate a count with zeroes added in for month and year where relevant
for(group in unique_groups) {
  # Filter the dataframe for the current group
  group_df <- filter(df, df$group.r == group)
  
  # Add zero data
  tm <- data.frame(year=c(1992:2022), month= rep(c(1:12), each = 31), n=c(0), sex = "M")
  tf <- data.frame(year=c(1992:2022), month= rep(c(1:12), each = 31), n=c(0), sex = "F")
  
  # Count the occurrences for the current group
  group_m <- count(group_df, year, month, sex) %>%
    filter(sex != "U")
  
  # Merge with zero data
  group_m <- merge(tm, group_m, all = TRUE)
  group_m <- merge(tf, group_m, all = TRUE)
  
  # Handle duplicates
  group_m$duplicate <- duplicated(group_m[c("year", "month", "sex")])
  group_m$duplicate <- ifelse(duplicated(group_m[,c("year", "month", "sex")], fromLast=T), yes=T, no = group_m$duplicate)
  group_result <- group_m[!(group_m$duplicate == TRUE & group_m$n == 0),] 
  
  # make sex a factor 
  group_result$sex <- as.factor(group_result$sex)
  
  # Store the result in the list with the group name
  sex_list[[group]] <- group_result
}

## AGE DATASETS ----
## initialise empty list to store results of each group
age_list <- list()

## loop through each group to generate a count with zeroes added in for month and year where relevant
for(group in unique_groups) {
  # Filter the dataframe for the current group
  group_df <- filter(df, df$group.r == group)
  
  # Add zero data
  tj <- data.frame(year=c(1992:2022), month= rep(c(1:12), each = 31), n=c(0), age.r = "juvenile")
  ta <- data.frame(year=c(1992:2022), month= rep(c(1:12), each = 31), n=c(0), age.r = "adult")
  tn <- data.frame(year=c(1992:2022), month= rep(c(1:12), each = 31), n=c(0), age.r = "neonate")
  
  
  # Count the occurrences for the current group
  group_m <- count(group_df, year, month, age.r) %>%
    filter(age.r != "unknown")
  
  # Merge with zero data
  group_m <- merge(tj, group_m, all = TRUE)
  group_m <- merge(ta, group_m, all = TRUE)
  group_m <- merge(tn, group_m, all = TRUE)
  
  # Handle duplicates
  group_m$duplicate <- duplicated(group_m[c("year", "month", "age.r")])
  group_m$duplicate <- ifelse(duplicated(group_m[,c("year", "month", "age.r")], fromLast=T), yes=T, no = group_m$duplicate)
  group_result <- group_m[!(group_m$duplicate == TRUE & group_m$n == 0),] 
  
  # Store the result in the list with the group name
  age_list[[group]] <- group_result
}

## MODEL COMPARISONS ----
# Initialize a list to store the AIC values for each model
sex_aic_results <- list()

# Loop through each dataframe in the group_list
for(group in names(sex_list)) {
  
  # Extract the current dataframe
  group_df <- sex_list[[group]]
  
  # Fit each model
  
  ## GLM
  sex_glm <- glm(n ~ month*sex + year*sex, data = group_df, family = poisson(link = "log"))
  
  ## GAM with no autocorrelation
  sex_gam_no_ac <- gam(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, 
                         knots = knots, method = "REML", data = group_df, 
                         family = poisson(link = "log"))
  
  ## GAMM with autocorrelation (AC 1)
  sex_gamm_ac1 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, 
                         knots = knots, method = "REML", data = group_df, 
                         family = poisson(link = "log"), 
                         correlation = corARMA(form = ~ 1|year, p = 1))
  
  ## GAMM with autocorrelation (AC 1) no interaction of sex
  sex_gamm_ac1_no_int <- gamm(n ~ s(month, bs = "cc", k = 12) + s(year) + sex, 
                       knots = knots, method = "REML", data = group_df, 
                       family = poisson(link = "log"), 
                       correlation = corARMA(form = ~ 1|year, p = 1))
  
  ## GAMM with random effect of year
  sex_gamm_re <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, bs = "re") + sex, 
                        knots = knots, method = "REML", data = group_df, 
                        family = poisson(link = "log"), 
                        correlation = corARMA(form = ~ 1|year, p = 1))
  
  ## GAMM with negative binomial (NB) with autocorrelation
  sex_gamm_nb_ac <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, 
                           knots = knots, method = "REML", data = group_df, 
                           family = negative.binomial(1), 
                           correlation = corARMA(form = ~ 1|year, p = 1))
  
  # Store the AIC values for each model in a list
  aic_values <- AIC(sex_glm, sex_gam_no_ac, 
                    sex_gamm_ac1$lme, sex_gamm_ac1_no_int$lme, 
                    sex_gamm_re$lme, sex_gamm_nb_ac$lme)
  
  # Add the AIC values to the results list with the group name
  sex_aic_results[[group]] <- aic_values
}

# Convert the results list to a data frame for easy viewing
sex_aic_table <- do.call(rbind, sex_aic_results)

## AIC results
## Baleen whale model = GAM no AC
## Common dolphin model = GAM no AC
## Deep diver model = GAM no AC
## Pelagic dolphin model = GAM no AC
## Harbour porpoise model = GAM AC, no interaction

## BALEEN WHALES ----
## SEX TRENDS ##

## get df with only annual counts. 
baleen.sex.y <- count(df, group.r, year, sex) %>%
  filter(group.r == "Baleen Whales")
baleen.sex.y <- filter(baleen.sex.y, baleen.sex.y$sex !="U")

#year
## Run model 
  ## GLMs of year and sex interaction
  baleen.sex.year <- glm(n ~ sex*year,data = baleen.sex.y, family = poisson(link = "log"))
  
  ## GLMs with interaction and re of year
  baleen.sex.year2 <- glm(n ~ sex*year + (1|year),data = baleen.sex.y, family = poisson(link = "log"))

#compare
AIC(baleen.sex.year, baleen.sex.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(baleen.sex.year) ##ok 

# results
summary(baleen.sex.year) ## No difference in sex increases per year. 

#month
## get df 
baleen.sex.s <- sex_list[["Baleen Whales"]]
baleen.sex.s$sex <- as.factor(baleen.sex.s$sex)

## Run models

  ## GAM with AC 
  baleen.sex.season <- gam(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                          data = baleen.sex.s, family = poisson(link="log"))
  
  ## GAMM with AC 
  baleen.sex.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                 data = baleen.sex.s, family = poisson(link="log"), 
                 correlation = corARMA(form = ~ 1|year, p = 1))
  
#compare
AIC(baleen.sex.season, baleen.sex.season2$lme)

# diagnostics
gam.check(sexmonth$gam) #not a great fit
overdispersion_gamm(sexmonth, baleen_countsm) #5.7, over dispersed. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(baleen.sex.season$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC plot.

  ## GAMM with AC and nb distribution 
  baleen.sex.season3 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                   data = baleen.sex.s, family = negative.binomial(1), 
                   correlation = corARMA(form = ~ 1|year, p = 1))

#diagnostics
gam.check(baleen.sex.season3$gam) #not a great fit but better. 

#compare 
AIC(baleen.sex.season2$lme, baleen.sex.season3$lme) ## model 1 is better but clear overdispersion so use nb

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(baleen.sex.season2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure. 

plot(res) #structure
summary(baleen.sex.season3$lme) # phi = -.1

## Results 
summary(baleen.sex.season3$gam)
plot(baleen.sex.season3$gam)
## Seasonality of males and females the same, no difference in males/females


## AGE ##

## get df with only annual counts. 
baleen.age.y <- count(df, group.r, year, age.r) %>%
  filter(group.r == "Baleen Whales")
baleen.age.y <- filter(baleen.age.y, baleen.age.y$age.r !="unknown")

#year
## Run model 
## GLMs of year and age interaction
baleen.age.year <- glm(n ~ age.r*year,data = baleen.age.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
baleen.age.year2 <- glm(n ~ age.r*year + (1|year),data = baleen.age.y, family = poisson(link = "log"))

#compare
AIC(baleen.age.year, baleen.age.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(baleen.age.year) ##good

# results
summary(baleen.age.year) ## juveniles increasing more than adults, neonates decreasing. 


#month
## get df 
baleen.age.s <- age_list[["Baleen Whales"]]
baleen.age.s$age.r <- as.factor(baleen.age.s$age.r)

## Run models

## GAM with AC 
baleen.age.season <- gam(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                         data = baleen.age.s, family = poisson(link="log"))

## GAMM with AC 
baleen.age.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = baleen.age.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(baleen.age.season, baleen.age.season2$lme)

# diagnostics
gam.check(baleen.age.season2$gam) #not a great fit
overdispersion_gamm(baleen.age.season2, baleen.age.s) #33.3, over dispersed. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(baleen.age.season$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC model

## GAMM with AC and nb distribution 
baleen.age.season3 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age, knots = knots, method = "REML", 
                           data = baleen.age.s, family = negative.binomial(1), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#diagnostics
gam.check(baleen.age.season3$gam) #not a great fit but better. 

#compare 
AIC(baleen.age.season2$lme, baleen.age.season3$lme) ## the same, but use nb as overdispersed and fit better

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(baleen.age.season3$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure. 

plot(res) #structure
summary(baleen.age.season3$lme) # phi = 0.04

## Results 
summary(baleen.age.season3$gam)
plot(baleen.age.season3$gam)
## Adults and juveniles in the summer, neonate no trend. 

## Visualise seasonal age trends
baleen.age.plot <- ggplot (baleen.age.s, aes(month, n, col=age.r)) +
  geom_smooth(method = "gam")  + 
  ylab ("") +
  xlab ("") + 
  labs(color = "Age Class") +
  ggtitle("Baleen Whales") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() 

## COMMON DOLPHINS ----
## SEX TRENDS ##

## get df with only annual counts. 
common.sex.y <- count(df, group.r, year, sex) %>%
  filter(group.r == "Common Dolphin")
common.sex.y <- filter(common.sex.y, common.sex.y$sex !="U")

#year
## Run model 
## GLMs of year and sex interaction
common.sex.year <- glm(n ~ sex*year,data = common.sex.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
common.sex.year2 <- glm(n ~ sex*year + (1|year),data = common.sex.y, family = poisson(link = "log"))

#compare
AIC(common.sex.year, common.sex.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(common.sex.year) ##ok 

# results
summary(common.sex.year) ## No difference in sex increases per year. 

#month
## get df 
common.sex.s <- sex_list[["Common Dolphin"]]
common.sex.s$sex <- as.factor(common.sex.s$sex)

## Run models

## GAM with AC 
common.sex.season <- gam(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                         data = common.sex.s, family = poisson(link="log"))

## GAMM with AC 
common.sex.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                           data = common.sex.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(common.sex.season, common.sex.season2$lme)

# diagnostics
gam.check(common.sex.season) #not a great fit
overdispersion_gamm(common.sex.season2, common.sex.s) #6.22, over dispersed. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(common.sex.season2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC plot.

## GAMM with AC and nb distribution 
common.sex.season3 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                           data = common.sex.s, family = negative.binomial(1), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#diagnostics
gam.check(common.sex.season3$gam) #not a great fit but better. 

#compare 
AIC(common.sex.season2$lme, common.sex.season3$lme) ## model 3 is best

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(common.sex.season3$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure. 

plot(res) #structure
summary(common.sex.season3$lme) # phi = -.02, fine

## Results 
summary(common.sex.season3$gam)
plot(common.sex.season3$gam)
## Seasonality of males and females the same


## AGE ##

## get df with only annual counts. 
common.age.y <- count(df, group.r, year, age.r) %>%
  filter(group.r == "Common Dolphin")
common.age.y <- filter(common.age.y, common.age.y$age.r !="unknown")

#year
## Run model 
## GLMs of year and age interaction
common.age.year <- glm(n ~ age.r*year,data = common.age.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
common.age.year2 <- glm(n ~ age.r*year + (1|year),data = common.age.y, family = poisson(link = "log"))

#compare
AIC(common.age.year, common.age.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(common.age.year) ##good

# results
summary(common.age.year) ## juveniles increasing more than adults, neonates increasing but non significantly.  


#month
## get df 
common.age.s <- age_list[["Common Dolphin"]]
common.age.s$age.r <- as.factor(common.age.s$age.r)

## Run models

## GAM with AC 
common.age.season <- gam(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                         data = common.age.s, family = poisson(link="log"))

## GAMM with AC 
common.age.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = common.age.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(common.age.season, common.age.season2$lme)

# diagnostics
gam.check(common.age.season) #not a great fit
overdispersion_gamm(common.age.season2, common.age.s) #6.52, over dispersed. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(common.age.season$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC model

## GAMM with AC and nb distribution 
common.age.season3 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = common.age.s, family = negative.binomial(1), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#diagnostics
gam.check(common.age.season3$gam) #not a great fit but better. 

#compare 
AIC(common.age.season2$lme, common.age.season3$lme) ## the same, but use nb as overdispersed and fit better

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(common.age.season3$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure. 

plot(res) #structure
summary(common.age.season3$lme) # phi = 0.04

## Results 
summary(common.age.season3$gam)
plot(common.age.season3$gam)
## Seasonality of adults and juveniles the same, neonates in the summer. 

## Visualise seasonal age trend

common.age.plot <- ggplot (common.age.s, aes(month, n, col=age.r)) +
  geom_smooth(method = "gam")  + 
  ylab ("") +
  xlab ("") + 
  ggtitle("Common Dolphins") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() 

## DEEP DIVERS ----
## SEX TRENDS ##

## get df with only annual counts. 
deep.sex.y <- count(df, group.r, year, sex) %>%
  filter(group.r == "Deep Divers")
deep.sex.y <- filter(deep.sex.y, deep.sex.y$sex !="U")

#year
## Run model 
## GLMs of year and sex interaction
deep.sex.year <- glm(n ~ sex*year,data = deep.sex.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
deep.sex.year2 <- glm(n ~ sex*year + (1|year),data = deep.sex.y, family = poisson(link = "log"))

#compare
AIC(deep.sex.year, deep.sex.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(deep.sex.year) ##ok 

# results
summary(deep.sex.year) ## No difference in sex increases per year. 

#month
## get df 
deep.sex.s <- sex_list[["Deep Divers"]]
deep.sex.s$sex <- as.factor(deep.sex.s$sex)

## Run models

## GAM with AC 
deep.sex.season <- gam(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                         data = deep.sex.s, family = poisson(link="log"))

## GAMM with AC 
deep.sex.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                           data = deep.sex.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(deep.sex.season, deep.sex.season2$lme)

# diagnostics
gam.check(deep.sex.season) #not a great fit
overdispersion_gamm(deep.sex.season2, deep.sex.s) #6.69, over dispersed. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(deep.sex.season2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC plot.

## GAMM with AC and nb distribution 
deep.sex.season3 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                           data = deep.sex.s, family = negative.binomial(1), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#diagnostics
gam.check(deep.sex.season3$gam) #not a great fit but better. 

#compare 
AIC(deep.sex.season2$lme, deep.sex.season3$lme) ## model 1 is better but clear overdispersion so use nb

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(deep.sex.season3$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure. 

plot(res) #structure
summary(deep.sex.season3$lme) # phi = 0.12

## Results 
summary(deep.sex.season3$gam)
plot(deep.sex.season3$gam)
## Seasonality of males and females the same, but male trend significant whereas female not. 


## AGE ##
## get df with only annual counts. 
deep.age.y <- count(df, group.r, year, age.r) %>%
  filter(group.r == "Deep Divers")
deep.age.y <- filter(deep.age.y, deep.age.y$age.r !="unknown")

#year
## Run model 
## GLMs of year and age interaction
deep.age.year <- glm(n ~ age.r*year,data = deep.age.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
deep.age.year2 <- glm(n ~ age.r*year + (1|year),data = deep.age.y, family = poisson(link = "log"))

#compare
AIC(deep.age.year, deep.age.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(deep.age.year) ##good

# results
summary(deep.age.year) ## juveniles increasing more than adults, neonates no trend - no neonates.


#month
## get df 
deep.age.s <- age_list[["Deep Divers"]]
deep.age.s$age.r <- as.factor(deep.age.s$age.r)

## Run models

## GAM with AC 
deep.age.season <- gam(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                         data = deep.age.s, family = poisson(link="log"))

## GAMM with AC 
deep.age.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = deep.age.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(deep.age.season, deep.age.season2$lme)

# diagnostics
gam.check(deep.age.season) #not a great fit
overdispersion_gamm(deep.age.season2, deep.age.s) #9.74, over dispersed. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(deep.age.season2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC model

## GAMM with AC and nb distribution 
deep.age.season3 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = deep.age.s, family = negative.binomial(1), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#diagnostics
gam.check(deep.age.season3$gam) #not a great fit but better. 

#compare 
AIC(deep.age.season2$lme, deep.age.season3$lme) ## mod2 better, but use nb as overdispersed and fit better

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(deep.age.season3$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure. 

plot(res) #structure
summary(deep.age.season3$lme) # phi = 0.03

## Results 
summary(deep.age.season3$gam)
plot(deep.age.season3$gam)
## Seasonality the same for all. 

## Visualise seasonal age trend

deep.age.plot <- ggplot (deep.age.s, aes(month, n, col=age.r)) +
  geom_smooth(method = "gam")  + 
  ylab ("") +
  xlab ("") + 
  ggtitle("Deep Divers") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() 

## HARBOUR PORPOISE ----
## SEX TRENDS ##

## get df with only annual counts. 
hp.sex.y <- count(df, group.r, year, sex) %>%
  filter(group.r == "Harbour porpoise")
hp.sex.y <- filter(hp.sex.y, hp.sex.y$sex !="U")

#year
## Run model 
## GLMs of year and sex interaction
hp.sex.year <- glm(n ~ sex*year,data = hp.sex.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
hp.sex.year2 <- glm(n ~ sex*year + (1|year),data = hp.sex.y, family = poisson(link = "log"))

#compare
AIC(hp.sex.year, hp.sex.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(hp.sex.year) ##ok 

# results
summary(hp.sex.year) ## No difference in sex increases per year. 

#month
## get df 
hp.sex.s <- sex_list[["Harbour porpoise"]]
hp.sex.s$sex <- as.factor(hp.sex.s$sex)

## Run models

## GAM with AC 
hp.sex.season <- gam(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                         data = hp.sex.s, family = poisson(link="log"))

## GAMM with AC 
hp.sex.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                           data = hp.sex.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(hp.sex.season, hp.sex.season2$lme) #model 2 better

# diagnostics
overdispersion_gamm(hp.sex.season2, hp.sex.s) #0.7, not overdispersed
gam.check(hp.sex.season2$gam)
summary(hp.sex.season2$lme) # phi = 0.06

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(hp.sex.season2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC plot.

## Results 
summary(hp.sex.season2$gam)
plot(hp.sex.season2$gam)
## Male smooth over spring, females have dip in May, no significant difference though. 


## AGE ##

## get df with only annual counts. 
hp.age.y <- count(df, group.r, year, age.r) %>%
  filter(group.r == "Harbour porpoise")
hp.age.y <- filter(hp.age.y, hp.age.y$age.r !="unknown")

#year
## Run model 
## GLMs of year and age interaction
hp.age.year <- glm(n ~ age.r*year,data = hp.age.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
hp.age.year2 <- glm(n ~ age.r*year + (1|year),data = hp.age.y, family = poisson(link = "log"))

#compare
AIC(hp.age.year, hp.age.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(hp.age.year) ##good

# results
summary(hp.age.year) ## no difference. 


#month
## get df 
hp.age.s <- age_list[["Harbour porpoise"]]
hp.age.s$age.r <- as.factor(hp.age.s$age.r)

## Run models

## GAM with AC 
hp.age.season <- gam(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                         data = hp.age.s, family = poisson(link="log"))

## GAMM with AC 
hp.age.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = hp.age.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(hp.age.season, hp.age.season2$lme)

# diagnostics
gam.check(hp.age.season) #good fit
overdispersion_gamm(hp.age.season2, hp.age.s) #6.51, over dispersed. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(hp.age.season2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC model

## GAMM with AC and nb distribution 
hp.age.season3 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = hp.age.s, family = negative.binomial(1), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#diagnostics
gam.check(hp.age.season3$gam) #not a great fit but better. 

#compare 
AIC(hp.age.season2$lme, hp.age.season3$lme) ## model 3 better

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(hp.age.season3$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure. 

plot(res) #structure
summary(hp.age.season3$lme) # phi = 0.02

## Results 
summary(hp.age.season3$gam)
plot(hp.age.season3$gam)
## Juvenile spring, neonate summer, adults no trend/ 

## Visualise seasonal age trend
hp.age.plot <- ggplot (hp.age.s, aes(month, n, col=age.r)) +
  geom_smooth(method = "gam")  + 
  ylab ("") +
  xlab ("") + 
  ggtitle("Harbour Porpoises") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() 

## PELAGIC DOLPHINS ----
## SEX TRENDS ##

## get df with only annual counts. 
pelagic.sex.y <- count(df, group.r, year, sex) %>%
  filter(group.r == "Pelagic Dolphins")
pelagic.sex.y <- filter(pelagic.sex.y, pelagic.sex.y$sex !="U")

#year
## Run model 
## GLMs of year and sex interaction
pelagic.sex.year <- glm(n ~ sex*year,data = pelagic.sex.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
pelagic.sex.year2 <- glm(n ~ sex*year + (1|year),data = pelagic.sex.y, family = poisson(link = "log"))

#compare
AIC(pelagic.sex.year, pelagic.sex.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(pelagic.sex.year) ##ok 

# results
summary(pelagic.sex.year) ## No difference in sex increases per year. 

#month
## get df 
pelagic.sex.s <- sex_list[["Pelagic Dolphins"]]
pelagic.sex.s$sex <- as.factor(pelagic.sex.s$sex)

## Run models

## GAM with AC 
pelagic.sex.season <- gam(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                         data = pelagic.sex.s, family = poisson(link="log"))

## GAMM with AC 
pelagic.sex.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = sex) + s(year, by = sex) + sex, knots = knots, method = "REML", 
                           data = pelagic.sex.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(pelagic.sex.season, pelagic.sex.season2$lme)

# diagnostics
gam.check(pelagic.sex.season) #not a great fit
overdispersion_gamm(pelagic.sex.season2, pelagic.sex.s) #1.27, fine. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(pelagic.sex.season2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC plot.


#diagnostics
gam.check(pelagic.sex.season2$gam) #not a great fit but better. 
summary(pelagic.sex.season2$lme) # phi = -.03

## Results 
summary(pelagic.sex.season2$gam)
plot(pelagic.sex.season2$gam)
## Seasonality of males and females the same, no difference in males/females


## AGE ##
## get df with only annual counts. 
pelagic.age.y <- count(df, group.r, year, age.r) %>%
  filter(group.r == "Pelagic Dolphins")
pelagic.age.y <- filter(pelagic.age.y, pelagic.age.y$age.r !="unknown")

#year
## Run model 
## GLMs of year and age interaction
pelagic.age.year <- glm(n ~ age.r*year,data = pelagic.age.y, family = poisson(link = "log"))

## GLMs with interaction and re of year
pelagic.age.year2 <- glm(n ~ age.r*year + (1|year),data = pelagic.age.y, family = poisson(link = "log"))

#compare
AIC(pelagic.age.year, pelagic.age.year2) ## the same, stick with simpler model. 

# diagnostics 
plot(pelagic.age.year) ##good

# results
summary(pelagic.age.year) ## juveniles increasing more than adults, neonates increasing the most.  


#month
## get df 
pelagic.age.s <- age_list[["Pelagic Dolphins"]]
pelagic.age.s$age.r <- as.factor(pelagic.age.s$age.r)

## Run models

## GAM with AC 
pelagic.age.season <- gam(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                         data = pelagic.age.s, family = poisson(link="log"))

## GAMM with AC 
pelagic.age.season2 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = pelagic.age.s, family = poisson(link="log"), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#compare
AIC(pelagic.age.season, pelagic.age.season2$lme)

# diagnostics
gam.check(pelagic.age.season) #ok fit
overdispersion_gamm(pelagic.age.season2, pelagic.age.s) #14.05, over dispersed. 

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(pelagic.age.season2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure, keep AC model

## GAMM with AC and nb distribution 
pelagic.age.season3 <- gamm(n ~ s(month, bs = "cc", k = 12, by = age.r) + s(year, by = age.r) + age.r, knots = knots, method = "REML", 
                           data = pelagic.age.s, family = negative.binomial(1), 
                           correlation = corARMA(form = ~ 1|year, p = 1))

#diagnostics
gam.check(pelagic.age.season3$gam) #not a great fit but better. 

#compare 
AIC(pelagic.age.season2$lme, pelagic.age.season3$lme) ## the same, but use nb as overdispersed and fit better

## Check ACF plots
layout(matrix(1:2, ncol = 2))
res <- resid(pelagic.age.season3$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(1) errors")
pacf(res, lag.max = 36, main = "pACF- AR(1) errors")
layout(1) #retains some correlation even with structure. 

plot(res) #structure
summary(pelagic.age.season3$lme) # phi = 0.04

## Results 
summary(pelagic.age.season3$gam)
plot(pelagic.age.season3$gam)
## All have summer peak, juveniles also have little peak in jan. 

## Visualise seasonal age trend
pelagic.age.plot <- ggplot (pelagic.age.s, aes(month, n, col=age.r)) +
  geom_smooth(method = "gam")  + 
  ylab ("") +
  xlab ("") + 
  ggtitle("Pelagic Dolphins") +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_y_continuous(breaks=pretty_breaks()) +
  theme_classic() 

## DEMOGRAPHIC FIGURES ----
## get legend from one plot
legend <- get_legend(baleen.age.plot)

## remove legends
baleen.age.plot <- baleen.age.plot + theme(legend.position = "none")
common.age.plot <- common.age.plot + theme(legend.position = "none")
deep.age.plot <- deep.age.plot + theme(legend.position = "none")
hp.age.plot <- hp.age.plot + theme(legend.position = "none")
pelagic.age.plot <- pelagic.age.plot + theme(legend.position = "none")

# Arrange the plots into a grid
combined_plots <- plot_grid(baleen.age.plot, common.age.plot, 
                            deep.age.plot, hp.age.plot, pelagic.age.plot, 
                            ncol = 2)  # You can adjust the ncol to your preference

# Combine the grid of plots with the common legend below
all_seasonal_age <- plot_grid(combined_plots, legend, 
                        ncol = 2,  
                        rel_widths = c(1, 0.1))



## all seasonal age plots
tiff('all_seasonal_age.tiff', units="in", width=9, height=5, res=1000)

all_seasonal_age

dev.off()

#####################################################################
# Spatiotemporal Trends ----
## SPATIAL DATASET ----
# set up a regional dataset
# clean up df of NAs in long/lat
spatial <- filter(df, latwgs84 != "NA")
spatial <- filter(df, longwgs84 != "NA")
spatial <- filter(df, longwgs84 != "#VALUE!")

# Convert data.frame to sf object to get points
cets_sf <- st_as_sf(spatial, coords = c("longwgs84", "latwgs84"),crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Download 12 regions shp file
regions <- st_read("C:/Users/Rachel Lennon/OneDrive - University of Glasgow/PhD/Projects/Trends/spatial/regions/Region shp/Normal/spatialregions.shp")

# Map data points within each region
Points_IN<-st_join(cets_sf, regions) 

# Assign datapoints to nearest region
Points_near<-st_nearest_feature(cets_sf, regions)
Points_IN$within <- Points_near

Points_IN <- cets_sf %>% mutate(
  within = as.integer(st_nearest_feature(cets_sf, regions))
  , area = if_else(is.na(within), '', regions$name[within])
) 

#Make the points a dataframe 
cets_spatial <- as.data.frame(Points_IN)

## SPATIAL COUNTS DATASETS ----
# Update columns
cets_spatial$ID <- cets_spatial$within
cets_spatial$name <- cets_spatial$area

## initialise empty list to store results of each group
spatial_list <- list()

for(group in unique_groups) {
  # Filter the dataframe for the current group
  spatial_df <- filter(cets_spatial, cets_spatial$group.r == group)
  
  # Count the occurrences for the current group
  group_s <- count(spatial_df, ID, name, year, month)
  group_s$ID <- as.factor(group_s$ID)
  
  # Create function to generate zero dataframes
  create_zero_df <- function(name, ID) {
    data.frame(year = 1992:2022, month = rep(1:12, each = 31), n = 0, name = name, ID = ID)
  }
  
  # List of areas and IDs
  areas <- list(
    c("Argyll", "1"), c("Clyde", "2"), c("Forth and Tay", "3"), c("Inner Minch", "10"),
    c("Inner Moray Firth", "9"), c("North Coast", "4"), c("North East", "5"), c("Orkney Islands", "6"),
    c("Outer Hebrides", "11"), c("Outer Moray Firth", "8"), c("Shetland Isles", "7"), c("Solway", "12")
  )
  
  # Create and merge dataframes in a loop
  for (area in areas) {
    t <- create_zero_df(area[1], area[2])
    group_s <- merge(group_s, t, all = TRUE)
  }
  
  # Handle duplicates
  group_s$duplicate <- duplicated(group_s[c("year", "month", "name", "ID")]) | 
    duplicated(group_s[c("year", "month", "name", "ID")], fromLast = TRUE)
  group_s <- group_s[!(group_s$duplicate & group_s$n == 0),]
  
  # Store the result in the list with the group name
  spatial_list[[group]] <- group_s
}
## SENSITIVITY DATASETS ----
# PER 1KM SHIFT
## Get 1km region shifted dataset 
one_regions <- st_read("C:/Users/Rachel Lennon/OneDrive - University of Glasgow/PhD/Projects/Trends/spatial/regions/Region shp/1km/1km.shp")
onepoints<-st_join(cets_sf, one_regions) 

onepoints_near <- onepoints %>% mutate(
  within = as.integer(st_nearest_feature(onepoints, one_regions))
  , area = if_else(is.na(within), '', one_regions$name[within])
) 

## Update column
onebal_r$name <- onebal_r$area

## initialise empty list to store results of each group
one_spatial_list <- list()

for(group in unique_groups) {
  # Filter the dataframe for the current group
  one_spatial_df <- filter(cets_spatial, cets_spatial$group.r == group)
  
  # Count the occurrences for the current group
  group_s_one <- count(one_spatial_df, ID, name, year, month)
  group_s_one$ID <- as.factor(group_s_one$ID)
  
  # Create function to generate zero dataframes
  create_zero_df <- function(name, ID) {
    data.frame(year = 1992:2022, month = rep(1:12, each = 31), n = 0, name = name, ID = ID)
  }
  
  # List of areas and IDs
  areas <- list(
    c("Argyll", "1"), c("Clyde", "2"), c("Forth and Tay", "3"), c("Inner Minch", "10"),
    c("Inner Moray Firth", "9"), c("North Coast", "4"), c("North East", "5"), c("Orkney Islands", "6"),
    c("Outer Hebrides", "11"), c("Outer Moray Firth", "8"), c("Shetland Isles", "7"), c("Solway", "12")
  )
  
  # Create and merge dataframes in a loop
  for (area in areas) {
    t <- create_zero_df(area[1], area[2])
    group_s_one <- merge(group_s_one, t, all = TRUE)
  }
  
  # Handle duplicates
  group_s_one$duplicate <- duplicated(group_s_one[c("year", "month", "name", "ID")]) | 
    duplicated(group_s_one[c("year", "month", "name", "ID")], fromLast = TRUE)
  group_s_one <- group_s_one[!(group_s_one$duplicate & group_s_one$n == 0),]
  
  # Store the result in the list with the group name
  one_spatial_list[[group]] <- group_s_one
}


# 5KM SHIFT dataset
## Get 5km region shifted dataset 
five_regions <- st_read("C:/Users/Rachel Lennon/OneDrive - University of Glasgow/PhD/Projects/Trends/spatial/regions/Region shp/5km/regions_5km_shift.shp")
fivepoints<-st_join(cets_sf, five_regions) 

fivepoints_near <- fivepoints %>% mutate(
  within = as.integer(st_nearest_feature(fivepoints, five_regions))
  , area = if_else(is.na(within), '', five_regions$name[within])
) 

## Update column
five_regions$name <- five_regions$area

## initialise empty list to store results of each group
five_spatial_list <- list()

for(group in unique_groups) {
  # Filter the dataframe for the current group
  five_spatial_df <- filter(cets_spatial, cets_spatial$group.r == group)
  
  # Count the occurrences for the current group
  group_s_five <- count(five_spatial_df, ID, name, year, month)
  group_s_five$ID <- as.factor(group_s_five$ID)
  
  # Create function to generate zero dataframes
  create_zero_df <- function(name, ID) {
    data.frame(year = 1992:2022, month = rep(1:12, each = 31), n = 0, name = name, ID = ID)
  }
  
  # List of areas and IDs
  areas <- list(
    c("Argyll", "1"), c("Clyde", "2"), c("Forth and Tay", "3"), c("Inner Minch", "10"),
    c("Inner Moray Firth", "9"), c("North Coast", "4"), c("North East", "5"), c("Orkney Islands", "6"),
    c("Outer Hebrides", "11"), c("Outer Moray Firth", "8"), c("Shetland Isles", "7"), c("Solway", "12")
  )
  
  # Create and merge dataframes in a loop
  for (area in areas) {
    t <- create_zero_df(area[1], area[2])
    group_s_five <- merge(group_s_five, t, all = TRUE)
  }
  
  # Handle duplicates
  group_s_five$duplicate <- duplicated(group_s_five[c("year", "month", "name", "ID")]) | 
    duplicated(group_s_five[c("year", "month", "name", "ID")], fromLast = TRUE)
  group_s_five <- group_s_five[!(group_s_five$duplicate & group_s_five$n == 0),]
  
  # Store the result in the list with the group name
  five_spatial_list[[group]] <- group_s_five
}

# 10KM SHIFT DATASET
## Get 10km region shifted dataset 
ten_regions <- st_read("C:/Users/Rachel Lennon/OneDrive - University of Glasgow/PhD/Projects/Trends/spatial/regions/Region shp/10km/regions_10km_shift.shp")
tenpoints<-st_join(cets_sf, ten_regions) 

tenpoints_near <- tenpoints %>% mutate(
  within = as.integer(st_nearest_feature(tenpoints, ten_regions))
  , area = if_else(is.na(within), '', ten_regions$name[within])
) 

## Update column
ten_regions$name <- ten_regions$area

## initialise empty list to store results of each group
ten_spatial_list <- list()

for(group in unique_groups) {
  # Filter the dataframe for the current group
  ten_spatial_df <- filter(cets_spatial, cets_spatial$group.r == group)
  
  # Count the occurrences for the current group
  group_s_ten <- count(ten_spatial_df, ID, name, year, month)
  group_s_ten$ID <- as.factor(group_s_ten$ID)
  
  # Create function to generate zero dataframes
  create_zero_df <- function(name, ID) {
    data.frame(year = 1992:2022, month = rep(1:12, each = 31), n = 0, name = name, ID = ID)
  }
  
  # List of areas and IDs
  areas <- list(
    c("Argyll", "1"), c("Clyde", "2"), c("Forth and Tay", "3"), c("Inner Minch", "10"),
    c("Inner Moray Firth", "9"), c("North Coast", "4"), c("North East", "5"), c("Orkney Islands", "6"),
    c("Outer Hebrides", "11"), c("Outer Moray Firth", "8"), c("Shetland Isles", "7"), c("Solway", "12")
  )
  
  # Create and merge dataframes in a loop
  for (area in areas) {
    t <- create_zero_df(area[1], area[2])
    group_s_ten <- merge(group_s_ten, t, all = TRUE)
  }
  
  # Handle duplicates
  group_s_ten$duplicate <- duplicated(group_s_ten[c("year", "month", "name", "ID")]) | 
    duplicated(group_s_ten[c("year", "month", "name", "ID")], fromLast = TRUE)
  group_s_ten <- group_s_ten[!(group_s_ten$duplicate & group_s_ten$n == 0),]
  
  # Store the result in the list with the group name
  ten_spatial_list[[group]] <- group_s_ten
}

## BALEEN WHALES ----
# RUN MODELS
## set up
## Make baleen spatial dataframe
bal_s <- spatial_list[["Baleen Whales"]]

## Create a neighboring matrix like above, specifying the names as the ones in the df
bal_nb <- poly2nb(bal_regions, row.names = bal_spat$ID)

## Attribute the same names to the name(nb)
names(bal_nb) <- attr(bal_nb, "region.id")

## Make the regions a factor
bal_s$ID <- as.factor(bal_s$ID)

## Models
  ## gam with mrf smoother to spatial element
  bal_spatialmod <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year, bs = "re") +
                        s(ID, bs="mrf", xt=list(nb=bal_nb)),
                      knots = knots, method = "REML", 
                      data = bal_s, family = poisson(link="log"))

  ## gam Without year as a random effect
  bal_spatialmod2 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=bal_nb)),
                       knots = knots, method = "REML", 
                       data = bal_s, family = poisson(link="log"))
#Compare
AIC(bal_spatialmod, bal_spatialmod2) #2 is better

  ## gam without the base variable of ID
  bal_spatialmod3 <- gam(n ~ s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=bal_nb)),
                       knots = knots, method = "REML", 
                       data = bal_s, family = poisson(link="log"))

#Compare
AIC(bal_spatialmod2, bal_spatialmod3) #2 is better, marginally

  ## GAM with ML instead of REML
  bal_spatialmod4 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=bal_nb)),
                       knots = knots, method = "ML", 
                       data = bal_s, family = poisson(link="log"))

#Compare
AIC(bal_spatialmod2, bal_spatialmod4) # They are the same

  ## GAM without the mrf smoother in ID
  bal_spatialmod5 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year),
                       knots = knots, method = "REML", 
                       data = bal_s, family = poisson(link="log"))

#Compare
AIC(bal_spatialmod2, bal_spatialmod5) # they are the same, no interdependence so having no effect?
#Compare fit
gam.check(bal_spatialmod2) #looks better
gam.check(bal_spatialmod5)

  ## GAM without year 
  bal_spatialmod7 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +
                         s(ID, bs="mrf", xt=list(nb=nb)),
                       knots = knots, method = "REML", 
                       data = bal_s, family = poisson(link="log"))

#Compare
AIC(bal_spatialmod2, bal_spatialmod7) #Better with year

  ## GAM with year as a random effect?
bal_spatialmod8 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) + s(year, bs = "re") +
                         s(ID, bs="mrf", xt=list(nb=nb)),
                       knots = knots, method = "REML", 
                       data = bal_s, family = poisson(link="log"))

#Compare
AIC(bal_spatialmod2, bal_spatialmod8) #Better as a fixed effect

  ## GAM with interaction of year & month
bal_spatialmod9 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) + s(month, bs = "cc", k = 12, by = as.factor(year)) + as.factor(year) +
                         s(ID, bs="mrf", xt=list(nb=nb)),
                       knots = knots, method = "REML", 
                       data = bal_s, family = poisson(link="log"))

#Compare
AIC(bal_spatialmod2, bal_spatialmod9) #without is better indicating no change over time.
##Bal_spatialmod2 is best, look at diagnostics

## Diagnostic checks for bal_spatialmod2

## summary and visualise
gam.check(bal_spatialmod2) 
check_overdispersion(bal_spatialmod2) # = 0.11, ok as poisson
## good fit. 

## read results
summary(bal_spatialmod2) ## R2 = 0.11
plot(bal_spatialmod2)

##
##
##

# SENSITIVITY CHECKS
## per species
bal_species<- count(bal_spat, area, year, month, spec.r)

#fin
fin_spat <- filter(bal_species, bal_species$spec.r == "Fin whale")
fin_spat$name <- fin_spat$area
fin_spat$name <- as.factor(fin_spat$name)

fin_spat_model <- glm(n ~ name*year*month, data = fin_spat, family = poisson(link = "log"))
summary(fin_spat_model)
## no trend, not enough data alone

#humpback
hump_spat <- filter(bal_species, bal_species$spec.r == "Humpback whale")
hump_spat$name <- hump_spat$area
hump_spat$name <- as.factor(hump_spat$name)

hump_spat_model <- glm(n ~ name*year*month, data = hump_spat, family = poisson(link = "log"))
summary(hump_spat_model)
## no trend, not enough data alone

#minke
minke_spat <- filter(bal_species, bal_species$spec.r == "Minke whale")
minke_spat$name <- minke_spat$area
minke_spat$name <- as.factor(minke_spat$name)

minke_spat_model <- gam(n ~ name*year*month, data = minke_spat, family = poisson(link = "log"))
summary(minke_spat_model)
## no trend, not enough data alone
## Need to combine all to see a trend.

# PER 1KM SHIFT
## Get 1km region shifted for baleens 
obal_s <- one_spatial_list[["Baleen Whales"]]

# Check trends are similar with basic model 
obal_s$name <- as.factor(obal_s$name)

  ## GAM with 1km shift model
  obal_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                       knots = knots, method = "REML", 
                       data = obal_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(obal_spatialmod) #ok
summary(obal_spatialmod) ## R2 = 0.07
plot(obal_spatialmod) ## Looks the same, trend wise.
check_overdispersion(obal_spatialmod) # = 0.08
## OK.

## PER 5km SHIFT
## Get 1km region shifted for baleens 
fivebal_s <- five_spatial_list[["Baleen Whales"]]

# Check trends are similar with basic model 
fivebal_s$name <- as.factor(fivebal_s$name)

## GAM with 1km shift model
fivebal_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                       knots = knots, method = "REML", 
                       data = fivebal_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(fivebal_spatialmod) #ok
summary(fivebal_spatialmod) ## R2 = 0.07
plot(fivebal_spatialmod) ## Looks the same, trend wise.
check_overdispersion(fivebal_spatialmod) # = 0.08, ok as poisson.
## OK.


## PER 10km SHIFT
## Get 1km region shifted for baleens 
tenbal_s <- ten_spatial_list[["Baleen Whales"]]

# Check trends are similar with basic model 
tenbal_s$name <- as.factor(tenbal_s$name)

## GAM with 1km shift model
tenbal_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                          knots = knots, method = "REML", 
                          data = tenbal_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(tenbal_spatialmod) #ok
summary(tenbal_spatialmod) ## R2 = 0.07
plot(tenbal_spatialmod) ## Looks the same, trend wise.
check_overdispersion(tenbal_spatialmod) # = 0.08, ok as poisson.
## Trends are the same, signif has been lost in some. 

## COMMON DOLPHINS ----
# AUTOCORRELATION TEST 
## Make common specific spatial dataframe
com_spat <- filter(cets_spatial, group.r == "Common Dolphin")

## make a common specific regions file
com_regions <- regions

## Generate a list of neighboring polygons 
com_nb <- poly2nb(regions, queen=TRUE) 

## Assign equal weight to neighbours
com_rw <- nb2listw(com_nb, style="W", zero.policy=TRUE)

## Calculate the stranding value for each 
com_num <- count(com_spat, area)
com_regions$value <- com_num$n

## Moran's I: Monte Carlo
com_MC<- moran.mc(com_regions$value, com_rw, nsim = 999, alternative = "greater")
com_MC
plot(com_MC) ## Accept hypothesis of random distribution, i.e. not spatial autocorrelation

# RUN MODELS
## set up
## Make common spatial dataframe
com_s <- spatial_list[["Common Dolphin"]]

## Create a neighboring matrix like above, specifying the names as the ones in the df
com_nb <- poly2nb(com_regions, row.names = com_spat$ID)

## Attribute the same names to the name(nb)
names(com_nb) <- attr(com_nb, "region.id")

## Make the regions a factor
com_s$ID <- as.factor(com_s$ID)

## Models
  ## gam with mrf smoother to spatial element
  com_spatialmod <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year, bs = "re") +
                        s(ID, bs="mrf", xt=list(nb=com_nb)),
                      knots = knots, method = "REML", 
                      data = com_s, family = poisson(link="log"))

  ## gam Without year as a random effect
  com_spatialmod2 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=com_nb)),
                       knots = knots, method = "REML", 
                       data = com_s, family = poisson(link="log"))

# Compare
AIC(com_spatialmod, com_spatialmod2) #1st is better, but very negative which is a concern?

  ## GAM without the base variable of ID
com_spatialmod3 <- gam(n ~ s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=com_nb)),
                       knots = knots, method = "REML", 
                       data = com_s, family = poisson(link="log"))

# Compare
AIC(com_spatialmod2, com_spatialmod3) #2 is better, marginally

  ## GAM with ML instead of REML
com_spatialmod4 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=com_nb)),
                       knots = knots, method = "ML", 
                       data = com_s, family = poisson(link="log"))

# Compare
AIC(com_spatialmod2, com_spatialmod4) # They are the same

  ## GAM without the mrf smoother in ID
com_spatialmod5 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year),
                       knots = knots, method = "REML", 
                       data = com_s, family = poisson(link="log"))

# Compare
AIC(com_spatialmod2, com_spatialmod5) # they are the same...no interdependence so having no effect?
gam.check(com_spatialmod5) 
gam.check(com_spatialmod2) #looks marginally better

## Run the diagnostics
check_overdispersion(com_spatialmod2) # 0.2

## Read results
summary(com_spatialmod2) ## R2 = 0.29
plot(com_spatialmod2)

##
##
##

# SENSITIVITY CHECKS

# PER 1KM SHIFT
## Get 1km region shifted for baleens 
ocom_s <- one_spatial_list[["Common Dolphin"]]

# Check trends are similar with basic model 
ocom_s$name <- as.factor(ocom_s$name)

## GAM with 1km shift model
ocom_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                       knots = knots, method = "REML", 
                       data = ocom_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(ocom_spatialmod) #ok
summary(ocom_spatialmod) ## R2 = 0.07
plot(ocom_spatialmod) ## Looks the same, trend wise.
check_overdispersion(ocom_spatialmod) # = 0.08
## OK.

## PER 5km SHIFT
## Get 5km region shifted for commons 
fivecom_s <- five_spatial_list[["Common Dolphins"]]

# Check trends are similar with basic model 
fivecom_s$name <- as.factor(fivecom_s$name)

## GAM with 1km shift model
fivecom_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                          knots = knots, method = "REML", 
                          data = fivecom_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(fivecom_spatialmod) #ok
summary(fivecom_spatialmod) ## R2 = 0.07
plot(fivecom_spatialmod) ## Looks the same, trend wise.
check_overdispersion(fivecom_spatialmod) # = 0.08, ok as poisson.
## OK.


## PER 10km SHIFT
## Get 10km region shifted for commons 
tencom_s <- ten_spatial_list[["Common Dolphin"]]

# Check trends are similar with basic model 
tencom_s$name <- as.factor(tencom_s$name)

## GAM with 1km shift model
tencom_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                         knots = knots, method = "REML", 
                         data = tencom_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(tencom_spatialmod) #ok
summary(tencom_spatialmod) ## R2 = 0.07
plot(tencom_spatialmod) ## Looks the same, trend wise.
check_overdispersion(tencom_spatialmod) # = 0.08, ok as poisson.
## Trends are the same, signif has been lost in some. 

## DEEP DIVERS ----
## Make deep diver spatial dataframe
deep_s <- spatial_list[["Deep Divers"]]

## Create a neighboring matrix like above, specifying the names as the ones in the df
deep_nb <- poly2nb(regions, row.names = deep_s$ID)

## Attribute the same names to the name(nb)
names(deep_nb) <- attr(deep_nb, "region.id")

## Make the regions a factor
deep_s$ID <- as.factor(deep_s$ID)
## Remove 2018 UME
deep_s <- filter(deep_s, year != "2018")

## Models
  ## gam with mrf smoother to spatial element
  deep_spatialmod <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year, bs = "re") +
                         s(ID, bs="mrf", xt=list(nb=deep_nb)),
                       knots = knots, method = "REML", 
                       data = deep_s, family = poisson(link="log"))

  ## gam as a ziP
  deep_spatialmod1 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year, bs = "re") +
                          s(ID, bs="mrf", xt=list(nb=deep_nb)),
                        knots = knots, method = "REML", 
                        data = deep_s, family = ziP())

#compare
AIC(deep_spatialmod, deep_spatialmod1) ## ziP is better. 

  ## gam Without year as a random effect
  deep_spatialmod3 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                          s(ID, bs="mrf", xt=list(nb=deep_nb)),
                        knots = knots, method = "REML", 
                        data = deep_s, family = ziP())
#compare
AIC(deep_spatialmod, deep_spatialmod3) ##mod3 better


  ## gam With ML instead of REML
  deep_spatialmod4 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                          s(ID, bs="mrf", xt=list(nb=deep_nb)),
                        knots = knots, method = "ML", 
                        data = deep_s, family = ziP())

#compare
AIC(deep_spatialmod3, deep_spatialmod4) # reml better

  ## gam without the mrf smoother in ID
deep_spatialmod5 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year),
                        knots = knots, method = "REML", 
                        data = deep_s, family = ziP())

#compare
AIC(deep_spatialmod3, deep_spatialmod5) # they are the same...no interdependence so having no effect?

gam.check(deep_spatialmod5)
gam.check(deep_spatialmod3) #somewhat better 

## Diagnostics
check_overdispersion(deep_spatialmod3) # 0.19

## Read results
summary(deep_spatialmod3) ## R2 = 0.506
plot(deep_spatialmod3)

##
##
##

# SENSITIVITY CHECKS

# PER 1KM SHIFT
## Get 1km region shifted for deep divers 
odeep_s <- one_spatial_list[["Deep Divers"]]

# Check trends are similar with basic model 
odeep_s$name <- as.factor(odeep_s$name)

## GAM with 1km shift model
odeep_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                       knots = knots, method = "REML", 
                       data = odeep_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(odeep_spatialmod) #ok
summary(odeep_spatialmod) ## R2 = 0.07
plot(odeep_spatialmod) ## Looks the same, trend wise.
check_overdispersion(odeep_spatialmod) # = 0.08
## OK.

## PER 5km SHIFT
## Get 5km region shifted for deep divers 
fivedeep_s <- five_spatial_list[["Deep Divers"]]

# Check trends are similar with basic model 
fivedeep_s$name <- as.factor(fivedeep_s$name)

## GAM with 1km shift model
fivedeep_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                          knots = knots, method = "REML", 
                          data = fivedeep_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(fivedeep_spatialmod) #ok
summary(fivedeep_spatialmod) ## R2 = 0.07
plot(fivedeep_spatialmod) ## Looks the same, trend wise.
check_overdispersion(fivedeep_spatialmod) # = 0.08, ok as poisson.
## OK.


## PER 10km SHIFT
## Get 10km region shifted for deep divers 
tendeep_s <- ten_spatial_list[["Deep divers"]]

# Check trends are similar with basic model 
tendeep_s$name <- as.factor(tendeep_s$name)

## GAM with 1km shift model
tendeep_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                         knots = knots, method = "REML", 
                         data = tendeep_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(tendeep_spatialmod) #ok
summary(tendeep_spatialmod) ## R2 = 0.07
plot(tendeep_spatialmod) ## Looks the same, trend wise.
check_overdispersion(tendeep_spatialmod) # = 0.08, ok as poisson.
## Trends are the same, signif has been lost in some. 

## HARBOUR PORPOISE ----
## Make harbour porpoise spatial dataframe
hp_s <- spatial_list[["Harbour porpoise"]]

## Create a neighboring matrix like above, specifying the names as the ones in the df
hp_nb <- poly2nb(regions, row.names = hp_s$ID)

## Attribute the same names to the name(nb)
names(hp_nb) <- attr(hp_nb, "region.id")

## Make the regions a factor
hp_s$ID <- as.factor(hp_s$ID)

## Models
  ## gam with mrf smoother to spatial element
  hp_spatialmod <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year, bs = "re") +
                       s(ID, bs="mrf", xt=list(nb=hp_nb)),
                     knots = knots, method = "REML", 
                     data = hp_s, family = poisson(link="log"))

  ## gam without year as a random effect
  hp_spatialmod2 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                        s(ID, bs="mrf", xt=list(nb=hp_nb)),
                      knots = knots, method = "REML", 
                      data = hp_s, family = poisson(link="log"))

#compare
AIC(hp_spatialmod, hp_spatialmod2) #2 is better

  ## gam without the base variable of ID
  hp_spatialmod3 <- gam(n ~ s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                        s(ID, bs="mrf", xt=list(nb=hp_nb)),
                      knots = knots, method = "REML", 
                      data = hp_s, family = poisson(link="log"))

#compare
AIC(hp_spatialmod2, hp_spatialmod3) #2 is better, marginally

  ## gam with ML instead of REML
  hp_spatialmod4 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                        s(ID, bs="mrf", xt=list(nb=hp_nb)),
                      knots = knots, method = "ML", 
                      data = hp_s, family = poisson(link="log"))

#compare
AIC(hp_spatialmod2, hp_spatialmod4) # They are the same

  ## gam without the mrf smoother in ID
hp_spatialmod5 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year),
                      knots = knots, method = "REML", 
                      data = hp_s, family = poisson(link="log"))

#compare
AIC(hp_spatialmod2, hp_spatialmod5) # they are the same...no interdependence so having no effect?
# look at fit between 
gam.check(hp_spatialmod5)
gam.check(hp_spatialmod2) # alittle bit better


  ## gam without year 
hp_spatialmod7 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +
                        s(ID, bs="mrf", xt=list(nb=hp_nb)),
                      knots = knots, method = "REML", 
                      data = hp_s, family = poisson(link="log"))

#compare
AIC(hp_spatialmod2, hp_spatialmod7) #2 is better

  ## gam with a negative binomial distribution 
hp_spatialmod8 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                        s(ID, bs="mrf", xt=list(nb=hp_nb)),
                      knots = knots, method = "REML", 
                      data = hp_s, family = nb(theta=NULL, link="log"))

#compare
AIC(hp_spatialmod2, hp_spatialmod8)
#spatial mod 2 is better. 

## Diagnostic checks 
gam.check(hp_spatialmod2) #not great
check_overdispersion(hp_spatialmod2) # = 0.45

## read results
summary(hp_spatialmod2) 
plot(hp_spatialmod8)

##
##
## 

## SENSITIVITY CHECKS

# PER 1KM SHIFT
## Get 1km region shifted for harbour porpoise
oporp_s <- one_spatial_list[["Harbour porpoise"]]

# Check trends are similar with basic model 
oporp_s$name <- as.factor(oporp_s$name)

## GAM with 1km shift model
oporp_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                       knots = knots, method = "REML", 
                       data = oporp_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(oporp_spatialmod) #ok
summary(oporp_spatialmod) ## R2 = 0.07
plot(oporp_spatialmod) ## Looks the same, trend wise.
check_overdispersion(oporp_spatialmod) # = 0.08
## OK.

## PER 5km SHIFT
## Get 5km region shifted for harbour porpoise 
fiveporp_s <- five_spatial_list[["Harbour porpoise"]]

# Check trends are similar with basic model 
fiveporp_s$name <- as.factor(fiveporp_s$name)

## GAM with 1km shift model
fiveporp_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                          knots = knots, method = "REML", 
                          data = fiveporp_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(fiveporp_spatialmod) #ok
summary(fiveporp_spatialmod) ## R2 = 0.07
plot(fiveporp_spatialmod) ## Looks the same, trend wise.
check_overdispersion(fiveporp_spatialmod) # = 0.08, ok as poisson.
## OK.


## PER 10km SHIFT
## Get 10km region shifted for harbour porpoise 
tenporp_s <- ten_spatial_list[["Harbour porpoise"]]

# Check trends are similar with basic model 
tenporp_s$name <- as.factor(tenporp_s$name)

## GAM with 1km shift model
tenporp_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                         knots = knots, method = "REML", 
                         data = tenporp_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(tenporp_spatialmod) #ok
summary(tenporp_spatialmod) ## R2 = 0.07
plot(tenporp_spatialmod) ## Looks the same, trend wise.
check_overdispersion(tenporp_spatialmod) # = 0.08, ok as poisson.
## Trends are the same, signif has been lost in some. 

## PELAGIC DOLPHINS ----
## Make pelaic spatial dataframe
pel_s <- spatial_list[["Pelagic Dolphins"]]

## Create a neighboring matrix like above, specifying the names as the ones in the df
pel_nb <- poly2nb(regions, row.names = pel_s$ID)

## Attribute the same names to the name(nb)
names(pel_nb) <- attr(pel_nb, "region.id")

## Make the regions a factor
pel_s$ID <- as.factor(pel_s$ID)

## Models

  ## gam with mrf smoother to spatial element
  pel_spatialmod <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year, bs = "re") +
                        s(ID, bs="mrf", xt=list(nb=pel_nb)),
                      knots = knots, method = "REML", 
                      data = pel_s, family = poisson(link="log"))

  ## gam Without year as a random effect
  pel_spatialmod2 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=pel_nb)),
                       knots = knots, method = "REML", 
                       data = pel_s, family = poisson(link="log"))

#compare
AIC(pel_spatialmod, pel_spatialmod2) #the same so go without

  ## GAM with ML instead of REML
  pel_spatialmod3 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=pel_nb)),
                       knots = knots, method = "ML", 
                       data = pel_s, family = poisson(link="log"))

#compare
AIC(pel_spatialmod2, pel_spatialmod3) # They are the same

  ## GAM as a nb
  pel_spatialmod4 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year) +
                         s(ID, bs="mrf", xt=list(nb=pel_nb)),
                       knots = knots, method = "REML", 
                       data = pel_s, family = nb(theta=NULL, link="log"))

#compare
AIC(pel_spatialmod2, pel_spatialmod4)
check_overdispersion(pel_spatialmod2) # = 0.14
# data is not overdispersed.

  # GAM without the mrf smoother in ID
pel_spatialmod5 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +  s(year),
                       knots = knots, method = "REML", 
                       data = pel_s, family = poisson(link="log"))

#compare
AIC(pel_spatialmod2, pel_spatialmod5) # they are the same...no interdependence so having no effect?
gam.check(pel_spatialmod5) 
gam.check(pel_spatialmod2) ## this is a better fit.

  ## GAM without year 
  pel_spatialmod6 <- gam(n ~ ID + s(month, bs = "cc", k = 12, by = ID) +
                         s(ID, bs="mrf", xt=list(nb=pel_nb)),
                       knots = knots, method = "REML", 
                       data = pel_s, family = poisson(link="log"))

#compare
AIC(pel_spatialmod2, pel_spatialmod6) #2 is better

## Diagnostic checks 
## pelspatialmod2
gam.check(pel_spatialmod2) #funky looking
check_overdispersion(pel_spatialmod2) # = 0.14

# read results
summary(pel_spatialmod2) #0.17
plot(pel_spatialmod2)

##
##
##

# SENSITIVITY CHECKS

# PER SPECIES
## make dataframe list of all species
pel_species<- count(cets_spatial, name, year, month, group.r, spec.r)
pel_species <- filter(pel_species, group.r == "Pelagic Dolphins")
pel_groups <- unique(pel_species$spec.r)

## initialise empty list 
pelagic_list <- list()

for(group in pel_groups) {
  # Filter the dataframe for the current group
  group_df <- filter(pel_species, pel_species$spec.r == group)
  
  # Create function to generate zero dataframes
  create_zero_df <- function(name) {
    data.frame(year = 1992:2022, month = rep(1:12, each = 31), n = 0, name = name)
  }
  
  # List of areas and IDs
  areas <- list(
    c("Argyll"), c("Clyde"), c("Forth and Tay"), c("Inner Minch"),
    c("Inner Moray Firth"), c("North Coast"), c("North East"), c("Orkney Islands"),
    c("Outer Hebrides"), c("Outer Moray Firth"), c("Shetland Isles"), c("Solway")
  )
  
  # Create and merge dataframes in a loop
  for (area in areas) {
    t <- create_zero_df(area[1])
    group_df <- merge(group_df, t, all = TRUE)
  }
  
  # Handle duplicates
  group_df$duplicate <- duplicated(group_df[c("year", "month", "name")]) | 
    duplicated(group_df[c("year", "month", "name")], fromLast = TRUE)
  group_df <- group_df[!(group_df$duplicate & group_df$n == 0),]
  
  group_df$name <- as.factor(group_df$name)
  
  # Store the result in the list with the group name
  pelagic_list[[group]] <- group_df
}

#lfpw
## get df for lfpw
lfpw_spat <- pelagic_list[["Long-finned pilot whale"]]

## run simple model to see if trend holds
lfpw_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year), 
                       knots = knots, method = "REML", 
                       data = lfpw_spat, family = poisson(link="log"))

## results
summary(lfpw_spatialmod)
layout(matrix(1:12, ncol = 4))
plot(lfpw_spatialmod)
## no trend 

#striped dolphin 
striped_spat <- pelagic_list[["Striped dolphin"]]

## run simple model to see if trend holds
sd_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year), 
                     knots = knots, method = "REML", 
                     data = striped_spat, family = poisson(link="log"))

## results
summary(sd_spatialmod)
plot(sd_spatialmod)
## Slight up in winter for OH

#risso's dolphin
rissos_spat <- pelagic_list[["Risso's dolphin"]]

## run simple model to see if trend holds
rissos_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year), 
                         knots = knots, method = "REML", 
                         data = rissos_spat, family = poisson(link="log"))

## results
summary(rissos_spatialmod)
plot(rissos_spatialmod)
#Summer bump seen in Argyll.

#white-beaked
white_spat <- pelagic_list[["White-beaked dolphin"]]

## run simple model to see if trend holds
white_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year), 
                        knots = knots, method = "REML", 
                        data = white_spat, family = poisson(link="log"))

## Results
summary(white_spatialmod)
plot(white_spatialmod)
## Inner minch, north coast, north east, orkney and outerhebrides matches except
# don't have the bump in winter in the OH. i.e. not vulnerable to fishing?

#awsd
awsd_spat <- pelagic_list[["Atlantic white-sided dolphin"]]

## run simple model to see if trend holds
awsd_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year), 
                       knots = knots, method = "REML", 
                       data = awsd_spat, family = poisson(link="log"))

## results
summary(awsd_spatialmod)
plot(awsd_spatialmod)
## Some differences here. Slight winter bump in OH


#bnd
bnd_spat <- pelagic_list[["Bottlenose dolphin"]]

## run simple model to see if trend holds
bnd_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year), 
                      knots = knots, method = "REML", 
                      data = bnd_spat, family = poisson(link="log"))

## results
summary(bnd_spatialmod)
plot(bnd_spatialmod)
## same as group trend for those that are signif. 

#killer
killer_spat <- pelagic_list[["Killer whale"]]

## run simple model to see if trend holds
killer_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year), 
                         knots = knots, method = "REML", 
                         data = killer_spat, family = poisson(link="log"))

## results
summary(killer_spatialmod)
plot(killer_spatialmod)
layout(matrix(1))
# no trends

# PER 1KM SHIFT
## Get 1km region shifted for pelagic dolphins
opela_s <- one_spatial_list[["Pelagic Dolphins"]]

# Check trends are similar with basic model 
opela_s$name <- as.factor(opela_s$name)

## GAM with 1km shift model
opela_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                       knots = knots, method = "REML", 
                       data = opela_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(opela_spatialmod) #ok
summary(opela_spatialmod) ## R2 = 0.07
plot(opela_spatialmod) ## Looks the same, trend wise.
check_overdispersion(opela_spatialmod) # = 0.08
## OK.

## PER 5km SHIFT
## Get 5km region shifted for pelagic 
fivepela_s <- five_spatial_list[["Pelagic Dolphins"]]

# Check trends are similar with basic model 
fivepela_s$name <- as.factor(fivepela_s$name)

## GAM with 1km shift model
fivepela_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                          knots = knots, method = "REML", 
                          data = fivepela_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(fivepela_spatialmod) #ok
summary(fivepela_spatialmod) ## R2 = 0.07
plot(fivepela_spatialmod) ## Looks the same, trend wise.
check_overdispersion(fivepela_spatialmod) # = 0.08, ok as poisson.
## OK.


## PER 10km SHIFT
## Get 10km region shifted for pelagic dolphins 
tenpela_s <- ten_spatial_list[["Pelagic Dolphin"]]

# Check trends are similar with basic model 
tenpela_s$name <- as.factor(tenpela_s$name)

## GAM with 1km shift model
tenpela_spatialmod <- gam(n ~ name + s(month, bs = "cc", k = 12, by = name) +  s(year),
                         knots = knots, method = "REML", 
                         data = tenpela_s, family = poisson(link="log"))

## Diagnostic checks 
## Visualise
gam.check(tenpela_spatialmod) #ok
summary(tenpela_spatialmod) ## R2 = 0.07
plot(tenpela_spatialmod) ## Looks the same, trend wise.
check_overdispersion(tenpela_spatialmod) # = 0.08, ok as poisson.
## Trends are the same, signif has been lost in some. 

