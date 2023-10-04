
#Load packages

library(randomForest)
library(feather)
library(tidyverse)
library(tictoc)
library(terra)
library(VIM)
library(groupdata2)
library(lsa)
#library(missForest)



#Specify Drive Path
drive_path <- "//worldpop.files.soton.ac.uk/worldpop/Projects/WP517763_GRID3/"
input_path <-  paste0(drive_path, "Working/GHA/Ortis/Output/")
raster_path <- paste0(drive_path, "Working/GHA/Ortis/Other_covariates/")
output_path <- paste0(drive_path, "Working/GHA/Ortis/Output/")
output_path1 <- paste0(drive_path, "Working/GHA/Ortis/Output/Predicted Population/")
output_path2 <- paste0(drive_path, "Working/GHA/Ortis/Output/Posterior Predictions/")

#Load Ghana population dataset

GHA_df <- read.csv(paste0(input_path, "/GHA_Data_df.csv"))

head(GHA_df[,1:5]) # only showing first five columns


#Define response variable as log of pop_density

GHA_df <-GHA_df %>% 
  mutate(pop_density = log(Pop_2021/Area)) %>% 
  select(-X) 

#Covariates selection
covs <- GHA_df %>% 
  select(starts_with("x"))

# Calculate mean and standard deviation of covariates
cov_stats <- data.frame(Covariate = colnames(covs),
                        Mean = apply(covs, 2, mean, na.rm = TRUE),
                        Std_Dev = apply(covs, 2, sd, na.rm = TRUE))

#Scaling function to scale covariates
stdize <- function(x)
{ stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
return(stdz) }

#apply scaling function
covs <- apply(covs, 2, stdize)   #z-score

head(covs[,1:2]) # only showing first two columns

#Select pop_density and cbind covs
GHA_df2 <-GHA_df %>% 
select(Pop_2021, pop_density, Dist_ID, Area) %>% 
cbind(covs) 


# Fit model to all the training data -------------------------------------------------

set.seed(4567)

# Search for best hyperparameters tunning

tuneRF(x = covs, y = GHA_df2$pop_density, na.action = na.omit,
       plot = T, trace = T, importance=TRUE,  sampsize=length(GHA_df2), replace=TRUE) 


#Fit model
model1 <- randomForest(x = covs, y = GHA_df2$pop_density, mtry = 16, na.action = na.omit, 
                      plot = T, trace = T, importance=TRUE, sampsize=length(GHA_df2), replace=TRUE) 


model1


#All data predictions
model1_predictions <- model1$predicted %>% as_tibble()

#cbind predicted data to original data
model1_predictions <- model1_predictions %>% 
  cbind(GHA_df2$pop_density) %>% 
  mutate(observed = exp(GHA_df2$pop_density), predicted = exp(value),
         residual = predicted - observed,
         model = "RandomForest")

write.csv(model1_predictions, paste0(output_path, "RandomForest model results.csv"))

# compute goodness-of-fit metrics
model1_predictions %>%
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            rmse = sqrt(mse),
            corr = cor(predicted, observed))



# Model plots -------------------------------------------------------------
#Plot Population Density Estimate

ggplot(model1_predictions) +
  geom_point(aes(x = observed, y = predicted), fill = 'grey50', color = 'grey70', shape = 21) +
  geom_smooth(aes(x = observed, y = predicted), method = lm) + 
  #geom_abline(slope = 1, intercept = 0, color = 'orange', linewidth = 1) +
  theme_minimal() +
  labs(title = '', x = 'Observed population density', y = 'Predicted population density')


#Histogram Plot

model1_predictions %>% 
  pivot_longer(cols = c(predicted, observed), names_to = "Density", 
               values_to = "predicted_density") %>%
 filter(predicted_density <0.2) %>% 
  ggplot() +aes(predicted_density, fill = Density) +geom_histogram(alpha = 0.5, bins = 100)+
  labs(title = 'Histogram plot', x='Predicted Density', y='Frequency')+ 
  theme_bw()+
  scale_fill_discrete(name="Density",
                      breaks=c("predicted", "observed"),
                      labels=c("Predicted Population Density", "Observed Population Density"))
#Density Plot

model1_predictions %>% 
  pivot_longer(cols = c(predicted, observed), names_to = "Density", 
               values_to = "predicted_density") %>%
  filter(predicted_density <0.2) %>% 
  ggplot() +aes(predicted_density, fill = Density) +geom_density(alpha = 0.5)+
  labs(title = 'Density plot', x='Predicted Population Density', y='Frequency')+ 
  theme_bw()+
  scale_fill_discrete(name="Density",
                      breaks=c("predicted", "observed"),
                      labels=c("Predicted Population Density", "Observed Population Density"))

#Variable importance
#Read variable names
original_names <- read.csv(paste0(output_path, "var_names.csv"))

varImpPlot(model1)

var_importance <- importance(model1)


var_importance_df <- var_importance %>% 
  as_tibble() %>% 
mutate(variable = row.names(var_importance)) %>% 
  rename(inc_prop = "%IncMSE")

#Join var_importance-df to original names
var_importance_df <- var_importance_df %>% 
  inner_join(original_names, by = c("variable" = "var_names2"))%>% 
  mutate(Model = "RandomForest")

write.csv(var_importance_df, paste0(output_path, "RandomForest_var_importance.csv"))

#plot variable importance 
ggplot(var_importance_df, aes(x = reorder(Original.Name, inc_prop), y = inc_prop, fill = Original.Name)) +
  geom_bar(stat = "identity")+
  geom_text(aes(label = round(inc_prop,3)), hjust=-0.2, size=5) + # add y values as labels to the bars
  coord_flip() +
  theme_bw()+
  labs(x = "Variables", y = "Variable Importance(%)")+
  scale_fill_manual(values = c(rep("#8c2981", 27))) + # use only one color for the fill
  theme(legend.position="none",  axis.text.y = element_text(size = 14))


# Perform Out-of-Sample Cross Validation ----------------------------------

#Cross Validation using Training and Test data

#Create training dataset
train <- GHA_df %>% 
  sample_frac(.70)

#Create test set
test <- anti_join(GHA_df, train, by = "Dist_ID")

#train covariates
covs_train <- train %>% 
  select(starts_with("x"))

# Calculate mean and standard deviation of covariates
covs_train_stats <- data.frame(Covariate = colnames(covs_train),
                               Mean = apply(covs_train, 2, mean, na.rm = TRUE),
                               Std_Dev = apply(covs_train, 2, sd, na.rm = TRUE))

#apply scaling function
covs_train <- apply(covs_train, 2, stdize) %>%    #z-score
  as_tibble()

head(covs_train[,1:2]) # only showing first two columns

#fit model to train dataset

model2 <- randomForest(x = covs_train, y = train$pop_density, mtry = 16, na.action = na.omit,
                       plot = T, trace = T, importance=TRUE,  sampsize=length(train), replace=TRUE) 


model2

model2_predictions <- model2$predicted %>% as_tibble()

#cbind predicted data to train data
model2_predictions <- model2_predictions %>% 
  cbind(train$pop_density) %>% 
  mutate(observed = exp(train$pop_density), predicted = exp(value),
         residual = predicted - observed)


# compute goodness-of-fit metrics
model2_predictions %>%
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            rmse = sqrt(mse),
            corr = cor(predicted, observed))



#Make predictions on the test data
covs_test <- test %>% 
  select(starts_with("x"))

#Scale covariates with train data mean and sd
for (var in names(covs_test)) {
  var_mean <- covs_train_stats$Mean[covs_train_stats$Covariate == var]
  var_sd <- covs_train_stats$Std_Dev[covs_train_stats$Covariate == var]
  covs_test[[var]] <- (covs_test[[var]] - var_mean) / var_sd
}


head(covs_test[,1:2]) 

#Make predictions on the test data

test_predicted <- predict(model2, newdata = covs_test)

#cbind to test data
test <- test %>% 
  select(pop_density) %>% 
  cbind(test_predicted)%>% 
  mutate(observed = exp(pop_density), predicted = exp(test_predicted),
         residual = predicted - observed)


# compute goodness-of-fit metrics
test %>% 
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            rmse = sqrt(mse),
            corr = cor(predicted, observed))


rm(covs_test, covs_train, covs, train, test)

# Weighting Layer Analysis (Predictions) ------------------------------------------------

settled_df <- read_feather(paste0(input_path, "settled_df.feather"))
#settled_df <- read.csv(paste0(input_path, "settled_df.csv"))

#Based on simulation study we will impute missing values using their mean value
#Select covariates
covs1 <- settled_df %>% 
  select(starts_with("x"))

#apply mean imputation
covs1 <- covs1 %>%
  mutate_all(~as.numeric(.x)) %>% 
  mutate_all(~replace_na(.x, mean(.x, na.rm = TRUE)))

head(covs1)

# Scale covariates using means and standard deviations from covs_stat

for (var in names(covs1)) {
  var_mean <- cov_stats$Mean[cov_stats$Covariate == var]
  var_sd <- cov_stats$Std_Dev[cov_stats$Covariate == var]
  covs1[[var]] <- (covs1[[var]] - var_mean) / var_sd
}

# Viewing the scaled dataframe
head(covs1)

#check for NA values
any(is.na(covs1))


#Make predictions on scaled covariates
predicted <- predict(model1, newdata = covs1)



# Disagregate Population Totals -------------------------------------------

#cbind predictions to settled_df

population_predictions <- settled_df %>% 
  select(Grid_ID, Dist_ID) %>% 
  cbind(predicted) %>% 
  mutate(predicted_exp = exp(predicted)) # back-transform predictions to natural scale
 
#Sum exponentiated predictions for all pixels in a given district
population_predictions <- population_predictions %>% 
  group_by(Dist_ID) %>% 
  mutate(predicted_exp_sum = sum(predicted_exp)) %>% 
  ungroup()

#Select Population totals from GHA_df and merge with predicted population df
district_total_pop <- GHA_df %>% 
  select(Pop_2021, Dist_ID)%>% 
  rename(district_pop = Pop_2021)

#Join to population
population_predictions <-population_predictions %>% 
  inner_join(district_total_pop, by = "Dist_ID") 


# calculate pixel-level population estimates

population_predictions <- population_predictions %>% 
  mutate(predicted_pop = (predicted_exp/predicted_exp_sum)*district_pop)


# Predicted Population Diagnostics 

#Sum each pixel population totals to see if it matches district totals

test <- population_predictions %>% 
  group_by(Dist_ID) %>% 
  summarise(dist_pop = sum(predicted_pop)) %>% 
  ungroup() %>% 
  inner_join(GHA_df, by = "Dist_ID") %>% 
  select(dist_pop, Pop_2021)

# test if estimates match district population totals
all(test$Pop_2021 == round(test$dist_pop))

# Mapping population estimates --------------------------------------------

# Rasterize Predictions 

# Load grid_ids
Grid_ID <-rast(paste0(raster_path, "Grid_ID.tif"))

#plot(Grid_ID)

Grid_Pop <- terra::values(Grid_ID, dataframe = T) %>% 
  filter(!is.na(Grid_ID)) 

#check for NA values
any(is.na(Grid_Pop))

#check for Grid_ID duplicates
any(duplicated(Grid_Pop$Grid_ID))


#Join estimated population  to right Grid ID

Pixel_Estimates <- Grid_Pop %>% 
  left_join(population_predictions, by = "Grid_ID")

#Assign predictions to Grid Raster

Grid_ID[]<-Pixel_Estimates$predicted_pop

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "RF_pop_predict.tif"), overwrite=T,
            names = "RF_Pop_population")


################### End of Script ######################################################






