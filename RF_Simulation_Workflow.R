library(randomForest)
library(feather)
library(tidyverse)
library(tictoc)
library(terra)
library(VIM)
library(groupdata2)
library(lsa)


#Specify Drive Path
drive_path <- "//worldpop.files.soton.ac.uk/worldpop/Projects/WP517763_GRID3/"
input_path <-  paste0(drive_path, "Working/GHA/Ortis/Output/")
output_path1 <- paste0(drive_path, "Working/GHA/Ortis/Output/Simulated Posterior/")
output_path2 <- paste0(drive_path, "Working/GHA/Ortis/Output/Simulated Posterior/Predicted Population/")
raster_path <- paste0(drive_path, "Working/GHA/Ortis/Other_covariates/")

#Load Simulated Data

simu_data <- read_feather(paste0(input_path, "simu_district_pop.feather"))

head(simu_data[,1:5]) # only showing first five columns


#Define response variable as log of pop_density

 simu_data<- simu_data %>% 
  mutate(pop_density = log(pop/Area)) 

#Covariates selection
covs <- simu_data %>% 
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
simu_data2<-simu_data %>% 
  select(pop, pop_density, dist_id, Area) %>% 
  cbind(covs)


# Fit model to all the training data -------------------------------------------------

set.seed(4567)

# Search for best hyperparameters tunning

tuneRF(x = covs, y = simu_data2$pop_density, na.action = na.omit,
plot = T, trace = T, importance=TRUE,  sampsize=length(simu_data2), replace=TRUE) 


#Fit model
model1 <- randomForest(x = covs, y = simu_data2$pop_density, mtry = 16, na.action = na.omit, 
                      plot = T, trace = T, importance=TRUE, sampsize=length(simu_data), replace=TRUE) 


model1


#train data predictions
model1_predictions <- model1$predicted %>% as_tibble()

#cbind predicted data to original data
model1_predictions <- model1_predictions %>% 
  cbind(simu_data2$pop_density) %>% 
  mutate(observed = exp(simu_data2$pop_density), predicted = exp(value),
         residual = predicted - observed,
         residual1 = value - simu_data2$pop_density,
         model = "RandomForest")

write.csv(model1_predictions, paste0(output_path1, "Simu RF model results.csv"))

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
original_names <- read.csv(paste0(input_path, "var_names.csv"))

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


write.csv(var_importance_df, paste0(output_path1, "Simu_RF_var_importance.csv"))

#plot variable importance 
ggplot(var_importance_df, aes(x = reorder(Original.Name, inc_prop), y = inc_prop, fill = Original.Name)) +
  geom_bar(stat = "identity")+
  #geom_text(aes(label = round(inc_prop,3)), hjust=-0.2, size=5) + # add y values as labels to the bars
  geom_text(aes(label = round(inc_prop,2)), position = position_stack(vjust = 0.8), size=4, color = "white") + # add y values as labels to the bars inside the bars
  coord_flip() +
  theme_bw()+
  labs(x = "Variables", y = "Variable Importance(%)")+
  scale_fill_manual(values = c(rep("#8c2981", 27))) + # use only one color for the fill
  theme(legend.position="none",  axis.text.y = element_text(size = 14))


# Perform Out-of-Sample Cross Validation ----------------------------------

#Cross Validation using Training and Test data

#Create training dataset
train <- simu_data %>% 
  sample_frac(.70)

#Create test set
test <- anti_join(simu_data, train, by = "dist_id")

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


# Covariate Processing ----------------------------------------------------

settled_df <- read_feather(paste0(input_path, "simu_data.feather"))

#We have missing values in the covariates. Hence need to do imputation
#Testing three imputations methods and comparing metrics


# First method - mean imputation -----------------------------------------------------------
#Select covariates
covs1 <- settled_df %>% 
  select(starts_with("x"))

#apply mean imputation
covs1 <- covs1 %>%
  mutate_all(~as.numeric(.x)) %>% 
  mutate_all(~replace_na(.x, mean(.x, na.rm = TRUE)))

head(covs1)


# Second Imputation - Replacing NA values with 0 ---------------------------

#Select covariates
covs1 <- settled_df %>% 
  select(starts_with("x"))

#Replace NA values with 0
covs1 <- covs1 %>% 
  mutate_all(~replace_na(.x, 0))


head(covs1)

# Third method KNN Imputation (Mean) --------------------------------------


#Select covariates
#covs1 <- settled_df %>% 
 # select(starts_with("x"))


#Create a grouping variable for subsetting data in chunks
#covs1 <- covs1 %>% 
  #group(n = 100000, method = "greedy", col_name = "Group_ID") %>% 
  #ungroup() 


# split the data into chunks based on the Group_ID
#covs1 <- covs1 %>% 
  #group_split(Group_ID)


# Function to impute missing values using KNN

#tic()

#for(dd in covs1){
  # get the ID of the current chunk being processed
 # typro <- unique(dd$Group_ID)
  #print(typro)
  
  #Impute NA values
  #dd <- kNN(dd, numFun = mean, k = 4)
  
  #Write each group to file
  
 # write_feather(dd, paste0(output_path1, "KNN_", unique(dd$Group_ID), ".feather"))
  
#} 

#toc()

#load imputed data back to memory

#specify pattern for file names
#pattern = "KNN.*\\.feather$"

#list all files that match the pattern
#myfiles <-dir(output_path1,pattern=pattern)
#myfiles


#read files and rbind them

#tic()

#covs1 <- myfiles %>% 
  #map(function(x) read_feather(file.path(output_path1, x))) %>% 
  #reduce(rbind) %>% 
  #select(x1:x24)

#toc() 

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

# Weighting Layer Analysis (Predictions) ------------------------------------------------

predicted <- predict(model1, newdata = covs1)

#cbind predictions to settled_df

population_predictions <- settled_df %>% 
  select(grid_id, dist_id) %>% 
  cbind(predicted) %>% 
  mutate(predicted_exp = exp(predicted)) # back-transform predictions to natural scale

#Sum exponentiated predictions for all pixels in a given district
population_predictions <- population_predictions %>% 
  group_by(dist_id) %>% 
  mutate(predicted_exp_sum = sum(predicted_exp)) %>% 
  ungroup()

#Select Population totals from simu_data and merge with predicted population df
district_total_pop <- simu_data %>% 
  select(pop, dist_id)%>% 
  rename(district_pop = pop)

#Join to population
population_predictions <-population_predictions %>% 
  inner_join(district_total_pop, by = "dist_id") 


# calculate pixel-level population estimates

population_predictions <- population_predictions %>% 
  mutate(predicted_pop = (predicted_exp/predicted_exp_sum)*district_pop)


# Predicted Population Diagnostics 

#Sum each pixel population totals to see if it matches district totals

test <- population_predictions %>% 
  group_by(dist_id) %>% 
  summarise(dist_pop = sum(predicted_pop)) %>% 
  ungroup() %>% 
  inner_join(simu_data, by = "dist_id") %>% 
  select(dist_pop, pop)

# test if estimates match district population totals
all(test$pop == round(test$dist_pop))

#########################################################################
#########################################################################

# Check predicted population and observed

population_predictions1 <- population_predictions %>% 
  select(grid_id, predicted_pop) %>% 
  left_join(settled_df, by = "grid_id") %>% 
  select(predicted_pop, pop, grid_id) %>% 
  rename(predicted = predicted_pop, observed = pop)

write_feather(population_predictions1, paste0(output_path1, "Simu_RF_predictions.feather"))

#Calculate goodness of fit metrics
population_predictions1 %>% 
  mutate(residual = predicted - observed) %>% 
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            rmse = sqrt(mse),
            corr = cor(predicted, observed), 
            cosine_sim = cosine(predicted, observed))


#Plot Predicted vrs observed population

#ggplot(population_predictions1) +
  #geom_point(aes(x = observed, y = predicted), fill = 'grey50', color = 'grey70', shape = 21) +
  #geom_smooth(aes(x = observed, y = predicted), method = lm) + 
  #geom_abline(slope = 1, intercept = 0, color = 'orange', linewidth = 1) +
  #theme_minimal() +
  #labs(title = '', x = 'Observed population density', y = 'Predicted population density')


#Histogram Plot

population_predictions1 %>% 
  pivot_longer(cols = c(predicted, observed), names_to = "Population", 
               values_to = "predicted_population") %>%
  # filter(predicted_population <0.2) %>% 
  ggplot() +aes(predicted_population, fill = Population) +
  geom_histogram(alpha = 0.5, bins = 100)+
  labs(title = 'Histogram plot', x='Predicted Population', y='Frequency')+ 
  theme_bw()+
  scale_fill_discrete(name="Population",
                      breaks=c("predicted", "observed"),
                      labels=c("Predicted Population ", "Observed Population"))
#Density Plot
population_predictions1 %>% 
  pivot_longer(cols = c(predicted, observed), names_to = "Population", 
               values_to = "predicted_population") %>%
  filter(predicted_population <30) %>% 
  ggplot() +aes(predicted_population, fill = Population) +
  geom_density(alpha = 0.5)+
  labs(title = 'Histogram plot', x='Predicted Population', y='Frequency')+ 
  theme_bw()+
  scale_fill_discrete(name="Population",
                      breaks=c("predicted", "observed"),
                      labels=c("Predicted Population ", "Observed Population"))

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
  left_join(population_predictions, by = c("Grid_ID" = "grid_id"))

#Assign predictions to Grid Raster

Grid_ID[]<-Pixel_Estimates$predicted_pop

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "simu_random_pop_predict.tif"), overwrite=T)


############ End ##################################################################
