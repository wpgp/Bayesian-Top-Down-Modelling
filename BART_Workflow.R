
options(java.parameters = "-Xmx100g")

#library(caret)
library(tidyverse)
library(bartMachine)
library(tictoc)
library(feather)
library(groupdata2)
library(terra)

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

head(covs[,1:2]) # only showing first two columns



# Fit model to all the training data -------------------------------------------------

# Search for best hyperparameters tunning

 #bartMachineCV(X = covs, y = GHA_df$pop_density, seed = 1234, use_missing_data = T)

#Fit model
model1 <- bartMachine(X = covs, y = GHA_df$pop_density,
                      k = 5, nu = 3, q = 0.9, num_trees = 200, seed = 1234,
                      use_missing_data = T)

model1

#Check for model convergence
plot_convergence_diagnostics(model1)

#model upper and lower CI

model1_CI<- calc_credible_intervals(model1, new_data = covs)
model1_CI <- model1_CI %>% 
  as_tibble() %>% 
  mutate(ci_lower_bd = exp(ci_lower_bd), ci_upper_bd = exp(ci_upper_bd))


#predicted values
model1_predictions <- model1$y_hat_train %>% as_tibble()


#cbind predicted posteriors to original data
model1_predictions <- model1_predictions %>% 
  cbind(GHA_df$pop_density, model1_CI) %>% 
  mutate(observed = exp(GHA_df$pop_density), predicted = exp(value),
         residual = predicted - observed,
         model = "BART")
#write.csv(model1_predictions, paste0(output_path, "BART model results.csv"))

# compute goodness-of-fit metrics
model1_predictions %>%
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            #rmse = sqrt(mse),
            corr = cor(predicted, observed),
            In_IC = mean(observed<ci_upper_bd & observed>ci_lower_bd)*100)


# Model plots -------------------------------------------------------------
#Plot Population Density Estimate
ggplot(model1_predictions) +
  geom_pointrange(aes(x=observed, y=predicted, ymin=ci_lower_bd, ymax=ci_upper_bd
  ),
  fill='grey50', color='grey70', shape=21
  )+
  geom_abline(slope=1, intercept = 0, color='orange', linewidth=1)+
  theme_minimal()+
  labs(title = '', x='Observed Population Density', y='Predicted density')


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


# Variable importance plot 

#Read variable names
original_names <- read.csv(paste0(output_path, "var_names.csv"))

#Variable importance
var_importance <- investigate_var_importance(model1)

var_names <- names(var_importance$avg_var_props)[grep("^x", names(var_importance$avg_var_props))]
var_importance_df <- data.frame(variable = var_names, inc_prop = var_importance$avg_var_props[var_names])

#Join var_importance-df to var_names
var_importance_df <- var_importance_df %>% 
  inner_join(original_names, by = c("variable" = "var_names2")) %>% 
  mutate(Model = "BART")

#write.csv(var_importance_df, paste0(output_path, "BART_var_importance.csv"))

#plot variable importance 
ggplot(var_importance_df, aes(x = reorder(Original.Name, inc_prop), y = inc_prop, fill = Original.Name)) +
  geom_bar(stat = "identity")+
  geom_text(aes(label = round(inc_prop,3)), hjust=-0.2, size=3) + # add y values as labels to the bars
  coord_flip() +
  theme_bw()+
  labs(x = "Variables", y = "Inclusion Proportion(%)")+
  scale_fill_manual(values = c(rep("#8c2981", 27))) + # use only one color for the fill
  theme(legend.position="none")


# Perform Out-of-Sample Cross Validation ----------------------------------

#Cross Validation using Training and Test data

set.seed(4567)

#Create training dataset
train <- GHA_df %>% 
  sample_frac(.70)

#Create test set
test <- anti_join(GHA_df, train, by = "Dist_ID")

#train covariates
covs <- train %>% 
  select(starts_with("x"))

#fit model to train dataset

model2 <- bartMachine(X = covs, y = train$pop_density,
                      k = 5, nu = 10, q = 0.75, num_trees = 200, seed = 1234)


model2

#Make predictions on the test data
covs <- test %>% 
  select(starts_with("x"))

predicted <- predict(model2, new_data = covs)

#cbind to test data
test <- test %>% 
  select(pop_density) %>% 
  cbind(predicted)%>% 
  mutate(observed = exp(pop_density), predicted = exp(predicted),
         residual = predicted - observed)


# compute goodness-of-fit metrics
test %>% 
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            rmse = sqrt(mse),
            corr = cor(predicted, observed))


rm(covs, train, test)

# Weighting Layer Analysis ------------------------------------------------

settled_df <- read_feather(paste0(input_path, "settled_df.feather"))

# Create a grouping variable for subsetting data
settled_df <- settled_df %>% 
  group(n = 100000, method = "greedy", col_name = "Group_ID") %>%
  ungroup()
        
# split the prediction data by Group_ID
covs_test_list <- settled_df %>% 
  group_split(Group_ID)

# create a function to make predictions 
make_predictions <- function(df) {
  # get the ID of the current region being processed
  typro <- unique(df$Group_ID)
  print(typro)
  
  # extract id of current area
  covs_id <- df %>% 
    select(Grid_ID, Dist_ID, Group_ID)
  
  covs <- df %>% 
    select(-c(Grid_ID, Dist_ID, Group_ID))
  
  # make predictions on the settled_df split
  covs$predicted <- predict(model1, new_data = covs)
  
  # bind predictions
  covs <- covs %>% 
    select(predicted) %>% 
    cbind(covs_id)
  
  write_feather(covs, paste0(output_path1, "Group_", unique(df$Group_ID), ".feather"))
  #return(covs)
  
}

tic()

# apply the function to the list of splitted dataframes
covs_results <- map(covs_test_list, make_predictions)

toc()

#Read files back into memory
tic()

myfiles <-dir(output_path1,pattern="*.feather")

covs_predictions <- myfiles %>% 
  map(function(x) read_feather(file.path(output_path1, x))) %>% 
  reduce(rbind) 

toc()


#Join predictions to Original test data
 covs_test_predictions <-settled_df %>% 
  select(Grid_ID) %>% 
  inner_join(covs_predictions, by = c("Grid_ID" = "Grid_ID"))%>% 
  mutate(predicted_exp = exp(predicted)) # back-transform predictions to natural scale


#Sum exponentiated predictions for all pixels in a given district
predicted_muni_totals <- covs_test_predictions %>% 
  group_by(Dist_ID) %>% 
  summarise(predicted_exp_sum = sum(predicted_exp)) %>% 
  ungroup()

#join total to covs_test_predictions
covs_test_predictions<- covs_test_predictions %>% 
  inner_join(predicted_muni_totals, by = c("Dist_ID" = "Dist_ID"))

#Select Population totals from GHA_df and merge with covs_test_predictions
municipal_total_pop <- GHA_df %>% 
  select(Pop_2021, Dist_ID)%>% 
  rename(pop_municipality = Pop_2021)

covs_test_predictions <-covs_test_predictions %>% 
  inner_join(municipal_total_pop, by = "Dist_ID") 


# calculate pixel-level population estimates

covs_test_predictions <-covs_test_predictions %>% 
  mutate(predicted_pop = (predicted_exp/predicted_exp_sum)*pop_municipality)


# Predicted Population Diagnostics 

#Sum each EA population totals to see if it matches municipal totals

test <- covs_test_predictions %>% 
  group_by(Dist_ID) %>% 
  summarise(muni_pop = sum(predicted_pop)) %>% 
  ungroup() %>% 
  inner_join(GHA_df, by = "Dist_ID") %>% 
  select(muni_pop, Pop_2021)

# test if estimates match municipality population totals
all(test$Pop_2021 == round(test$muni_pop))


# Method to get posteriors------------------

#Select test data covariates for making posterior predictions

# extract the columns needed for prediction
covs_test <- settled_df %>% 
  select(Grid_ID, Dist_ID, starts_with("x")) %>% 
  inner_join(municipal_total_pop, by = "Dist_ID") %>% 
  group(n = 100000, method = "greedy", col_name = "Group_ID") %>% 
  ungroup()

# split the prediction data by Group_ID
covs_test_list <- covs_test %>% 
  group_split(Group_ID)


# create a function to calculate the posterior predictions and bind the results with covs_id

get_posteriors <- function(df) {
  # get the ID of the current region being processed
  typro <- unique(df$Group_ID)
  print(typro)
  
  # extract id of current area
  covs_id <- df %>% 
    select(Grid_ID, Dist_ID, Group_ID, pop_municipality)
  
  covs <- df %>% 
    select(-c(Grid_ID, Dist_ID, Group_ID, pop_municipality))
  
  # make predictions on the test data
  get_posteriors <- bart_machine_get_posterior(model1, covs)
  
  #Back transform predicted posteriors 
  get_posteriors <- exp(get_posteriors$y_hat_posterior_samples) %>% 
    as_tibble()
  
  # bind covs_id to posteriors 
  get_posteriors <- get_posteriors %>% 
   cbind(covs_id)
  
  write_feather(get_posteriors, paste0(output_path2, "Group_", unique(df$Group_ID), ".feather"))
 # return(get_posteriors)
}

tic()

# apply the function to the list of splitted dataframes
posteriors_results <- map(covs_test_list, get_posteriors)

toc()

#Read files back into memory
tic()

myfiles <-dir(output_path2,pattern="*.feather")

posterior_predictions <- myfiles %>% 
  map(function(x) read_feather(file.path(output_path2, x))) %>% 
  reduce(rbind) 

toc()

## Combine the posteriors into a single dataframe using cbind
#posterior_predictions <- do.call(rbind, posteriors_results) 

# Function to calculate posteriors population estimates for each pixel ---------------------------------------------

# create an empty list to store the predictions
results <- list()

#get id and mean predicted population
pred_id <- covs_test_predictions %>% 
  select(Grid_ID, Dist_ID, predicted_pop)

tic()

# Loop through each posterior prediction variable (V)
for (i in 1:1000) {
  
 data <- posterior_predictions %>% 
    select(Grid_ID, Dist_ID, pop_municipality, !!paste0("V", i))   # Select the geocode and V variable for the current iteration
  
  
  # Group by geo_code and calculate the total for the current V variable
  total <- data %>% 
    group_by(Dist_ID) %>% 
    summarize(total = as.numeric(sum(!!sym(paste0("V", i)))))
  
  #Join total to data
  data <- data %>% 
    inner_join(total, by = "Dist_ID")
  
  # Calculate the prediction for the current V variable
  data <- data %>% 
    mutate(!!paste0("prediction_", i) := (!!sym(paste0("V", i))) / total * pop_municipality) %>%
    select(!!sym(paste0("prediction_", i))) 
  # Add the prediction results to the list
  results[[i]] <- data
  
}

# Combine the predictions into a single dataframe using cbind
results_df <- do.call(cbind, results)

toc()

tic() 

#mean_population <- rowMeans(results_df, na.rm = T)
#median_population <- apply(results_df, 1, FUN = function(x) quantile(x, probs = 0.5, na.rm = T))
std_population <- apply(results_df, 1, sd)
lower_quantile <- apply(results_df, 1, FUN = function(x) quantile(x, probs = 0.025, na.rm=T))
upper_quantile <- apply(results_df, 1, FUN = function(x) quantile(x, probs = 0.975, na.rm = T))

toc()

#cbind results to ids
predicted_population <- cbind(pred_id, lower_quantile, upper_quantile, std_population) %>% 
                        select(Grid_ID, Dist_ID, lower_quantile, predicted_pop, upper_quantile, everything()) %>% 
  mutate(uncertainty =(upper_quantile - lower_quantile)/predicted_pop, 
         coe_var = std_population/predicted_pop) 

#predicted_population <- format(predicted_population, scientific = F)

#write.csv(predicted_population, paste0(output_path2, "predicted_population.csv"))
predicted_population<- read.csv(paste0(output_path2, "predicted_population.csv"))

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
  left_join(predicted_population, by = "Grid_ID")

#Assign predictions to Grid Raster

Grid_ID[]<-Pixel_Estimates$predicted_pop

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "pop_predict.tif"), overwrite=T)

# Upper bound ------------------------------------------------
Upper_Bound <- Grid_Pop %>% 
  left_join(predicted_population, by = "Grid_ID")


Grid_ID[]<- Upper_Bound$upper_quantile

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "upper_pop_predict.tif"), overwrite=T)


# Lower Bound ------------------------------------------------------

Lower_Bound <- Grid_Pop %>% 
  left_join(predicted_population, by = "Grid_ID")

#Assign predictions to Grid Raster

Grid_ID[]<- Lower_Bound$lower_quantile

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "lower_pop_predict.tif"), overwrite=T)


# Coefficient of Variation ------------------------------------------------

coe_var <- Grid_Pop %>% 
  left_join(predicted_population, by = "Grid_ID")

#Assign predictions to Grid Raster

Grid_ID[]<- coe_var$coe_var

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "coe_pop_predict.tif"), overwrite=T)

# Uncertainty------------------------------------------------

uncertainty <- Grid_Pop %>% 
  left_join(predicted_population, by = "Grid_ID")

#Assign predictions to Grid Raster

Grid_ID[]<- uncertainty$uncertainty

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "uncertainty_pop_predict.tif"), overwrite=T)

#We can plot the other variables in same way
#####################End of Predictions########################################





