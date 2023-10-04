
options(java.parameters = "-Xmx100g")

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

# Calculate mean and standard deviation of covariates
cov_stats <- data.frame(Covariate = colnames(covs),
                        Mean = apply(covs, 2, mean, na.rm = TRUE),
                        Std_Dev = apply(covs, 2, sd, na.rm = TRUE))

#Scaling function to scale covariates
stdize <- function(x)
{ stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
return(stdz) }

#apply scaling function
covs <- apply(covs, 2, stdize) %>%    #z-score
  as_tibble()

head(covs[,1:2]) # only showing first two columns

#Select pop_density and cbind covs
GHA_df2 <-GHA_df %>% 
  select(Pop_2021, pop_density, Dist_ID, Area) %>% 
  cbind(covs) 

# Fit model to all the training data -------------------------------------------------

set.seed(4567)

# Search for best hyperparameters tunning

#bartMachineCV(X = covs, y = GHA_df2$pop_density, use_missing_data = T)

#Fit model
model1 <- bartMachine(X = covs, y = GHA_df2$pop_density,
                      k = 5, nu = 10, q = 0.75, num_trees = 200, use_missing_data = T)

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
  cbind(GHA_df2$pop_density, model1_CI) %>% 
  mutate(observed = exp(GHA_df2$pop_density), predicted = exp(value),
         residual = predicted - observed,
         model = "BART")
write.csv(model1_predictions, paste0(output_path, "BART model results.csv"))

# compute goodness-of-fit metrics
model1_predictions %>%
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            rmse = sqrt(mse),
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
  mutate(Model = "BART", inc_prop = 100*inc_prop)

write.csv(var_importance_df, paste0(output_path, "BART_var_importance.csv"))

#plot variable importance 
ggplot(var_importance_df, aes(x = reorder(Original.Name, inc_prop), y = inc_prop, fill = Original.Name)) +
  geom_bar(stat = "identity")+
  geom_text(aes(label = round(inc_prop,3)), hjust=-0.2, size=5) + # add y values as labels to the bars
  coord_flip() +
  theme_bw()+
  labs(x = "Variables", y = "Variable Importance(%)")+
  scale_fill_manual(values = c(rep("#8c2981", 27))) + # use only one color for the fill
  theme(legend.position="none", axis.text.y = element_text(size = 14))


#plot variable importance 
ggplot(var_importance_df, aes(x = reorder(Original.Name, inc_prop), y = inc_prop, fill = Original.Name)) +
  geom_bar(stat = "identity")+
  #geom_text(aes(label = round(inc_prop,3)), hjust=-0.2, size=5) + # add y values as labels to the bars
  geom_text(aes(label = round(inc_prop,2)), position = position_stack(vjust = 0.5), size=4, color = "white") + # add y values as labels to the bars inside the bars
  coord_flip() +
  theme_bw()+
  labs(x = "Variables", y = "Variable Importance(%)")+
  scale_fill_manual(values = c(rep("#8c2981", 27))) + # use only one color for the fill
  theme(legend.position="none", axis.text.y = element_text(size = 14))

#remove variables
rm(covs, model1_CI, model1_predictions, original_names, var_importance_df, var_importance); gc()


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

model2 <- bartMachine(X = covs_train, y = train$pop_density,
                      k = 5, nu = 10, q = 0.75, num_trees = 200, use_missing_data = T)


model2

#model upper and lower CI
model2_CI<- calc_credible_intervals(model2, new_data = covs_train)
model2_CI <- model2_CI %>% 
  as_tibble() %>% 
  mutate(ci_lower_bd = exp(ci_lower_bd), ci_upper_bd = exp(ci_upper_bd))


#predicted values
model2_predictions <- model2$y_hat_train %>% as_tibble()


#cbind predicted posteriors to original data
model2_predictions <- model2_predictions %>% 
  cbind(train$pop_density, model2_CI) %>% 
  mutate(observed = exp(train$pop_density), predicted = exp(value),
         residual = predicted - observed)

# compute goodness-of-fit metrics on In-sample
model2_predictions %>%
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            rmse = sqrt(mse),
            corr = cor(predicted, observed),
            In_IC = mean(observed<ci_upper_bd & observed>ci_lower_bd)*100)

#Make predictions on the test data
covs_test <- test %>% 
  select(starts_with("x"))

#Scale covariates
for (var in names(covs_test)) {
  var_mean <- covs_train_stats$Mean[covs_train_stats$Covariate == var]
  var_sd <- covs_train_stats$Std_Dev[covs_train_stats$Covariate == var]
  covs_test[[var]] <- (covs_test[[var]] - var_mean) / var_sd
}


head(covs_test[,1:2]) 


predicted <- predict(model2, new_data = covs_test)

test_CI<- calc_credible_intervals(model2, new_data = covs_test)
test_CI <- test_CI %>% 
  as_tibble() %>% 
  mutate(ci_lower_bd = exp(ci_lower_bd), ci_upper_bd = exp(ci_upper_bd))

#cbind to test data
test <- test %>% 
  select(pop_density) %>% 
  cbind(predicted, test_CI)%>% 
  mutate(observed = exp(pop_density), predicted = exp(predicted),
         residual = predicted - observed)


# compute goodness-of-fit metrics

test %>% 
  summarise(Bias= mean(residual),
            Imprecision = sd(residual),
            Inaccuracy = mean(abs(residual)),
            mse = mean((residual)^2),
            rmse = sqrt(mse),
            corr = cor(predicted, observed),
            In_IC = mean(observed<ci_upper_bd & observed>ci_lower_bd)*100)

rm(covs_test, covs_train, train, test, model2, model2_CI, model2_predictions, test_CI, covs_train_stats);

# Weighting Layer Analysis ------------------------------------------------

settled_df <- read_feather(paste0(input_path, "settled_df.feather"))


# Scale covariates using means and standard deviations from covs_stat

covs1 <- settled_df %>% 
  select(starts_with("x"))

#Scale covariates
for (var in names(covs1)) {
  var_mean <- cov_stats$Mean[cov_stats$Covariate == var]
  var_sd <- cov_stats$Std_Dev[cov_stats$Covariate == var]
  covs1[[var]] <- (covs1[[var]] - var_mean) / var_sd
}

# Viewing the scaled dataframe
head(covs1)

#check for NA values
any(is.na(covs1))


#Add a grouping variable to covs1 and split data
covs1 <- covs1 %>% 
  group(n = 100000, method = "greedy", col_name = "Group_ID") %>% 
  ungroup() %>% 
  group_split(Group_ID)


# Method to Get Posteriors ------------------------------------------------

#Function to calculate the posteriors

get_posteriors <- function(df) {
  # get the ID of the current region being processed
  typro <- unique(df$Group_ID)
  print(typro)
  
  df <- df %>% 
   select(-Group_ID)
  
  # make predictions on the test data
  get_posteriors <- bart_machine_get_posterior(model1, df)
  
  #Back transform predicted posteriors 
  get_posteriors <- exp(get_posteriors$y_hat_posterior_samples) %>% 
    as_tibble()
  
  #Write predictions to file
    write_feather(get_posteriors, paste0(output_path2, "Group", typro, ".feather"))
    }

tic()

# apply the function to the list of splitted dataframes
posteriors_results <- map(covs1, get_posteriors)

toc()

#Read files back into memory

#specify pattern for file names
pattern = "Group.*\\.feather$"

tic()

myfiles <-dir(output_path2,pattern= pattern)

posterior_predictions <- myfiles %>% 
  map(function(x) read_feather(file.path(output_path2, x))) %>% 
  reduce(rbind) 

toc() #14 min

#Select Pop from Ghana_Data
Dist_Pop <- GHA_df %>% 
  select(Dist_ID, Pop_2021)
  

#get grid ids and join Dist_Pop data
id <- settled_df %>% 
  select(Dist_ID, Grid_ID) %>% 
  full_join(Dist_Pop, by = "Dist_ID")

         
#Cbind ids to the posteriors
posterior_predictions <- posterior_predictions %>% 
  cbind(id)


# calculate posteriors population estimates for each pixel 

# create an empty list to store the predictions
results <- list()

tic()

# Loop through each posterior prediction variable (V)
for (i in 1:1000) {
  
  data <- posterior_predictions %>% 
    select(Grid_ID, Dist_ID, Pop_2021, !!paste0("V", i))   # Select  current iteration (V)
  
  
  # Group by Dist_ID and calculate the total for the current V variable
  total <- data %>% 
    group_by(Dist_ID) %>% 
    summarize(total = as.numeric(sum(!!sym(paste0("V", i)))))
  
  #Join total to data
  data <- data %>% 
    inner_join(total, by = "Dist_ID")
  
  # Calculate the prediction for the current V variable
  data <- data %>% 
    mutate(!!paste0("prediction_", i) := (!!sym(paste0("V", i))) / total * Pop_2021) %>%
    select(!!sym(paste0("prediction_", i))) 
  # Add the prediction results to the list
  results[[i]] <- data
  
}

# Combine the predictions into a single dataframe using cbind
results_df <- do.call(cbind, results)

toc() #82.98 sec elapsed


#Summarize posteriors

tic() 

mean_population <- rowMeans(results_df, na.rm = T)
median_population <- apply(results_df, 1, FUN = function(x) quantile(x, probs = 0.5, na.rm = T))
std_population <- apply(results_df, 1, sd)
lower_quantile <- apply(results_df, 1, FUN = function(x) quantile(x, probs = 0.025, na.rm=T))
upper_quantile <- apply(results_df, 1, FUN = function(x) quantile(x, probs = 0.975, na.rm = T))

toc()  # 15 min

#cbind results to summarized posteriors
predicted_population <- cbind(id, lower_quantile, mean_population, median_population,
                              upper_quantile, std_population) %>% 
  select(Grid_ID, Dist_ID,Pop_2021, lower_quantile, mean_population, upper_quantile, everything()) %>% 
  mutate(uncertainty =(upper_quantile - lower_quantile)/mean_population, 
         coe_var = std_population/mean_population) 

#Check if predicted population falls between lower and upper intervals
with(predicted_population, all(mean_population >= lower_quantile & mean_population <= upper_quantile))

write_feather(predicted_population, paste0(output_path2, "predicted_population.feather"))
#predicted_population<- read_feather(paste0(output_path2, "predicted_population.feather"))


# Validate Predicted District Totals 

#Sum each district population totals to see if it matches predicted totals

test <- predicted_population %>% 
  group_by(Dist_ID, Pop_2021) %>% 
  summarise(dist_total_pop = sum(mean_population)) %>% 
  ungroup() 


# test if Pop_2021 match predicted population totals
all(test$Pop_2021 == round(test$dist_total_pop))

rm(cov_stats, covs1, data, Dist_Pop, GHA_df, GHA_df2, covs_predictions, predicted_muni_totals); gc()


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

Grid_ID[]<-Pixel_Estimates$mean_population

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "BART_predict.tif"), 
            overwrite=T, names = "BART_mean_population")

# Upper bound ------------------------------------------------
Upper_Bound <- Grid_Pop %>% 
  left_join(predicted_population, by = "Grid_ID")


Grid_ID[]<- Upper_Bound$upper_quantile

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "upper_BART_predict.tif"), 
            overwrite=T, names = "BART_upper_population")


# Lower Bound ------------------------------------------------------

Lower_Bound <- Grid_Pop %>% 
  left_join(predicted_population, by = "Grid_ID")

#Assign predictions to Grid Raster

Grid_ID[]<- Lower_Bound$lower_quantile

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "lower_BART_predict.tif"), 
            overwrite=T, names = "BART_lower_population")


# Coefficient of Variation ------------------------------------------------

coe_var <- Grid_Pop %>% 
  left_join(predicted_population, by = "Grid_ID")

#Assign predictions to Grid Raster

Grid_ID[]<- coe_var$coe_var

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "coe_BART_predict.tif"),
            overwrite=T, names = "BART_COE_population")

# Uncertainty------------------------------------------------

uncertainty <- Grid_Pop %>% 
  left_join(predicted_population, by = "Grid_ID")

#Assign predictions to Grid Raster

Grid_ID[]<- uncertainty$uncertainty

plot(Grid_ID)


writeRaster(Grid_ID, paste0(output_path2, "uncertainty_BART_predict.tif"), 
            overwrite=T, names = "BART_uncertainty_population")

#We can plot the other variables in same way
#####################End of Predictions########################################

# Visualizations ----------------------------------------------------------

#geom_bar populations
 predicted_population %>% 
    filter(Pop_2021 == 443981) %>% 
    sample_n(10) %>%
  arrange() %>% 
  mutate(Predictions = paste0("prediction ", 1:n())) %>% 
    
  ggplot()+
    aes(mean_population, reorder(Predictions, mean_population)) + 
    
    geom_errorbarh(aes(xmin = lower_quantile, xmax = upper_quantile, color = Predictions),
                   height = 0, size = 3) +
    geom_errorbarh(aes(xmin = lower_quantile, xmax = upper_quantile, color = Predictions),
                   height = 0.5) + 
    geom_point(size = 3, color = "#D55E00")+ 
  labs(x = "Mean Population", y = " ") +
  theme_minimal()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
  
    

#Geom_ribbon

predicted_population %>% 
  filter(Pop_2021 == 443981) %>% 
  
  ggplot()+
  aes(Grid_ID, mean_population)+
  geom_ribbon(aes(ymin = lower_quantile, ymax = upper_quantile), fill = "lightblue", alpha = 0.5) + # Plot the credible band as a ribbon
  geom_line(color = "blue") + # Plot the posterior estimate as a line
  labs(x = "X variable", y = "Posterior estimate", title = "Plot with geom_ribbon and geom_line")






################### End #############################################################################








