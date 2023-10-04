
options(java.parameters = "-Xmx100g")


library(tidyverse)
library(bartMachine)
library(tictoc)
library(feather)
library(groupdata2)
library(terra)
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
covs <- apply(covs, 2, stdize) %>%    #z-score
  as_tibble()

head(covs[,1:2]) # only showing first two columns

#Select pop_density and cbind covs
simu_data2<-simu_data %>% 
  select(pop, pop_density, dist_id, Area) %>% 
  cbind(covs) 

# Fit model to all the training data -------------------------------------------------
set.seed(4567)

# Search for best hyperparameters tunning

#bartMachineCV(X = covs, y = simu_data2$pop_density, use_missing_data = T)

#Fit model
model1 <- bartMachine(X = covs, y = simu_data2$pop_density,
                      k = 5, nu = 10, q = 0.75, num_trees = 200,
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
  cbind(simu_data$pop_density, model1_CI) %>% 
  mutate(observed = exp(simu_data$pop_density), predicted = exp(value),
         residual = predicted - observed,
         model = "BART")
write.csv(model1_predictions, paste0(output_path1, "Simu BART model results.csv"))

##Check if predicted value falls between lower and upper intervals
with(model1_predictions, all(predicted >= ci_lower_bd & predicted <= ci_upper_bd))

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
  filter(predicted_density <0.25) %>% 
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
  filter(predicted_density <0.25) %>% 
  ggplot() +aes(predicted_density, fill = Density) +geom_density(alpha = 0.5)+
  labs(title = 'Density plot', x='Predicted Population Density', y='Frequency')+ 
  theme_bw()+
  scale_fill_discrete(name="Density",
                      breaks=c("predicted", "observed"),
                      labels=c("Predicted Population Density", "Observed Population Density"))


# Variable importance plot 
original_names <- read.csv(paste0(input_path, "var_names.csv"))

var_importance <- investigate_var_importance(model1)

var_names <- names(var_importance$avg_var_props)[grep("^x", names(var_importance$avg_var_props))]
var_importance_df <- data.frame(variable = var_names, inc_prop = var_importance$avg_var_props[var_names])

#Join var_importance-df to var_names
var_importance_df <- var_importance_df %>% 
  inner_join(original_names, by = c("variable" = "var_names2")) %>% 
  mutate(Model = "BART", inc_prop = 100*inc_prop)

write.csv(var_importance_df, paste0(output_path1, "Simu_BART_var_importance.csv"))

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

rm(covs, model1_CI, model1_predictions, original_names, var_importance_df, var_importance); gc()


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

#Scale covariates with train data mean and sd
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

settled_df <- read_feather(paste0(input_path, "simu_data.feather"))


#Select covariates
covs1 <- settled_df %>% 
  select(starts_with("x"))

# Scale covariates using means and standard deviations from covs_stat

#Scale covariates
for (var in names(covs1)) {
  var_mean <- cov_stats$Mean[cov_stats$Covariate == var]
  var_sd <- cov_stats$Std_Dev[cov_stats$Covariate == var]
  covs1[[var]] <- (covs1[[var]] - var_mean) / var_sd
}

# Viewing the scaled dataframe
head(covs1)

#cbind scaled covariates to data & Create a grouping variable for subsetting data
settled_df1 <- settled_df %>% 
  select(dist_id, grid_id) %>% 
  cbind(covs1) %>% 
  group(n = 100000, method = "greedy", col_name = "Group_ID") %>%  
  ungroup()

# split the prediction data by Group_ID
covs_test_list <- settled_df1 %>% 
  group_split(Group_ID)

# create a function to make predictions 
make_predictions <- function(df) {
  # get the ID of the current region being processed
  typro <- unique(df$Group_ID)
  print(typro)
  
  # extract id of current area
  covs_id <- df %>% 
    select(grid_id, dist_id, Group_ID)
  
  covs <- df %>% 
    select(-c(grid_id, dist_id, Group_ID))
  
  # make predictions on the settled_df split
  covs$predicted <- predict(model1, new_data = covs)
  
  # bind predictions
  covs <- covs %>% 
    select(predicted) %>% 
    cbind(covs_id)
  
  write_feather(covs, paste0(output_path1, "Group_", unique(df$Group_ID), ".feather"))
  #return(covs)
  
}

rm(covs_test_list); gc();

tic()

# apply the function to the list of splitted dataframes
covs_results <- map(covs_test_list, make_predictions)

toc()

#Read files back into memory

#specify pattern for file names
pattern = "Group_.*\\.feather$"

tic()

myfiles <-dir(output_path1,pattern= pattern)

covs_predictions <- myfiles %>% 
  map(function(x) read_feather(file.path(output_path1, x))) %>% 
  reduce(rbind) 

toc()


#Join predictions to Original test data
covs_test_predictions <-settled_df1 %>% 
  select(grid_id) %>% 
  inner_join(covs_predictions, by = c("grid_id" = "grid_id"))%>% 
  mutate(predicted_exp = exp(predicted)) # back-transform predictions to natural scale


#Sum exponentiated predictions for all pixels in a given district
predicted_muni_totals <- covs_test_predictions %>% 
  group_by(dist_id) %>% 
  summarise(predicted_exp_sum = sum(predicted_exp, na.rm = T)) %>% 
  ungroup()

#join total to covs_test_predictions
covs_test_predictions<- covs_test_predictions %>% 
  inner_join(predicted_muni_totals, by = c("dist_id" = "dist_id"))

#Select Population totals from simu_data and merge with covs_test_predictions
municipal_total_pop <- simu_data %>% 
  select(pop, dist_id)%>% 
  rename(pop_municipality = pop)

covs_test_predictions <-covs_test_predictions %>% 
  inner_join(municipal_total_pop, by = "dist_id") 


# calculate pixel-level population estimates

covs_test_predictions <-covs_test_predictions %>% 
  mutate(predicted_pop = (predicted_exp/predicted_exp_sum)*pop_municipality) %>% 
  select(-predicted)


# Predicted Population Diagnostics 

#Sum each district population totals to see if it matches municipal totals

test <- covs_test_predictions %>% 
  group_by(dist_id) %>% 
  summarise(muni_pop = sum(predicted_pop)) %>% 
  ungroup() %>% 
  inner_join(simu_data, by = "dist_id") %>% 
  select(muni_pop, pop)

# test if estimates match municipality population totals
all(test$pop == round(test$muni_pop))

rm(cov_stats, covs_predictions, predicted_muni_totals); gc()


# compare predicted vrs observed predictions 
#get simulated pop
pred_id <- settled_df %>% 
  select(pop)

#cbind simulated Pop to predictions
covs_test_predictions <- covs_test_predictions %>% 
  cbind(pred_id)

predicted_population1 <- covs_test_predictions %>% 
  rename(predicted = predicted_pop, observed = pop) %>% 
  mutate(residual = predicted - observed)

write_feather(predicted_population1, paste0(output_path1, "Simu_BART_predictions.feather"))

#Calculate goodness of fit metrics
predicted_population1 %>% 
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

predicted_population1 %>% 
  pivot_longer(cols = c(predicted, observed), names_to = "Population", 
               values_to = "predicted_population") %>%
  filter(predicted_population <250) %>% 
  ggplot() +aes(predicted_population, fill = Population) +
  geom_histogram(alpha = 0.5, bins = 100)+
  labs(title = 'Histogram plot', x='Predicted Population', y='Frequency')+ 
  theme_bw()+
  scale_fill_discrete(name="Population",
                      breaks=c("predicted", "observed"),
                      labels=c("Predicted Population ", "Observed Population"))
#Density Plot
predicted_population1 %>% 
  pivot_longer(cols = c(predicted, observed), names_to = "Population", 
               values_to = "predicted_population") %>%
  filter(predicted_population <250) %>% 
  ggplot() +aes(predicted_population, fill = Population) +
  geom_density(alpha = 0.5)+
  labs(title = 'Histogram plot', x='Predicted Population', y='Frequency')+ 
  theme_bw()+
  scale_fill_discrete(name="Population",
                      breaks=c("predicted", "observed"),
                      labels=c("Predicted Population ", "Observed Population"))

################ End ##################################################################

