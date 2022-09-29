suppressWarnings(if(!require("pacman")) install.packages("pacman"))
pacman::p_load(tidyverse, tidymodels)

library(tidyverse)

library(tidymodels)

# Import the data set
diabetes <- read_table2(file = "https://www4.stat.ncsu.edu/~boos/var.select/diabetes.rwrite1.txt")


# Get a glimpse and dimensions of the data
glimpse(diabetes)


# Select the first 5 rows of the data
diabetes %>% 
  slice(1:5)

# Select predictor feature `bmi` and outcome `y`
diabetes_select <- diabetes %>% 
  select(c(bmi, y))

# Print the first 5 rows
diabetes_select %>% 
  slice(1:10)


set.seed(2056)
# Split 67% of the data for training and the rest for tesing
diabetes_split <- diabetes_select %>% 
  initial_split(prop = 0.67)

# Extract the resulting train and test sets
diabetes_train <- training(diabetes_split)
diabetes_test <- testing(diabetes_split)

# Print the first 3 rows of the training set
diabetes_train %>% 
  slice(1:10)


# Build a linear model specification
lm_spec <- 
  # Type
  linear_reg() %>% 
  # Engine
  set_engine("lm") %>% 
  # Mode
  set_mode("regression")


# Print the model specification
lm_spec

# Build a linear model specification
lm_spec <- linear_reg() %>% 
  set_engine("lm") %>%
  set_mode("regression")


# Train a linear regression model
lm_mod <- lm_spec %>% 
  fit(y ~ ., data = diabetes_train)

# Print the model
lm_mod

# Make predictions for the test set
predictions <- lm_mod %>% 
  predict(new_data = diabetes_test)

# Print out some of the predictions
predictions %>% 
  slice(1:5)

# Combine the predictions and the original test set
results <- diabetes_test %>% 
  bind_cols(predictions)


results %>% 
  slice(1:5)

library(ggplot2)
# Set a theme for the plot
theme_set(theme_light())
# Create a scatter plot
results %>% 
  ggplot(aes(x = bmi)) +
  # Add a scatter plot
  geom_point(aes(y = y), size = 1.6) +
  # Add a line plot
  geom_line(aes(y = .pred), color = "blue", size = 1.5)

