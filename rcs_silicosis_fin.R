library(ggplot2)
library(ggrepel)
library(tidyverse)
library(dosresmeta)
library(rms)
library(readxl)
library(survival)
library(flexsurv)
library(gridExtra)
library(ggpubr)
library(writexl)

setwd("/Users/phowlett/Library/CloudStorage/Dropbox/PhD/BSc project")

#############
# Load data #
#############

# For the cumulative risk analysis
risk <- read_excel("risk_complete.xlsx")

# For the drma
risk_drma <- read_excel("risk_drma_2.xlsx")

##########################################
# Comparison of life table cumulative risks#
##########################################

# Cum risk may be given by: 1-exp[-sum of (hazards*interval width)]

# hazards is given by: no. cases/(width*(no. entering category-(0.5*no. cases)-(0.5*no.withdrawals))
# Given that we do no know the number of withdrawals, this is assumed to be 0 
# We have already calculated the hazards (rate = cases / (width * (cum_cases - (cases * 0.5))), # rate within category i.e. cases/exposed person years
# Thus we can now calculate the cumulative risk for each study

# Split data into cross sectional analysis and longitudinal analysis
risk_cat <- risk %>% subset(., type == "cc")
table(risk_cat$author)

risk_py <- risk %>% subset(., type == "ir")
table(risk_py$author)

# For cumulative incidence
risk_cat <- risk_cat %>% 
  group_split(id) %>% 
  map(~ mutate(.x, 
               width = d_hi - d_low, # width of the category
               rate = cases / (width * (cum_cases-(n*0.5))), # rate within category i.e. cases/exposed person years, where exposed years are CUMULATIVE
               w_risk = cumsum(rate * width), # Calculation of cumulative hazard rate
               cum_risk = 1-exp(-w_risk),
               na = cases / (cum_cases-(n*0.5)),
               na_haz = cumsum(na), 
               prop = cases / n)) %>% # Calculation of cumulative risk up to dose-time
  bind_rows()

# For person years data
risk_py <- risk_py %>% 
  group_split(id) %>% 
  map(~ mutate(.x, 
               width = d_hi - d_low, 
               rate = cases / n, # incidence rate within a band
               rate2 = cases / (width * (control - (n*0.5))), # rate within category i.e. cases/exposed person years. Again, exposed person years are the py up to and including that category
               w_risk = cumsum(rate2 * width), # think of as a cumulative rate??
               cum_risk = 1-exp(-w_risk),
               na = cases / (cum_cases-(n*0.5)),
               na_haz = cumsum(na))) %>% # Calculation of cumulative risk up to dose-time
  bind_rows()

# Combine two data frames 'risk_py' and 'risk_cat' using the 'bind_rows' function from the 'dplyr' package.
risk2 <- bind_rows(risk_py, risk_cat)

# Add a new column "pop" based on "id"
risk2 <- risk2 %>%
  mutate(pop = ifelse(id %in% c("usf", "cpo"), "non-mine", "mine"))


# Create a new column 'risk_diff' in the 'risk2' data frame that represents the difference between 'cum_risk' and 'rep_cum_risk'.
risk2$risk_diff <- risk2$cum_risk - risk2$rep_cum_risk

# Calculate summary statistics for the 'risk_diff' column using the 'describe' function.
describe(risk2$risk_diff)

# Plot that compared to the reported and calculated risks
raw_risk2 <- subset(risk2, !(id %in% c("det", "usf", "mpo", "scc")))

# Convert to long format
df_long <- raw_risk2 %>%
  pivot_longer(cols = c(cum_risk, rep_cum_risk),
               names_to = "variable",
               values_to = "cum_risk2")


# Investigate a little deeper into reported vs calculated risks

# Calculate SAS formula 

rep_calc <- subset(risk2, !(id %in% c("det")))

rep_calc <- rep_calc %>% 
  group_by(id) %>% 
  mutate(q = cases / (control),
         p = 1 - q,
         St = 1) %>%
  bind_rows()

# Create a new column St and initialize it with 1
rep_calc$St <- 1

# Calculate St using a for loop within each group
for (group_id in unique(rep_calc$id)) {
  group_data <- rep_calc %>% filter(id == group_id)
  St <- group_data$St
  p <- group_data$p
  
  for (i in 2:length(St)) {
    St[i] <- St[i - 1] * p[i - 1]
  }
  
  rep_calc$St[rep_calc$id == group_id] <- St
}

# Calculate S(t) and add it to the data frame
rep_calc$strisk <- 1 - rep_calc$St

# Select relevant variables
select <- rep_calc %>% select(id, author, d_low, dose, rep_cum_risk,cum_risk, strisk, pop)
glimpse(rep_calc)

# Example of S(t) calculation for supplement
example <- rep_calc %>% select(id, author, dose, cases, control, St, p, q, strisk, rep_cum_risk) %>% filter(author == "Chen et al. (tungsten)")

# Example of Steenland for supplement 
example_ht <- risk_py %>% select(id, author, dose, cases, control, width, rate2, w_risk, cum_risk, rep_cum_risk) %>% filter(author == "Steenland and Brown")

# Write to excel
write_xlsx(example, "example_st.xlsx")
write_xlsx(example_ht, "example_ht.xlsx")

# Convert to long format
sel_long <- select %>%
  pivot_longer(cols = c(rep_cum_risk, cum_risk, strisk),
               names_to = "var",
               values_to = "risk")
glimpse(sel_long)

# If we wanted to move the cumulative risk estimate to the left hand side of the dose range for studies in which SAS formula is used 
# sel_long <- sel_long %>%
  # mutate(dose_2 = ifelse(id %in% c("cpo", "ctu", "cti", "sag_2") & var == "strisk", d_low, dose))

# Plot of CALCULATED cumulative risks - comparison betweeen studies - with the STEENLAND FORMULA (cum_risk)
raw_risk <- subset(sel_long, (var %in% c("cum_risk"))) # Remove hughes et al - not possible to calculate

glimpse(raw_risk)

#geom_label(data = summarized_data[!(summarized_data$id %in% c("usg", "ctu")), ], aes(x = max_dose, y = max_cum_risk, label = id), alpha = 0.8) +
#geom_label(data  = summarized_data[summarized_data$id == "ctu", ], aes(x = max_dose, y = max_cum_risk, label = id), nudge_x = -0.5, alpha = 0.8) +
#geom_label(data = summarized_data[summarized_data$id == "usg", ], aes(x = max_dose, y = max_cum_risk, label = id),alpha = 0.8, nudge_x = 0) 

raw_risk <- raw_risk %>% 
  dplyr::filter(!id %in% c("sag_2"))

p1 <- ggplot(raw_risk, aes(x=dose, y=risk, color=id, linetype = pop)) +
  geom_line(linewidth = 0.8) +
  theme_pubr() +
  labs(x = expression(RCS~mg/m^3~"-years"), y = "Cumulative risk (%)", colour = "Author", linetype = "Population") + 
  scale_color_manual(
    values = c("ush" = "#ff9040", "sag" = "#a8ae3e", "usg" = "#6971d7", "scc" = "#3c8350","cti" = "#45c097", "mpo" = "#be0034", 
               "ctu" = "#5b3787", "cpo" = "#b64773", "usf" = "#bc7e36"),
    breaks = c("ush", "sag", "usg","scc",  "cti", "mpo", "ctu", "cpo", "usf"),
    labels = c("Kreiss and Zhen", "Hnizdo and Sluis-Cremer", "Steenland and Brown", "Miller et al.","Chen et al. (tin)", "Liu et al.", 
                "Chen et al. (tungsten)","Chen et al. (pottery)", "Rosenman et al.")) +
scale_linetype_manual(
    values = c("mine" = "solid", "non-mine" = "dashed"),
    labels = c("mine" = "Mine", "non-mine" = "Non-mine")
  ) 

p1

# Create Loess smoothed cum_risk plot 
ggplot(risk2, aes(x = dose, y = (cum_risk*100), color = id)) +
  geom_point() +  
  theme_pubr() +
  geom_smooth(se = FALSE, span = 0.75) +  # Adjust the span value (default is 0.75)
  labs(x = expression(RCS~mg/m^3~"-years"), y = "Silicosis (%)", colour = "Study", linetype = "Population") + 
  scale_color_manual(
    values = c("sag" = "#a8ae3e","usg" = "#6971d7", "scc" = "#3c8350","cti" = "#45c097", "mpo" = "#be0034", 
               "ctu" = "#5b3787", "cpo" = "#b64773", "usf" = "#bc7e36", "ush" = "#ff9040"),
    breaks = c("sag", "usg","scc",  "cti", "mpo", "ctu", "cpo", "usf"),
    labels = c("Hnizdo and Sluis-Cremer", "Steenland and Brown", "Miller et al.","Chen et al. (tin)", "Liu et al.", 
               "Chen et al. (tungsten)","Chen et al. (pottery)", "Rosenman et al."))

# Create point estimates from Loess line OF CUMULATIVE RISK
risk4 <- subset(risk2, !(id %in% c("det"))) # Remove hughes et al - not possible to calculate

# Get unique id values
unique_ids <- unique(risk4$id)

# Initialize an empty data frame to store results
result_df <- data.frame()

# Loop through each unique id
for (id_value in unique_ids) {
  # Create a loess model for the subset of data with the current id value
  subset_data <- risk4 %>% filter(id == id_value)
  loess_model <- stats::loess(cum_risk ~ dose, data = subset_data, span = 0.75)
  
  # Specify the dose values where you want to extract predictions
  dose_values <- c((0.025*40), (0.05*40), (0.1*40), (0.25*40), (0.5*40))
  
  # Create a data frame with the specified dose values
  data_to_predict <- data.frame(dose = dose_values)
  
  # Use the loess model to predict na_haz values at the specified dose values
  predicted_values <- predict(loess_model, newdata = data_to_predict)
  
  # Multiply the predicted cum_risk values by 100
  predicted_values <- predicted_values * 100
  
  # Create a data frame with the results for the current id
  id_result <- data.frame(id = id_value, dose = dose_values, cum_risk = predicted_values)
  
  # Append the current id's results to the overall result data frame
  result_df <- bind_rows(result_df, id_result)
}

# Reshape the result_df to have dose values horizontally and study names vertically
wide_result_df <- result_df %>%
  pivot_wider(names_from = dose, values_from = cum_risk)

wide_result_df <- wide_result_df %>%
  mutate(pop = ifelse(id %in% c("usf", "cpo"), "non-mine", "mine"))

# Print the reshaped data frame
print(wide_result_df)
write_xlsx(wide_result_df, "rr_comparison.xlsx")

########################################
# Comparison of fitted cumulative risks#
########################################

# Now try and draw all of the cumulative risk fitted lines on one graph

chen2 <- read.csv("chen2.csv", header = TRUE)

# Create df for fitted distributions of different variables

# Hughes hi risk
df <- data.frame(group = "deh - Hughes et al. (high)", 
                 rcs = seq(0, 20, length.out = 200))
df$risk <- 1 - (1 / (1 + exp((-3.250 + 0.731) / 0.597) * (df$rcs^(1 / 0.597))))

# Hughes lo risk
df1 <- data.frame(group = "del - Hughes et al. (low)", 
                  rcs = seq(0, 20, length.out = 200))
df1$risk <- 1 - (1 / (1 + exp((-3.250) / 0.597) * (df1$rcs^(1 / 0.597))))

# Miller 
df2 <- data.frame(group = "scc - Miller et al.", 
                  rcs = seq(0, 20, length.out = 200))
df2$risk <- exp(-4.32 + 0.416*df2$rcs) / (1+exp(-4.32 + 0.416*df2$rcs))

# Hnizdo
df3 <- data.frame(group = "sag - Hnizdo and Sluis-Cremer",
                  rcs = seq(0, 20, length.out = 200)) %>%
  mutate(cde = rcs * (10/3),
         risk = 1 - (1 / (1 + exp((-2.439)/0.2199) * (cde^(1/0.2199)))))

# Combine models
models <- bind_rows(chen2, df, df1, df3)
str(models)

# Labels for graphs 
models <- models %>%
  mutate( id = case_when(
    group == "cpo - Chen et al. (Pottery)" ~ "cpo",
    group == "cti - Chen et al. (Tin)" ~ "cti",
    group == "ctu - Chen et al. (Tungsten)" ~ "ctu",
    group == "deh - Hughes et al. (high)" ~ "deh",
    group == "del - Hughes et al. (low)" ~ "del",
    group == "sag - Hnizdo and Sluis-Cremer" ~ "sag"))
table(models$group)

# Combine all the dfs

plot_lab <- models %>%
  group_by(id) %>%
  summarise(max_dose = quantile(rcs,0.60), max_cum_risk = quantile(risk,0.60))

# Labels mining and non-mining
models <- models %>%
  mutate(pop = ifelse(id %in% c("cpo"), "non-mine", "mine"))

p2 <- ggplot(models, aes(x = rcs, y = risk, color = id, linetype = pop)) +
  geom_line(linewidth = 0.8) +
  theme_pubr() + 
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Cumulative risk (%)", color = "Author", linetype = "Population") +
  scale_color_manual(
    values = c("sag" = "#a8ae3e", "cti" = "#45c097", "ctu" = "#5b3787", "deh" = "#9d32a8", 
               "cpo" = "#b64773", "del" = "#464a4d"),  
    breaks = c("sag", "cti", "ctu", "deh", "cpo", "del"),
    labels = c("Hnizdo and Sluis-Cremer", "Chen et al. (tin)",
               "Chen et al. (tungsten)", "Hughes et al. (high)", "Chen et al. (pottery)",
               "Hughes et al. (low)")
  ) +
  scale_linetype_manual(
    values = c("mine" = "solid", "non-mine" = "dashed"),
    labels = c("mine" = "Mine", "non-mine" = "Non-mine")
  )
p2

ggsave("p2.png", p2, dpi = 600, width = 14, height = 8)


# Create a shared plot
fig_2 <- ggpubr::ggarrange(p1 + theme(legend.position="none"), 
                  p2 + theme(legend.position="none"), # list of plots
                  labels = "AUTO", # labels
                  align = "hv", # Align them both, horizontal and vertical
                  ncol = 2, 
                  common.legend = T)  + 
  theme(legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16))  # Adjust the font size as needed
# number of rows
fig_2
ggsave("fig_2.png", fig_2, dpi = 600, width = 14, height = 8)


# Predicted formula values from the fitted curves

# Get unique id values
chen3 <- chen2 %>%
  mutate( id = case_when(
    group == "cpo - Chen et al. (Pottery)" ~ "cpo_fit",
    group == "cti - Chen et al. (Tin)" ~ "cti_fit",
    group == "ctu - Chen et al. (Tungsten)" ~ "ctu_fit"))
ids <- unique(chen3$id)

# Initialize an empty data frame to store results
fitted_df <- data.frame()

# Loop through each unique id
for (id_value in ids) {
  # Create a loess model for the subset of data with the current id value
  sub_data <- chen3 %>% filter(id == id_value)
  loess_model <- stats::loess(risk ~ rcs, data = sub_data, span = 0.7)
  
  # Specify the dose values where you want to extract predictions
  dose_values <- c((0.025*40), (0.05*40), (0.1*40), (0.25*40), (0.5*40))
  # values for OSHA (0.05*45) (0.1*45) (0.25*45) (0.5*45)
  
  # Create a data frame with the specified dose values
  df_fit <- data.frame(rcs = dose_values)
  
  # Use the loess model to predict risk values at the specified dose values
  pred_fit <- predict(loess_model, newdata = df_fit)
  
  # Multiply the predicted risk values by 100
  val_fit <- pred_fit * 100
  
  # Create a data frame with the results for the current id
  id_fit <- data.frame(id = id_value, rcs = dose_values, risk = val_fit)
  
  # Append the current id's results to the overall result data frame
  fitted_df <- bind_rows(fitted_df, id_fit)
}

# Remove duplicates based on id and rcs
fitted_df <- fitted_df %>%
  distinct(id, rcs, .keep_all = TRUE)

# Reshape the result_df to have dose values horizontally and study names vertically
fit_wide <- fitted_df %>%
  pivot_wider(names_from = rcs, values_from = risk)

# Print the reshaped data frame
print(fit_wide)
print(wide_result_df)

# For other studies

# Hughes hi risk
deh_fit <- data.frame(id = "deh_fit", 
                 rcs = dose_values)
deh_fit$risk <- 1 - (1 / (1 + exp((-3.250 + 0.731) / 0.597) * (deh_fit$rcs^(1 / 0.597))))
deh_fit$risk <- deh_fit$risk *100

# Hughes lo risk
del_fit <- data.frame(id = "del_fit", 
                  rcs = dose_values)
del_fit$risk <- 1 - (1 / (1 + exp((-3.250) / 0.597) * (del_fit$rcs^(1 / 0.597))))
del_fit$risk <- del_fit$risk *100

# Hnizdo
sag_fit <- data.frame(id = "sag_fit",
                  rcs = dose_values) %>%
  mutate(cde = rcs * (10/3),
         risk = 1 - (1 / (1 + exp((-2.439)/0.2199) * (cde^(1/0.2199)))))
sag_fit$risk <- sag_fit$risk * 100

# Combine models (assuming deh_fit, deh_fit, sag_fit are already defined)
fitted_df2 <- bind_rows(del_fit, deh_fit, sag_fit) 
fitted_df2
# Remove the "cde" column from fitted_df2
fitted_df2$cde <- NULL

# Reshape the resulting data frame fitted_df2
fit_wide2 <- fitted_df2 %>%
  pivot_wider(names_from = rcs, values_from = risk)

# Chen cohorts
fit_wide
# Other fitted values 
fit_wide2

# Non-parametric estimates
wide_result_df

# Combine data frames into a single data frame
df_sum <- bind_rows(wide_result_df, fit_wide, fit_wide2)

# Define the order of levels for the 'id' variable
order <- c("cti", "cti_fit", "ctu", "ctu_fit", "sag", "sag_fit", "sag_2", "deh_fit", "del_fit", "mpo", "scc", "usg", "usf", "cpo",
           "cpo_fit", "ush") 

# Mutate and arrange the data frame based on the defined order of levels
df_sum2 <- df_sum %>%
  mutate(id = factor(id, levels = order)) %>% 
  arrange(id)

# For chen studies, compare raw, fitted, and reported data

# Rename columns in df_long
sel_long$rcs <- sel_long$dose
sel_long$variable <- sel_long$var

# Select relevant columns and create a new data frame 'comp'
comp <- sel_long %>%
  dplyr::select(id, risk, rcs, variable, author)

# Create a new column 'id' with the value 'sag' in df3
df3$id <- 'sag'

# Append df3 to chen3 data frame
chen4 <- bind_rows(chen3, df3)

# Set the 'variable' column to "fitted" in chen3
chen4$variable <- "fitted"

# Select relevant columns and create a new data frame 'comp2' from chen3
comp2 <- chen4 %>%
  dplyr::select(id, risk, rcs, variable)

# Rename 'id' values in comp2 based on specified conditions
comp2 <- comp2 %>%
  mutate(id = case_when(
    id == "cpo_fit" ~ "cpo",
    id == "ctu_fit" ~ "ctu",
    id == "cti_fit" ~ "cti",
    TRUE ~ id  # Keep other values as they are
  ))

# Rename 'author' values in comp2 based on specified conditions=
comp2 <- comp2 %>%
  mutate( author = case_when(
    id == "cpo" ~ "Chen et al. (pottery)",
    id == "cti" ~ "Chen et al. (tin)",
    id == "ctu" ~ "Chen et al. (tungsten)",
    id == "sag" ~ "Hnizdo and Sluis-Cremer"))

# Combine data from 'comp' and 'comp2' into 'compare' data frame
compare <- rbind(comp, comp2)

# Create a line plot comparing risk vs. rcs with color and linetype differentiation
ggplot(compare, aes(x = rcs, y = risk, color = id, linetype = variable)) + 
  geom_line(linewidth = 1) +  # Line plot
  scale_color_manual(
    values = c("sag" = "#a8ae3e","usg" = "#6971d7", "cti" = "#45c097", "ctu" = "#5b3787", "cpo" = "#b64773"),
    breaks = c("sag", "usg", "cti", "ctu", "cpo"),
    labels = c("Hnizdo and Sluis-Cremer", "Steenland and Brown", "Chen et al. (tin)",
               "Chen et al. (tungsten)","Chen et al. (pottery)")) +
  scale_linetype_manual(
    values = c("cum_risk" = "solid", "rep_cum_risk" = "dashed", "fitted" = "dotted"),
    labels = c("cum_risk" = "Calculated", "rep_cum_risk" = "Reported", "fitted" = "Fitted")
  ) +
  theme_pubr() +  
  labs(x = expression(RCS~mg/m^3~"-years"), y = "Cumulative risk (%)", colour = "Study", linetype = "Risk")  # Add labels


# Make a nice facet plot
compare2 <- compare %>% 
  dplyr::filter(!id %in% c("sag_2"))
p3 <- ggplot(compare2, aes(x = rcs, y = risk, color = id, linetype = variable)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(
    values = c("sag" = "#a8ae3e","sag_2" = "#808080", "usg" = "#6971d7", "scc" = "#3c8350","cti" = "#45c097", "mpo" = "#be0034", 
               "ctu" = "#5b3787", "cpo" = "#b64773", "usf" = "#bc7e36", "ush" = "#ff9040"),
    guide = FALSE) +  # Set guide to FALSE to remove the color legend
  scale_linetype_manual(
    values = c("rep_cum_risk" = "solid", "cum_risk" = "dashed",  "strisk" = "dotted","fitted" = "dotdash"),
    labels = c("rep_cum_risk" = "Reported", "cum_risk" = "Steenland", "strisk" = "SAS","fitted" = "Fitted")) +
  theme_pubr() +
  labs(x = expression(RCS~mg/m^3~"-years"), y = "Cumulative risk", linetype = "Linetype") + 
  facet_wrap(~ author)
p3

# Create a DF to compare all other estaimtes of cumulative to St risk
# Tricky, as need to create a df of all the same risk values

# Pivot the 'comp' data frame to wider format using 'risk' values as columns
comp3 <- pivot_wider(comp, names_from = variable, values_from = risk)

# Filter out rows where 'id' is "sag_2"
comp3 <- comp3 %>% 
  dplyr::filter(!id %in% c("sag_2"))

# Initialize an empty data frame to store results
fitted_df <- data.frame()

# Create a loess model for the subset of data with the current id value
cpo_data <- chen3 %>% filter(id == "cpo_fit")
cpo_model <- stats::loess(risk ~ rcs, data = cpo_data, span = 0.7)
ctu_data <- chen3 %>% filter(id == "ctu_fit")
ctu_model <- stats::loess(risk ~ rcs, data = ctu_data, span = 0.7)
cti_data <- chen3 %>% filter(id == "cti_fit")
cti_model <- stats::loess(risk ~ rcs, data = cti_data, span = 0.7)

# Specify the dose values where you want to extract predictions
# Use mutate along with case_when to calculate 'fitted' values based on 'id'
comp3 <- comp3 %>%
  mutate(fitted = case_when(
    id == "cpo" ~ predict(cpo_model, newdata = data.frame(rcs = rcs)),
    id == "cti" ~ predict(cti_model, newdata = data.frame(rcs = rcs)),
    id == "ctu" ~ predict(ctu_model, newdata = data.frame(rcs = rcs)),
    id == "sag" ~ (1 - (1 / (1 + exp((-2.439)/0.2199) * ((rcs * (10/3))^(1/0.2199))))),
    TRUE ~ NA_real_  # Set other cases to NA
  ))

# Calculate new columns 'cum_fit' and 'cum_st' - the difference between fitted and Steenland risk
comp3$cum_fit <- comp3$fitted - comp3$cum_risk
comp3$cum_st <- comp3$strisk - comp3$cum_risk

# Summary of fitted values vs Steenland 
summary(comp3$cum_fit)

# Summary of SAS vs Steenland values 
summary(comp3$cum_st)

# Convert to long format for facet plot of difference in risk, compared to SAS formula
comp4 <- comp3 %>%
  pivot_longer(cols = c("cum_fit", "cum_st"), names_to = "diff_type", values_to = "risk_diff")
table(comp4$id)

# Make a nice facet plot of differences between Steenland and other formulas
p4 <- ggplot(comp4, aes(x = rcs, y = risk_diff, color = id, linetype = diff_type)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(
    values = c("sag" = "#a8ae3e","sag_2" = "#808080", "usg" = "#6971d7", "scc" = "#3c8350","cti" = "#45c097", "mpo" = "#be0034", 
               "ctu" = "#5b3787", "cpo" = "#b64773", "usf" = "#bc7e36", "ush" = "#ff9040"),
    guide = FALSE) +  # Set guide to FALSE to remove the color legend
  scale_linetype_manual(
    values = c("cum_fit" = "solid", "cum_st" = "dashed"),
    labels = c("cum_fit" = "Steenland (ref) vs Fitted", "cum_st" = "Steenland (ref) vs SAS ")) +
  theme_pubr() +
  labs(x = expression(RCS~mg/m^3~"-years"), y = "Cumulative risk difference", linetype = "Linetype") + 
  facet_wrap(~ author)
  
# Create a shared plot
fig_3 <- ggpubr::ggarrange(p3, p4, # list of plots
                  labels = "AUTO", # labels
                  align = "hv", # Align them both, horizontal and vertical
                  ncol = 2)  # number of rows

ggsave("fig_3.png", fig_3, width = 10, height = 5.5, dpi = 600)

##############################
# Dose response meta-analysis#
##############################

# Split data into cross sectional analysis and longitudinal analysis
risk_cat <- risk_drma %>% subset(., type == "cc")
table(risk_cat$author)

risk_py <- risk_drma %>% subset(., type == "ir")
table(risk_py$author)

# For cumulative incidence
risk_cat <- risk_cat %>% 
  group_split(id) %>% 
  map(~ mutate(.x, 
               width = d_hi - d_low, # width of the category
               # rate = cases / (width * (cum_cases-(n*0.5))), # rate within category i.e. cases/exposed person years, where exposed years are CUMULATIVE
               # w_risk = cumsum(rate * width), # Calculation of cumulative hazard rate
               # cum_risk = (cases/n), # pg Kirkwood 146-7
               or = cases / (n-cases))) %>% 
  bind_rows()


# Calculate odds ratios (rate ratio for dosresmeta) for categorical data
risk_cat <- risk_cat %>% 
  group_by(id) %>% 
  mutate(rr = or / or[1],
         logrr = log(rr), 
         selogrr = sqrt(1/cases[1] + (1/n[1]) + (1/cases + 1/n))) %>% # pg 164 Kirkwood and Sterne
  mutate(selogrr = ifelse(row_number() == 1, NA, selogrr)) %>% 
  ungroup()

# For person years data
risk_py <- risk_py %>% 
  group_split(id) %>% 
  map(~ mutate(.x, 
               width = d_hi - d_low, 
               rate = cases / n)) %>% # incidence rate within a band
               # rate2 = cases / (width * (control - (cases*0.5))), # rate within category i.e. cases/exposed person years. Again, exposed person years are the py up to and including that category
               # w_risk = cumsum(rate2 * width), # think of as a cumulative rate??
               # cum_risk = 1-exp(-w_risk))) 
  bind_rows()

# Calculate rate ratio, compared to baseline for py_data
risk_py <- risk_py %>% 
  group_by(id) %>% 
  mutate(rr = rate / rate[1], 
         logrr = log(rr), 
         selogrr = sqrt((1/cases[1]) + (1/cases))) %>% # p 242 Kirwood and Sterne 
  mutate(selogrr = ifelse(row_number() == 1, NA, selogrr)) %>% 
  ungroup()

risk2 <- bind_rows(risk_cat, risk_py)

# Create df to work with
risk3 <- risk2 %>%  dplyr::select(id, rr, logrr, dose, type, selogrr, cases, n)

# Separate mine and non-mine groups
mine <- risk3 %>% 
  dplyr::filter(!id %in% c("cpo", "usf"))
table(mine$id)

nonmine <- risk3 %>% 
  dplyr::filter(id %in% c("cpo", "usf"))
table(nonmine$id)


##############
# Miner DRMA #
#############

# Plot log RR's against  dose
sup_1 <- ggplot(mine, aes(dose, logrr, size = selogrr, color = id)) +
  geom_point(shape = 19, alpha = 0.8) +
  scale_size_area(max_size = 9) +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Log relative risk", color = "Author", size = "SE") +
  theme_pubr() + 
  scale_color_manual(
    values = c("sag_2" = "#a8ae3e", "usg" = "#6971d7", "scc" = "#3c8350","cti" = "#45c097", "mpo" = "#be0034", 
               "ctu" = "#5b3787","ush" = "#ff9040", "det" = "#9d32a8"),
    breaks = c("sag_2", "usg","scc",  "cti", "mpo", "ctu", "ush", "det"),
    labels = c("Hnizdo and Sluis-Cremer", "Steenland and Brown", "Miller et al.","Chen et al. (tin)", "Liu et al.", 
               "Chen et al. (tungsten)", "Kreiss and Zhen", "Hughes et al."))

ggsave("sup_1.png", sup_1, width = 14, height = 7, dpi = 600)

# Linear model
# Fit a linear model using dosresmeta function
lin_bin <- dosresmeta(formula = logrr ~ dose, id = id, type = type, 
                      se = selogrr, cases = cases, n = n, data = mine)

# Display a summary of the linear model
summary(lin_bin)

# RCS model
# Define knots for the RCS model based on quantiles of the 'dose' variable
knots <- quantile(mine$dose, c(.05, .35, .65, .95))
knots

# Fit an RCS model using dosresmeta function with specified knots
spl_bin_mine <- dosresmeta(formula = logrr ~ rcs(dose, knots), type = type, 
                           id = id, se = selogrr, cases = cases, n = n, data = mine)

# Quadratic model
# Fit a quadratic model using dosresmeta function
quad_bin <- dosresmeta(formula = logrr ~ dose + (dose^2), type = type, 
                       id = id, se = selogrr, cases = cases, n = n, data = mine)

# Display summaries for the quadratic and RCS models
summary(quad_bin)
summary(spl_bin_mine)

# Wald test for the RCS model
waldtest(b = coef(spl_bin_mine), Sigma = vcov(spl_bin_mine), Terms = 2:3)

# Set a refence cum RCS; in this case 4 mg/mg3-years i.e. 0.1 for 40 years
xref <- 4

# Predict risk between 0 and 20 range
pred_mine <- data.frame(dose = c(xref, seq(0, 20, 0.001))) %>%
  predict(spl_bin_mine, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2)
head(pred_mine)

# Plot that risk
p5 <- ggplot(pred_mine, aes(dose, pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_ribbon(alpha = .05) +
  geom_smooth(color = "black", linewidth = 0.5) + 
  scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50, 100), limits = c(0.1,100)) +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Relative risk of silicosis (log-scale)") +
  theme_pubr() +
  theme(text = element_text(size = 18)) +
  ggtitle("Miners")
p5

# Table by increases of 1 mg/m3
xref <- 4
pred_mine_tab <- data.frame(dose = c(xref, c(1,2,4,10))) %>%
  predict(spl_bin_mine, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2) %>% arrange(dose)
pred_mine_tab

##################
# Non-miner DRMA #
#################

# Plot log RR's against  dose
sup_2 <- ggplot(nonmine, aes(dose, logrr, size = selogrr, color = id)) +
  geom_point(shape = 19, alpha = 0.8) +
  scale_size_area(max_size = 9) +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Log relative risk", color = "Author", size = "SE") +
  theme_pubr() + 
  scale_color_manual(
      values = c("cpo" = "#b64773", "usf" = "#bc7e36"),
      breaks = c("cpo", "usf"),
      labels = c("Chen et al. (pottery)", "Steenland and Brown"))

ggsave("sup_2.png", sup_2, width = 14, height = 7, dpi = 600)

# RCS model
# Define knots for the RCS model based on quantiles of the 'dose' variable
knots <- quantile(nonmine$dose, c(.1, .5, .9)) 
knots

# Fit an RCS model using dosresmeta function with specified knots
spl_bin_non_mine <- dosresmeta(formula = logrr ~ rcs(dose, knots), type = type, 
                               id = id, se = selogrr, cases = cases, n = n, data = nonmine)

# Display a summary of the RCS model
summary(spl_bin_non_mine)
xref <- 4

# Predictions df
pred_non_mine <- data.frame(dose = c(xref, seq(0, 20, 0.1))) %>%
  predict(spl_bin_non_mine, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`) %>% 
  unique() %>% round(2)
head(pred_non_mine)

# Plot of risk
p6 <- ggplot(pred_non_mine, aes(dose, pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_ribbon(alpha = .05) +
  geom_smooth(color = "black", linewidth = 0.5) + 
  scale_y_continuous(trans = "log", breaks = c(0,1,5,10,50,100), limits = c(0.1,100)) +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Relative risk of silicosis (log-scale)") +
  theme_pubr() +
  theme(text = element_text(size = 18)) +
  ggtitle("Non-miners")
p6
xref <- 4

# Table by increase of 1 mg/m3
pred_non_mine_tab <- data.frame(dose = c(xref, c(1,2,4,10))) %>%
  predict(spl_bin_non_mine, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`) %>% 
  unique() %>% round(2)  %>% arrange(dose)
pred_non_mine_tab

## Compare mine and non-mine risks 

# Set a refence cum RCS; in this case 2 mg/mg3-years i.e. 0.1 for 20 years
xref <- 0.5

# Predict risk mine
pred_mine <- data.frame(dose = c(xref, seq(0.495, 21.995, 0.01)))%>%
  predict(spl_bin_mine, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2)
head(pred_mine)
pred_mine

# Predict risk non-mine
pred_non_mine <- data.frame(dose = c(xref, seq(0.495, 21.995, 0.01))) %>%
  predict(spl_bin_non_mine, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`) %>% 
  unique() %>% round(2)
head(pred_non_mine)

# Assign a label indicating "Miner morbidity" to the predictions for mine data
pred_mine$miner <- "Miner morbidity"
# Assign a label indicating "Non-miner morbidity" to the predictions for non-mine data
pred_non_mine$miner <- "Non-miner morbidity"

# Combine predictions for mine and non-mine data into a single data frame
mine_comp <- bind_rows(pred_mine, pred_non_mine)

# Create a ggplot for visualizing and comparing predicted relative risks
ggplot(mine_comp, aes(dose, pred, color = miner, ymin = ci.lb, ymax = ci.ub)) +
  geom_ribbon(alpha = 0.05, size = 0) +
  geom_smooth(linewidth = 0.5) + 
  scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50, 100, 500)) +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Relative Risk") + 
  theme_bw()

# Compare study to mortality study Mannetje et al 
# Create Mannetje df
mann <- data.frame(dose = c(0.495, 1.482, 2.423, 3.65, 5.725, 8.35, 11.395, 14.552, 21.995),
                   pred = c(1.00, 3.39, 6.22, 9.40, 13.69, 22.64, 23.97, 40.25, 25.11),
                   miner = "Mortality")

# Select relevant columns from mine_comp data frame
mine_comp <- select(mine_comp, dose, pred, miner)

# Combine predictions for mine and Mannetje et al into a single data frame
mine_comp <- bind_rows(mine_comp, mann)

# Define the order of levels for the 'miner' variable
mine_comp$miner <- factor(mine_comp$miner, levels = c("Miner morbidity", "Non-miner morbidity", "Mortality"))

# The plot is designed to compare relative risks for different groups at various dose levels
sup_3 <- ggplot(mine_comp, aes(dose, pred, color = miner)) +
  geom_line(size = 1) +
  scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50, 100, 500)) +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Relative Risk (log-scale)", color = "Group") + 
  theme_pubr()

ggsave("sup_3.png", sup_3, width = 8, height = 5, dpi = 600)

# Plot of the predicted relative risks for miners and non-miners
fig_4 <- ggpubr::ggarrange(p5, p6, # list of plots
                  labels = "AUTO", # labels
                  align = "hv", # Align them both, horizontal and vertical
                  ncol = 2)  # number of rows

ggsave("fig_4.png", fig_4, width = 14, height = 8, dpi = 600)

#######################
# Sensitivity analysis #
#######################

# No SCC as 2/1 outcome and USG as death certificate AND only COHORT studies 
mine_ilo <- risk3 %>% 
  dplyr::filter(!id %in% c("scc", "usg", "cpo", "usf"))
table(mine_ilo$id)

# Define knots for the RCS model based on quantiles of the 'dose' variable
knots <- quantile(mine_ilo$dose, c(.05, .35, .65, .95))
knots

# Fit an RCS model using dosresmeta function with specified knots for sensitivity analysis
spl_bin_ilo <- dosresmeta(formula = logrr ~ rcs(dose, knots), type = type, 
                          id = id, se = selogrr, cases = cases, n = n, data = mine_ilo)

# Display a summary of the RCS model for sensitivity analysis
summary(spl_bin_ilo)

# Set a refence cum RCS; in this case 2 mg/mg3-years i.e. 0.1 for 20 years
xref <- 4

# Predict risk between 0 and 20 range
pred <- data.frame(dose = c(xref, seq(0, 20, 0.001))) %>%
  predict(spl_bin_ilo, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2)
head(pred)

# Plot that risk
ggplot(pred, aes(dose, pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_ribbon(alpha = .05) +
  geom_smooth(color = "black", linewidth = 0.5) + 
  scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50, 100)) +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Relative risk of silicosis (log-scale)") +
  theme_bw() +
  theme(text = element_text(size = 18))

# Table by increases of 1 mg/m3
xref <- 4
pred_ilo <- data.frame(dose = c(xref, c(1,2,4,10))) %>%
  predict(spl_bin_ilo, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2) %>% arrange(dose)
pred_ilo 

# Low (<0.1) and high intensity miners
# Filter low intensity miners (usg, sag_2, ush)
low <- mine %>% 
  dplyr::filter(id %in% c("usg", "sag_2", "ush"))

# Filter high intensity miners (cti, ctu, scc, det)
high <- mine %>% 
  dplyr::filter(id %in% c("cti", "ctu", "scc", "det"))

# RCS model for low intensity miners
# Define knots for the RCS model based on quantiles of the 'dose' variable for low intensity miners
knots <- quantile(low$dose, c(.05, .35, 0.65, 0.95))
spl_bin_low <- dosresmeta(formula = logrr ~ rcs(dose, knots), type = type, 
                          id = id, se = selogrr, cases = cases, n = n, data = low)
summary(spl_bin_low)

# RCS model for high intensity miners
# Define knots for the RCS model based on quantiles of the 'dose' variable for high intensity miners
knots <- quantile(high$dose, c(.05, .35, 0.65, 0.95))
spl_bin_high <- dosresmeta(formula = logrr ~ rcs(dose, knots), type = type, 
                           id = id, se = selogrr, cases = cases, n = n, data = high)
summary(spl_bin_high)

# Plot
xref <- 0

# Predict risk between 0 and 20 range
pred_low <- data.frame(dose = c(xref, seq(0, 20, 0.001))) %>%
  predict(spl_bin_low, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`) %>% 
  unique() %>% round(2)

xref <- 4
pred_l_tab <- data.frame(dose = c(xref, c(1,2,4,10))) %>%
  predict(spl_bin_low, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`) %>% 
  unique() %>% round(2)
pred_l_tab

# Predict risk between 0 and 20 range
xref <- 0
pred_high <- data.frame(dose = c(xref, seq(0, 20, 0.001))) %>%
  predict(spl_bin_high, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`) %>% 
  unique() %>% round(2)

xref <- 4
pred_h_tab <- data.frame(dose = c(xref, c(1,2,4,10))) %>%
  predict(spl_bin_high, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`) %>% 
  unique() %>% round(2)
pred_h_tab 

# Create an 'intensity' variable to distinguish between low and high intensity miners
pred_low$intensity <- "low"
pred_high$intensity <- "high"

# Combine predictions for low and high intensity miners into a single data frame
intensity <- bind_rows(pred_low, pred_high)

# Create a ggplot for visualizing and comparing predicted relative risks by intensity
ggplot(intensity, aes(dose, pred, color = intensity, ymin = ci.lb, ymax = ci.ub)) +
  geom_ribbon(alpha = 0.05, size = 0) +  # Add ribbons for confidence intervals
  geom_smooth(linewidth = 0.5) +  # Add smoothed lines
  scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50, 100, 500)) +  
  # Apply log transformation to the y-axis with specified breaks
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Relative Risk") + 
  theme_bw()  # Apply a white and black theme for aesthetics

# No det or mpo as some non-miners 

mine_x <- mine %>% 
  dplyr::filter(!id %in% c("det", "mpo", "ush"))
table(mine_x$id)

knots <- quantile(mine_x$dose, c(.05, .35, .65, .95))
knots

spl_bin_x <- dosresmeta(formula = logrr ~ rcs(dose, knots), type = type, 
                          id = id, se = selogrr, cases = cases, n = n, data = mine_x)
summary(spl_bin_x)

# Set a refence cum RCS; in this case 2 mg/mg3-years i.e. 0.1 for 20 years
xref <- 4

# Predict risk between 0 and 20 range
pred_x <- data.frame(dose = c(xref, seq(0, 20, 0.001))) %>%
  predict(spl_bin_x, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2)
head(pred_x)

# Plot that risk
ggplot(pred_x, aes(dose, pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_ribbon(alpha = .05) +
  geom_smooth(color = "black", linewidth = 0.5) + 
  scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50, 100)) +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Relative risk of silicosis (log-scale)") +
  theme_bw() +
  theme(text = element_text(size = 18))

# Table by increases of 1 mg/m3
xref <- 4
pred_x2 <- data.frame(dose = c(xref, c(1,2,4,10))) %>%
  predict(spl_bin_x, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2) %>% arrange(dose)
pred_x2 


# Final sensitivity analysis; doses <8 mg/m3-years

mine_low <- subset(mine_x, dose <= 8)

knots <- quantile(mine_low$dose, c(.05, .35, .65, .95))
knots
spl_bin_mine_low <- dosresmeta(formula = logrr ~ rcs(dose, knots), type = type, 
                           id = id, se = selogrr, cases = cases, n = n, data = mine_low)
summary(spl_bin_mine_low)

# Set a refence cum RCS; in this case 4 mg/mg3-years i.e. 0.1 for 40 years
xref <- 4

# Predict risk between 0 and 20 range
pred_mine_low <- data.frame(dose = c(xref, seq(0, 5, 0.001))) %>%
  predict(spl_bin_mine_low, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2)
head(pred_mine_low)

# Plot that risk
p7 <- ggplot(pred_mine_low, aes(dose, pred, ymin = ci.lb, ymax = ci.ub)) +
  geom_ribbon(alpha = .05) +
  geom_smooth(color = "black", linewidth = 0.5) + 
  scale_y_continuous() +
  labs(x = expression(RCS ~ mg/m^3 ~ "-years"), y = "Relative risk of silicosis (log-scale)") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  ggtitle("Miners")
p7

# Table by increases of 1 mg/m3
xref <- 4
pred_mine_tab_low <- data.frame(dose = c(xref, c(1,2,4))) %>%
  predict(spl_bin_mine_low, newdata = ., expo = T, level = 0.95) %>% 
  mutate(dose = `rcs(dose, knots)dose`) %>% 
  dplyr::select(-`rcs(dose, knots)dose`,-`rcs(dose, knots)dose'`, -`rcs(dose, knots)dose''`) %>% 
  unique() %>% round(2) %>% arrange(dose)
pred_mine_tab_low


###################
# Impact analysis #
###################

# Remove ush and sag from wide_result_df
mine_prop <- wide_result_df %>% 
  filter(!id %in% c("ush", "sag_2")) %>% 
  subset(pop == "mine")

# Create non-mine data frame
non_mine_prop <- subset(wide_result_df, pop != "mine")

# Summary statistics for mine_prop and non_mine_prop
summary(mine_prop)
summary(non_mine_prop)

# Display pred_non_mine_tab and pred_mine_tab data
pred_non_mine_tab
pred_mine_tab

# Calculate medians for mine_prop and non_mine_prop
med_mine <- median(mine_prop$`4` / 100)
med_mine
med_non_mine <- median(non_mine_prop$`4` / 100)
med_non_mine

# Calculate absolute reductions for mine_prop and non_mine_prop
abs_red_mine <- data.frame(rr = c(0.18, 0.23, 0.29))
abs_red_mine$arr <- 1000 * med_mine * (1 - abs_red_mine$rr)
round(abs_red_mine, 2)

abs_red_nonmine <- data.frame(rr = c(0.36, 0.55, 0.83))
abs_red_nonmine$arr <- 1000 * med_non_mine * (1 - abs_red_nonmine$rr)
round(abs_red_nonmine, 2)
