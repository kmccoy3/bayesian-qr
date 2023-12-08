
###############################################################################

library(quantreg)
library(Brq)
library(bayesQR)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(gridExtra)

###############################################################################

# Data Sources
# https://stats.oecd.org/Index.aspx?DataSetCode=BLI
# https://www.ncdrisc.org/data-downloads-adiposity.html

# Load data
bmi_data <- read_csv("NCD_RisC_Lancet_2017_BMI_female_age_specific_country.csv")
qol_data <- read_csv("BLI_10112023062552408.csv")

# Clean data
bmi_data %>% 
  filter(Year==2016) %>%
  rename(BMI = "Mean body mass index") %>%
  select(Country, BMI, "Age group") -> bmi_data
qol_data %>%
  filter(Inequality %in% c("Total", "Women")) %>%
  select(Country, Indicator, Value) %>% 
  pivot_wider(names_from = Indicator, values_from=Value, values_fn=last ) -> qol_data



df_list <- list(bmi_data, qol_data)

#merge all data frames in list
df_list %>% reduce(inner_join, by='Country') %>% 
  mutate(across(where(is.numeric), scale)) %>%
  select(-Country) -> df

set.seed(2023)

nn <- 1500

model = bayesQR(BMI ~ ., data=df, quantile=c(0.05, 0.25, 0.5, 0.75, 0.95), alasso=TRUE, ndraw=nn)


summary(model)


###############################################################################



library("RColorBrewer")

quants <- c("0.05", "0.25", "0.5", "0.75", "0.95")
quants_num <- c(0.05, 0.25, 0.5, 0.75, 0.95)

df_hists2 <- data.frame()
for (j in 1:5){
  df_temp <- data.frame(name=quants[j], val=model[[j]][["betadraw"]][, 1])
  df_hists2 <- rbind(df_hists2, df_temp)
  
}
plot1 <- ggplot(df_hists2, aes(x=val, group=name, fill=name, after_stat(density))) +
  geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + 
  ggtitle(TeX('$\\beta_{0}$ Posterior Distribution')) + xlab(TeX("$\\beta_0$ (Intercept)")) +
  ylab("density") + labs(fill=TeX('Quantile')) + 
  scale_fill_brewer(palette = 'RdYlGn', direction=-1)

###############################################################################

df_hists2 <- data.frame()
for (j in 1:5){
  df_temp <- data.frame(name=quants[j], val=model[[j]][["betadraw"]][, 20])
  df_hists2 <- rbind(df_hists2, df_temp)
  
}
plot2 <- ggplot(df_hists2, aes(x=val, group=name, fill=name, after_stat(density))) +
  geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + 
  ggtitle(TeX('$\\beta_{20}$ Posterior Distribution')) + xlab(TeX("$\\beta_{20}$ (Feeling Safe Walking at Night)")) +
  ylab("density") + labs(fill=TeX('Quantile')) + 
  scale_fill_brewer(palette = 'RdYlGn', direction=-1)

###############################################################################

df_hists2 <- data.frame()
for (j in 1:5){
  df_temp <- data.frame(name=quants[j], val=model[[j]][["betadraw"]][, 22])
  df_hists2 <- rbind(df_hists2, df_temp)
  
}
plot3 <- ggplot(df_hists2, aes(x=val, group=name, fill=name, after_stat(density))) +
  geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + 
  ggtitle(TeX('$\\beta_{22}$ Posterior Distribution')) + xlab(TeX("$\\beta_{22}$ (Disposable Income)")) +
  ylab("density") + labs(fill=TeX('Quantile')) + 
  scale_fill_brewer(palette = 'RdYlGn', direction=-1)

###############################################################################

df_hists2 <- data.frame()
for (j in 1:5){
  df_temp <- data.frame(name=quants[j], val=model[[j]][["betadraw"]][, 34])
  df_hists2 <- rbind(df_hists2, df_temp)
  
}
plot4 <- ggplot(df_hists2, aes(x=val, group=name, fill=name, after_stat(density))) +
  geom_histogram(color='#e9ecef', alpha=0.6, position='identity') + 
  ggtitle(TeX('$\\beta_{34}$ Posterior Distribution')) + xlab(TeX("$\\beta_{34}$ (Life Expectancy)")) +
  ylab("density") + labs(fill=TeX('Quantile')) + 
  scale_fill_brewer(palette = 'RdYlGn', direction=-1)

###############################################################################

whole_plot <- grid.arrange(plot1, plot2, plot3, plot4, ncol=2)

ggsave('real_data.pdf', whole_plot)
  








