# TODO:
# Eventually try Bayesian to account for alive adults over time and eggs on leaf bottoms
# For now, split egg column into top and bottom, otherwise egg rates jump incorrectly if bottoms counted
# Mean Dead over count period or just lag-1 dead so that # adults doesn't equal 0 for division of oviprate?

library(tidyverse)
library(lubridate)
library(viridis)

dat <- readxl::read_xlsx("C:/Users/Tyson/Downloads/Aphalara_GC_2018.xlsx")

dat <- dat %>% 
  group_by(Population, Temperature, Photoperiod, EggAccumStartDate) %>% 
  arrange(Date) %>% 
  mutate(NewEggs = c(NA, diff(Eggs)),
         Dead = as.numeric(Dead)) %>%
  mutate(NewEggs = ifelse(NewEggs < 0, 0, NewEggs)) %>% 
  mutate(OvipRate = NewEggs / 
           as.numeric(difftime(Date, lag(Date, 1), units = "days")) / 
           (StartCount - ifelse(is.na(Dead), 0, Dead)),
         OvipRateDD = NewEggs / 
           c(NA, diff(AccumDD)) / 
           (StartCount - ifelse(is.na(Dead), 0, Dead))) %>% 
  filter(is.na(NewEggs) == FALSE)


plt <- ggplot(dat, aes(x = AccumDD, y = OvipRate, group = Photoperiod, color = Photoperiod)) +
  geom_point(size = 2.5, alpha = .5) +
  scale_color_viridis() +
  facet_wrap(Population~Temperature) +
  theme_bw() +
  ggtitle("Oviposition rate (new eggs / # days / # adults)")
plt


plt <- ggplot(dat, aes(x = Photoperiod, y = OvipRate, color = AccumDD)) +
  geom_point(size = 2.5, alpha = .5) +
  scale_color_viridis() +
  facet_wrap(Population~Temperature) +
  theme_bw() +
  ggtitle("Oviposition rate (new eggs / # days / # adults)")
plt
