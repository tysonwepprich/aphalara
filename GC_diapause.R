# GC diapause

library(tidyverse)
library(lubridate)
library(viridis)

dat <- readxl::read_xlsx("C:/Users/wepprict/Downloads/Aphalara_GC_2018_diapause.xlsx",
                         na = "NA")

dat <- dat %>% 
  rowwise() %>% 
  mutate(Eggs = sum(c(Eggs_2, Eggs_3, Eggs_4), na.rm = TRUE)) %>%
  group_by(Population, Temperature, Photoperiod) %>% 
  mutate(NumMales = sum(Males, na.rm = TRUE),
         NumFemales = sum(Females, na.rm = TRUE)) %>% 
  filter(Females >= 1) %>% 
  mutate(Decision = ifelse(Eggs > 0, "Reproductive", "Diapause"))
  

theme_set(theme_bw(base_size = 20)) 

plt <- ggplot(dat, aes(x = Photoperiod, fill = Decision)) +
  geom_bar(width = 0.9) +
  # scale_color_viridis() +
  facet_wrap(Population~Temperature) +
  xlab("Photoperiod") +
  ylab("Count of diapause decision") +
  ggtitle("Proportion of females reproductive at end of test")
plt
  