library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggpmisc)
library(cowplot)
library(ggplot2)

colors <- list(
  tone = c("#f6f2ed","#dedbca","#c4c0a5","#a59f83","#878264","#5e5948"),
  grey = c("#e5e4e8","#c5c9d6","#959fb3","#6e788c","#425468","#1b2942"),
  olive = c("#f2edb2","#dbdc63","#c4c400","#95a008","#637314","#304215"),
  green = c("#d7e5c5","#9fc978","#5db342","#41912f","#1d6e29","#0e3716"),
  teal = c("#c9e4ef","#96ced3","#48bcbc","#00959f","#006479","#003547"),
  blue = c("#c5e4fb","#9bc9e8","#5495ce","#006eae","#01478c","#002259"),
  purple = c("#e9d3e7","#d1a9ce","#b678b3","#a44990","#792373","#430b4d"),
  red = c("#f5cfc9","#e9a0a4","#db6463","#c5373c","#9c241b","#730c0d"),
  orange = c("#fbdcbc","#f9bd7b","#f29741","#e96700","#b34a00","#832a00"),
  yellow = c("#ffedc1","#f7dc86","#e8c54d","#ca9a23","#9b730a","#685409"),
  skin = c("#f6e5d3","#dcbc9f","#bc9778","#906852","#734d3d","#422a17"))

#comparison of biological age
data <- readxl::read_xlsx("./data.xlsx")
data %>% filter(Timepoint=="L-45") %>% dplyr::select(Crew.ID,Female,Age)
age <- data %>% dplyr::select(Patient.ID,Crew.ID,Age,Female,Timepoint) %>% left_join(reshape2::melt(data[,c(2,which(colnames(data)%in%clocks))]) %>% set_names("Patient.ID","Clock","BA") %>% mutate(Patient.ID=as.character(Patient.ID), Clock = as.character(Clock)))
colnames(age) <- c("patient","crew","age","sex","timepoint","clock","ba")
age$timepoint <- factor(age$timepoint, levels = c("L-45","FD+4","FD+7","R+1","R+7"))
age <- age %>% group_by(crew, age, sex, timepoint) %>% summarise(n = n(), mean_ba = mean(ba), sd_ba = sd(ba), se_ba = sd_ba / sqrt(n))

day4 <- age %>% filter(timepoint%in%c("L-45")) %>% ungroup %>% dplyr::select(crew,mean_ba) %>% set_names("crew","mean_pre") %>% 
  left_join(age %>% filter(timepoint%in%c("FD+4")) %>% ungroup %>% dplyr::select(crew,mean_ba) %>% set_names("crew","mean_post")) %>% 
  mutate(dif = round(mean_post-mean_pre,2)) %>% arrange(desc(dif))
day4
mean(day4$dif)

day7 <- age %>% filter(timepoint%in%c("L-45")) %>% ungroup %>% dplyr::select(crew,mean_ba) %>% set_names("crew","mean_pre") %>% 
  left_join(age %>% filter(timepoint%in%c("FD+7")) %>% ungroup %>% dplyr::select(crew,mean_ba) %>% set_names("crew","mean_post")) %>% 
  mutate(dif = round(mean_post-mean_pre,2)) %>% arrange(desc(dif))
day7
mean(day7$dif)

ret1 <- age %>% filter(timepoint%in%c("FD+7")) %>% ungroup %>% dplyr::select(crew,mean_ba) %>% set_names("crew","mean_pre") %>% 
  left_join(age %>% filter(timepoint%in%c("R+1")) %>% ungroup %>% dplyr::select(crew,mean_ba) %>% set_names("crew","mean_post")) %>% 
  mutate(dif = round(mean_post-mean_pre,2)) %>% arrange(desc(dif))
ret1
mean(ret1$dif)

permutation_test_dependent <- function(tp1, tp2, n_perm = 1000000) {
  diffs <- tp2 - tp1
  obs_stat <- mean(diffs)
  
  perm_stats <- replicate(n_perm, {
    signs <- sample(c(1, -1), length(diffs), replace = TRUE)
    mean(diffs * signs)
  })
  mean(abs(perm_stats) >= abs(obs_stat))
}

comparison <- data %>% dplyr::select(Patient.ID,Crew.ID,Age,Female,Timepoint) %>% left_join(reshape2::melt(data[,c(2,which(colnames(data)%in%clocks))]) %>% set_names("Patient.ID","Clock","BA") %>% mutate(Patient.ID=as.character(Patient.ID), Clock = as.character(Clock)))
colnames(comparison) <- c("patient","crew","age","sex","timepoint","clock","ba")
comparison <- comparison %>% filter(timepoint=="L-45") %>% dplyr::select(crew,clock,ba) %>% set_names("crew","clock","pre") %>% 
  left_join(comparison %>% filter(timepoint!="L-45") %>% dplyr::select(crew,clock,timepoint,ba) %>% set_names("crew","clock","timepoint","post"))
set.seed(1234)
comparison <- comparison %>% group_by(crew,timepoint) %>% summarise(p = permutation_test_dependent(pre,post))
comparison$timepoint <- factor(comparison$timepoint, levels = c("FD+4","FD+7","R+1","R+7"))

p2star <- function(x){
  symnum(
    x,
    corr = FALSE,
    cutpoints = c(0,   0.001, 0.01, 0.05, 1),
    symbols   = c("***","**",  "*",  "ns")
  )
}

pdf(file = "./figure1.pdf", width=6, height=5) 
ggplot(age, aes(timepoint, mean_ba, group = crew, fill = crew))+
  geom_hline(aes(yintercept = age, color = crew), lty = 2)+
  geom_rect(
    data = group_colors_ax,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = c('grey80','white','grey80'),
    inherit.aes = FALSE, alpha = 0.5
  ) +
  geom_line(color = "black")+
  geom_errorbar(aes(ymin = mean_ba - se_ba, ymax = mean_ba + se_ba), width = 0.1) +
  geom_point(shape = 21, size = 3)+
  
  geom_text(data = age %>% filter(crew%in%c("A1")), aes(label = c("",p2star(comparison %>% filter(crew=="A1") %>% pull(p)))), vjust = -1.5, size = 5, color = "black") +
  geom_text(data = age %>% filter(crew%in%c("A2")), aes(label = c("",p2star(comparison %>% filter(crew=="A2") %>% pull(p)))), vjust = 3, size = 5, color = "black") +
  geom_text(data = age %>% filter(crew%in%c("A3")), aes(label = c("",p2star(comparison %>% filter(crew=="A3") %>% pull(p)))), vjust = 3, size = 5, color = "black") +
  geom_text(data = age %>% filter(crew%in%c("A4")), aes(label = c("",p2star(comparison %>% filter(crew=="A4") %>% pull(p)))), vjust = -1.5, size = 5, color = "black") +
  
  scale_fill_manual(values = c(colors$green[2],
                               colors$orange[2],
                               colors$blue[2],
                               colors$purple[2]))+
  scale_color_manual(values = c(colors$green[2],
                                colors$orange[2],
                                colors$blue[2],
                                colors$purple[2]))+
  guides(color = "none")+
  theme_pubr(border = TRUE)+
  theme(legend.position = "right")+
  labs(x = "Time point", y = "Epigenetic age", fill = "Crew")
dev.off()

#linear mixed model
age <- data %>% dplyr::select(Patient.ID,Crew.ID,Age,Female,Timepoint) %>% left_join(reshape2::melt(data[,c(2,which(colnames(data)%in%clocks))]) %>% set_names("Patient.ID","Clock","BA") %>% mutate(Patient.ID=as.character(Patient.ID), Clock = as.character(Clock)))
colnames(age) <- c("patient","crew","age","sex","timepoint","clock","ba")
age <- age %>% filter(timepoint=="L-45") %>% dplyr::select(crew,age,sex,clock,ba) %>% set_names("crew","age","sex","clock","pre") %>% 
  left_join(age %>% filter(timepoint=="FD+7") %>% dplyr::select(crew,clock,ba) %>% set_names("crew","clock","during"))
age$ba_change <- age$during-age$pre
model <- lmer(ba_change ~ age + (1|clock) + (1|crew), data = age)
summary(model)

