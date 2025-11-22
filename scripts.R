library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggpmisc)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(gt)
library(lme4)
library(lmerTest)
library(dominanceanalysis)

get_group <- function(clock) {
  switch(clock,
         "Horvath" = "Chronological age",
         "Hannum" = "Chronological age",
         "PCHorvath1" = "Chronological age",
         "PCHorvath2" = "Chronological age",
         "PCHannum" = "Chronological age",
         
         "PhenoAge" = "Mortality",
         "OMICmAge" = "Mortality",
         "PCPhenoAge" = "Mortality",
         "PCGrimAge" = "Mortality",
         
         "DNAmFitAge" = "Physical fitness",
         
         "AdaptAge" = "Causal factors",
         "CausAge" = "Causal factors",
         "DamAge" = "Causal factors",
         
         "IntrinClock" = "Intrinsic age",
         
         "Stochastic.Zhang" = "Stochasticity",
         "Stochastic.Horvath" = "Stochasticity",
         "Stochastic.PhenoAge" = "Stochasticity",
         
         "Retroclock" = "Retroelements",
         "Retroclockv2" = "Retroelements",
         
         "Blood" = "Organ aging",
         "Brain" = "Organ aging",
         "Inflammation" = "Organ aging",
         "Heart" = "Organ aging",
         "Hormone" = "Organ aging",
         "Immune" = "Organ aging",
         "Kidney" = "Organ aging",
         "Liver" = "Organ aging",
         "Metabolic" = "Organ aging",
         "Lung" = "Organ aging",
         "MusculoSkeletal" = "Organ aging",
         "SystemsAge" = "Organ aging",
         
         "Unknown"  # default
  )
}

permutation_test_dependent <- function(tp1, tp2, n_perm = 1e4) {
  diffs <- tp2 - tp1
  obs_stat <- mean(diffs)
  set.seed(1234)
  perm_stats <- replicate(n_perm, {
    signs <- sample(c(1, -1), length(diffs), replace = TRUE)
    mean(diffs * signs)
  })
  mean(abs(perm_stats) >= abs(obs_stat)) #two-sided
}

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


data <- read_csv("./data/final_data/Project_EpigeneticClocks_EPICv2custom.csv")
data <- data[grepl("^A",data$Crew.ID),]
data$Crew.ID <- gsub("00","",data$Crew.ID)
data$Timepoint <- gsub("^Flight Day ","FD+",data$Timepoint)
sel_timepoints <- colSums(table(data$Crew.ID,data$Timepoint))
sel_timepoints <- names(sel_timepoints)[sel_timepoints==4]
data <- data %>% filter(Timepoint%in%sel_timepoints)
data %>% filter(Timepoint=="L-45") %>% dplyr::select(Crew.ID,Female,Age)
colnames(data)[1] <- "Patient.ID"
colnames(data)[12] <- "Decimal.Chronological.Age"
colnames(data)[14] <- "Age"

################## Age acceleration calculations ##################
clocks <- c("Age",
            "Horvath", "Hannum", "PCHorvath1", "PCHorvath2", "PCHannum", #chronological age
            "PhenoAge","OMICmAge","PCPhenoAge", "PCGrimAge", #mortality
            "DNAmFitAge", #physical fitness
            "AdaptAge", "DamAge", "CausAge", #causal factors
            "IntrinClock", #intrinsic age
            "Stochastic.Zhang", "Stochastic.PhenoAge", "Stochastic.Horvath", #stochasticity
            "Retroclock", "Retroclockv2", #retroelements
            "Blood", "Brain", "Inflammation", "Heart", "Hormone", "Immune", "Kidney", "Liver", "Metabolic", "Lung", "MusculoSkeletal", "SystemsAge") #Organ aging

#method 1
all_estimates <- tibble()
clocks <- clocks[2:length(clocks)]
temp <- reshape2::melt(data[,c(1,which(colnames(data)%in%clocks))]) %>% set_names("Patient.ID","Clock","BA") %>% mutate(Patient.ID=as.character(Patient.ID), Clock = as.character(Clock), BA = as.numeric(BA))
age <- data %>% dplyr::select(Patient.ID,Crew.ID,Age,Female,Timepoint) %>% left_join(temp) 
colnames(age) <- c("patient","crew","age","sex","timepoint","clock","ba")
age$ba <- age$ba-age$age
age$system <- "Age difference"
all_estimates <- rbind(all_estimates,age)

#method 2
temp_data <- data
for (clock in clocks) {
  f <- as.formula(paste0(clock, " ~ Decimal.Chronological.Age + Female"))
  model <- lm(f, data = temp_data)
  temp_data[[clock]] <- residuals(model)
}
temp <- reshape2::melt(temp_data[,c(1,which(colnames(temp_data)%in%clocks))]) %>% set_names("Patient.ID","Clock","BA") %>% mutate(Patient.ID=as.character(Patient.ID), Clock = as.character(Clock), BA = as.numeric(BA))
age <- temp_data %>% dplyr::select(Patient.ID,Crew.ID,Age,Female,Timepoint) %>% left_join(temp) 
colnames(age) <- c("patient","crew","age","sex","timepoint","clock","ba")
age$system <- "Age acceleration"
all_estimates <- rbind(all_estimates,age)

#method 3
cells <- c("CD4Tnv","CD4Tmem","CD8Tnv","CD8Tmem","Bnv","Bmem","NK","Mono","Eos","Baso","Neu","Treg")
temp_data <- data
for (clock in clocks) {
  f <- as.formula(paste0(clock, " ~ Decimal.Chronological.Age  + Female + ",paste(cells, collapse = " + ")))
  model <- lm(f, data = temp_data)
  temp_data[[clock]] <- residuals(model)
}
temp <- reshape2::melt(temp_data[,c(1,which(colnames(temp_data)%in%clocks))]) %>% set_names("Patient.ID","Clock","BA") %>% mutate(Patient.ID=as.character(Patient.ID), Clock = as.character(Clock), BA = as.numeric(BA))
age <- temp_data %>% dplyr::select(Patient.ID,Crew.ID,Age,Female,Timepoint) %>% left_join(temp) 
colnames(age) <- c("patient","crew","age","sex","timepoint","clock","ba")
age$system <- "Intrinsic age acceleration"
all_estimates <- rbind(all_estimates,age)

#adjust order of labels
all_estimates$timepoint <- factor(all_estimates$timepoint, levels = c("L-45","FD+4","FD+7","R+1","R+7"))
all_estimates$system <- factor(all_estimates$system, levels = c("Age difference","Age acceleration","Intrinsic age acceleration"))
all_estimates_temp <- all_estimates
all_estimates <- all_estimates %>% group_by(crew, age, sex, timepoint, system) %>% summarise(n = n(), mean_ba = mean(ba), sd_ba = sd(ba), se_ba = sd_ba / sqrt(n))

################## Figure Supplementary 1 ##################

correlations_age <- cor(data[,c("Age",clocks)])[,"Age"] %>% enframe %>% arrange(desc(value))
selected <- reshape2::melt(data[,c("Crew.ID",clocks)]) %>% left_join(correlations_age %>% set_names("variable","r"))
selected$variable <- factor(selected$variable, levels = correlations_age$name %>% rev)
selected <- selected %>% left_join(data %>% group_by(Crew.ID) %>% summarise(real_age = mean(Decimal.Chronological.Age)) %>% left_join(selected %>% group_by(Crew.ID) %>% summarise(epi_age = mean(value))))
selected$Crew.ID <- factor(selected$Crew.ID, levels = c("A1","A2","A3","A4"))

pdf("./output/Supplementary_Figure1.pdf", width=7.54, height=7.99)
ggplot(selected, aes(variable,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha = 0.1)+
  #geom_hline(aes(yintercept = real_age, color = Crew.ID), lwd = 1)+
  #geom_hline(aes(yintercept = epi_age, color = Crew.ID), lwd = 1, lty = 2)+
  geom_text(y = 90, aes(label = round(r,2)), stat = "unique", color = "red")+
  scale_color_manual(values = c(colors$green[3],
                               colors$orange[3],
                               colors$blue[3],
                               colors$purple[3]))+
  geom_hline(yintercept = 83, lwd = 0.25)+
  annotate("text", 
             label = "Age correlation", 
             x = 15,     
             y = 85,    
             angle = 90,  
             vjust = 0.5, 
             color = "black", 
             size = 4)+      
  coord_flip()+
  ylim(0,90)+
  theme_pubr(border = TRUE)+
  labs(x = "Epigenetic clock", y = "Age (years)", color = "Crew")
dev.off()

################## Figure Supplementary 2 ##################

average_timepoints <- all_estimates_temp %>% mutate(patient=interaction(crew,timepoint)) %>% group_by(system,clock,patient) %>%  summarise(ba = mean(ba))
avg_method1 <- reshape2::acast(average_timepoints %>% filter(system=="Age difference"), clock~patient, value.var = "ba")
avg_method2 <- reshape2::acast(average_timepoints %>% filter(system=="Age acceleration"), clock~patient, value.var = "ba")
avg_method3 <- reshape2::acast(average_timepoints %>% filter(system=="Intrinsic age acceleration"), clock~patient, value.var = "ba")

cor_method1_individuals <- cor(avg_method1)
cor_method2_individuals <- cor(avg_method2)
cor_method3_individuals <- cor(avg_method3)

col_fun <- colorRamp2(seq(-1,1,0.2), c(colors$blue[5:1],"white",colors$red[1:5]))

ht1_individuals <- Heatmap(
  cor_method1_individuals,
  rect_gp = gpar(col = "black", lwd = 1, lty = 3),
  border = TRUE,
  name = "R",
  column_title = "Age difference",
  row_dend_side = "right",
  column_dend_side = "top",
  row_names_side = "left",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_heatmap_legend = FALSE,
  col = col_fun
)

ht2_individuals <- Heatmap(
  cor_method2_individuals,
  rect_gp = gpar(col = "black", lwd = 1, lty = 3),
  border = TRUE,
  name = "R",
  column_title = "Age acceleration",
  row_dend_side = "right",
  column_dend_side = "top",
  row_names_side = "left",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_heatmap_legend = FALSE,
  col = col_fun
)

ht3_individuals <- Heatmap(
  cor_method3_individuals,
  rect_gp = gpar(col = "black", lwd = 1, lty = 3),
  border = TRUE,
  name = "R",
  column_title = "Intrinsic age acceleration",
  row_dend_side = "right",
  column_dend_side = "top",
  row_names_side = "left",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_heatmap_legend = TRUE,
  col = col_fun
)

p1_individuals <- grid.grabExpr({draw(ht1_individuals)})
p2_individuals <- grid.grabExpr({draw(ht2_individuals)})
p3_individuals <- grid.grabExpr({draw(ht3_individuals)})

inter_intra_method1 <- reshape2::melt(cor_method1_individuals) %>% 
  filter(value!=1) %>%
  separate(col = "Var1", into = c("sample1","timepoint1"), sep = "\\.") %>%
  separate(col = "Var2", into = c("sample2","timepoint2"), sep = "\\.")
inter_intra_method1 <- rbind(inter_intra_method1 %>% filter(sample1==sample2) %>% mutate(individual = "Same\nindividual", method = "Age difference"),
                             inter_intra_method1 %>% filter(sample1!=sample2) %>% mutate(individual = "Different\nindividual", method = "Age difference"),
                             inter_intra_method1 %>% filter(timepoint1==timepoint2) %>% mutate(individual = "Same\ntimepoint", method = "Age difference"),
                             inter_intra_method1 %>% filter(timepoint1!=timepoint2) %>% mutate(individual = "Different\ntimepoint", method = "Age difference"))

inter_intra_method2 <- reshape2::melt(cor_method2_individuals) %>% 
  filter(value!=1) %>%
  separate(col = "Var1", into = c("sample1","timepoint1"), sep = "\\.") %>%
  separate(col = "Var2", into = c("sample2","timepoint2"), sep = "\\.")
inter_intra_method2 <- rbind(inter_intra_method2 %>% filter(sample1==sample2) %>% mutate(individual = "Same\nindividual", method = "Age acceleration"),
                             inter_intra_method2 %>% filter(sample1!=sample2) %>% mutate(individual = "Different\nindividual", method = "Age acceleration"),
                             inter_intra_method2 %>% filter(timepoint1==timepoint2) %>% mutate(individual = "Same\ntimepoint", method = "Age acceleration"),
                             inter_intra_method2 %>% filter(timepoint1!=timepoint2) %>% mutate(individual = "Different\ntimepoint", method = "Age acceleration"))

inter_intra_method3 <- reshape2::melt(cor_method3_individuals) %>% 
  filter(value!=1) %>%
  separate(col = "Var1", into = c("sample1","timepoint1"), sep = "\\.") %>%
  separate(col = "Var2", into = c("sample2","timepoint2"), sep = "\\.")
inter_intra_method3 <- rbind(inter_intra_method3 %>% filter(sample1==sample2) %>% mutate(individual = "Same\nindividual", method = "Intrinsic age acceleration"),
                             inter_intra_method3 %>% filter(sample1!=sample2) %>% mutate(individual = "Different\nindividual", method = "Intrinsic age acceleration"),
                             inter_intra_method3 %>% filter(timepoint1==timepoint2) %>% mutate(individual = "Same\ntimepoint", method = "Intrinsic age acceleration"),
                             inter_intra_method3 %>% filter(timepoint1!=timepoint2) %>% mutate(individual = "Different\ntimepoint", method = "Intrinsic age acceleration"))
inter_intra <- rbind(inter_intra_method1,inter_intra_method2,inter_intra_method3)
inter_intra$method <- factor(inter_intra$method, levels = c("Age difference", "Age acceleration", "Intrinsic age acceleration"))

p123_comparison <- ggplot(inter_intra, aes(individual, value))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, color = "grey", alpha = 0.5)+
  facet_wrap(.~method, nrow = 1)+
  geom_hline(yintercept = 0, lty = 2)+
  theme_pubr(border = TRUE)+
  labs(x = "", y = "Average correlation coefficient")

library(cowplot)
pdf(file = "./output/Supplementary_Figure2.pdf", width=18.40, height=12.69)
plot_grid(plot_grid(p1_individuals,p2_individuals,p3_individuals, nrow = 1, rel_widths = c(1,1,1.1)),
          p123_comparison, nrow = 2)
dev.off()

################## Figure 1a ##################

group_colors_ax <- data.frame(
  xmin = c(0.5, 1.5, 3.5),
  xmax = c(1.5, 3.5, 5.5),
  group = c('L', 'FD', 'R')
)

p1 <- ggplot(all_estimates, aes(timepoint, mean_ba, group = crew, fill = crew))+
  geom_hline(yintercept = 0, color = "black", lty = 2)+
  geom_rect(
    data = group_colors_ax,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = c('grey80','white','grey80','grey80','white','grey80','grey80','white','grey80'),
    inherit.aes = FALSE, alpha = 0.5
  ) +
  geom_line(color = "black")+
  facet_wrap(system~., ncol = 1)+
  geom_point(shape = 21, size = 3)+
  scale_fill_manual(values = c(colors$green[2],
                                colors$orange[2],
                                colors$blue[2],
                                colors$purple[2]))+
  guides(color = "none")+
  theme_pubr(border = TRUE)+
  labs(x = "Time point", y = "Average value", fill = "Crew")

################## Figure 1b ##################

comparison <- all_estimates_temp %>% filter(timepoint=="L-45") %>% ungroup %>% dplyr::select(system,crew,clock,ba) %>% set_names("system","crew","clock","pre") %>% group_by(system,clock) %>% summarise(pre = mean(pre)) %>% 
  left_join(all_estimates_temp %>% filter(timepoint!="L-45") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","timepoint","clock","post") %>% group_by(system,clock,timepoint) %>% summarise(post = mean(post)))

all_clocks <- comparison %>% mutate(dif = post - pre) 
all_clocks <- all_clocks %>% dplyr::select(system,timepoint,clock,dif) %>% mutate(column = interaction(timepoint,system))
mat_clocks <- reshape2::acast(all_clocks, clock~column, value.var = "dif")

col_fun <- colorRamp2(seq(-7,7,14/10), c(colors$blue[5:1],"white",colors$red[1:5]))

cs_idx <- rep(1:3, each = 4)
sub  <- rep(rep(c("EAD","EAA\nA1","IEAA"), each = 4), times = 1)
col_split <- data.frame(Sub = factor(sub,  levels = c("EAD","EAA\nA1","IEAA")))
gaps <- rep(unit(1, "mm"), 12)

ha = HeatmapAnnotation(
  Timepoint = anno_simple(
    rep(c("FD+4","FD+7","R+1","R+7"), 3),
    col   = cols_time <- c("FD+4" = colors$green[2],"FD+7" = colors$green[3],"R+1"  = colors$orange[2],"R+7"  = colors$orange[3]),        # fill color by timepoint
    pch   = rep(c("FD+4","FD+7","R+1","R+7"), 3),       # text drawn inside each square
    pt_size = unit(1, "snpc")*0.55,
    border = FALSE
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 11, fontface = "bold")
)
#mat_clocks <- mat_clocks[rowMeans(mat_clocks[,c(1:2,5:6,9:10),drop = FALSE]) %>% enframe %>% arrange(desc(value)) %>% pull(name),]

hb = rowAnnotation(Category = sapply(rownames(mat_clocks), get_group),
                       col = list(Category = setNames(c(colors$tone[3],colors$grey[3],colors$olive[3],colors$teal[3],colors$purple[3],colors$yellow[3],colors$skin[3],colors$blue[3]), sapply(rownames(mat_clocks), get_group) %>% unique)))


ht <- Heatmap(
  mat_clocks,
  rect_gp = gpar(col = "black", lwd = 1, lty = 3),
  border = TRUE,
  name = "Age\ndifference",
  top_annotation = ha,
  left_annotation = hb,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(label = round(mat_clocks[i,j],2), x = x, y = y, gp = gpar(fontsize = 7.5))
  },
  column_split = cs_idx,
  cluster_columns = FALSE,
  column_title = c("Age\ndifference","Age\nacceleration","Intrinsic age\nacceleration"),
  row_dend_side = "right",
  row_names_side = "left",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_heatmap_legend = TRUE,
  col = col_fun
)
draw(ht)
p2 <- grid.grab()

################## Plot 1a & 1b ##################

library(cowplot)
pdf(file = "./output/Figure1.pdf", width=11, height=10)
plot_grid(p1,p2, rel_widths = c(1,2), labels = c("a","b"))
dev.off()

################## Supplementary Table 3 ##################

comparison1 <- all_estimates_temp %>% filter(timepoint=="L-45") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","t1","clock","pre") %>% 
  left_join(all_estimates_temp %>% filter(timepoint!="L-45") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","t2","clock","post"))

comparison2 <- all_estimates_temp %>% filter(timepoint=="FD+4") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","t1","clock","pre") %>% 
  left_join(all_estimates_temp %>% filter(timepoint=="FD+7") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","t2","clock","post"))

comparison3 <- all_estimates_temp %>% filter(timepoint=="FD+7") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","t1","clock","pre") %>% 
  left_join(all_estimates_temp %>% filter(timepoint=="R+1") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","t2","clock","post"))

comparison4 <- all_estimates_temp %>% filter(timepoint=="R+1") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","t1","clock","pre") %>% 
  left_join(all_estimates_temp %>% filter(timepoint=="R+7") %>% ungroup %>% dplyr::select(system,crew,timepoint,clock,ba) %>% set_names("system","crew","t2","clock","post"))

comparison <- rbind(comparison1,comparison2,comparison3,comparison4)

results <- comparison %>%
  group_by(system,crew,t1,t2) %>%
  mutate(diff = post - pre) %>%
  summarise(
    n = n(),
    mean = mean(diff),
    sd = sd(diff),
    se = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se,
    wilcox = formatC(wilcox.test(post, pre, paired = TRUE)$p.value, format = "e", digits = 2),
    perm = formatC(permutation_test_dependent(pre,post), format = "e", digits = 2),
  )

results <- results %>%
  mutate(mean = round(mean,2),
         se = round(se,2),
         ci = paste0("[",sprintf("%.2f", round(ci_lower,2)),", ",sprintf("%.2f", round(ci_upper,2)),"]")) %>%
  dplyr::select(system,crew,t1,t2,mean,se,ci,wilcox,perm)

m1 <- results %>% ungroup %>% filter(system=="Age difference") %>% dplyr::select(crew,t1,t2,mean,se,ci,wilcox,perm) %>% set_names("crew","timepoint_pre","timepoint_post","mean_1","se_1","ci_1","wilcox_1","perm_1")
m2 <- results %>% ungroup %>% filter(system=="Age acceleration") %>% dplyr::select(crew,t1,t2,mean,se,ci,wilcox,perm) %>% set_names("crew","timepoint_pre","timepoint_post","mean_2","se_2","ci_2","wilcox_2","perm_2")
m3 <- results %>% ungroup %>% filter(system=="Intrinsic age acceleration") %>% dplyr::select(crew,t1,t2,mean,se,ci,wilcox,perm) %>% set_names("crew","timepoint_pre","timepoint_post","mean_3","se_3","ci_3","wilcox_3","perm_3")

combined <- m1 %>% left_join(m2) %>% left_join(m3)
writexl::write_xlsx(combined, path = "./output/Supplementary_Table_3.xlsx")

#L-45 vs FD+4
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Age difference"), all_estimates_temp %>% filter(timepoint=="FD+4"&system=="Age difference"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Age acceleration"), all_estimates_temp %>% filter(timepoint=="FD+4"&system=="Age acceleration"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Intrinsic age acceleration"), all_estimates_temp %>% filter(timepoint=="FD+4"&system=="Intrinsic age acceleration"))))

#FD+4 vs FD+7
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="FD+4"&system=="Age difference"), all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Age difference"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="FD+4"&system=="Age acceleration"), all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Age acceleration"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="FD+4"&system=="Intrinsic age acceleration"), all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Intrinsic age acceleration"))))

#L-45 vs FD+7
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Age difference"), all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Age difference"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Age acceleration"), all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Age acceleration"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Intrinsic age acceleration"), all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Intrinsic age acceleration"))))

#FD+7 vs R+1
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Age difference"), all_estimates_temp %>% filter(timepoint=="R+1"&system=="Age difference"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Age acceleration"), all_estimates_temp %>% filter(timepoint=="R+1"&system=="Age acceleration"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="FD+7"&system=="Intrinsic age acceleration"), all_estimates_temp %>% filter(timepoint=="R+1"&system=="Intrinsic age acceleration"))))

#L-45 vs R+1
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Age difference"), all_estimates_temp %>% filter(timepoint=="R+1"&system=="Age difference"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Age acceleration"), all_estimates_temp %>% filter(timepoint=="R+1"&system=="Age acceleration"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Intrinsic age acceleration"), all_estimates_temp %>% filter(timepoint=="R+1"&system=="Intrinsic age acceleration"))))

#L-45 vs R+7
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Age difference"), all_estimates_temp %>% filter(timepoint=="R+7"&system=="Age difference"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Age acceleration"), all_estimates_temp %>% filter(timepoint=="R+7"&system=="Age acceleration"))))
summary(lmer(ba ~ timepoint + (1|clock) + (1|crew), data = rbind(all_estimates_temp %>% filter(timepoint=="L-45"&system=="Intrinsic age acceleration"), all_estimates_temp %>% filter(timepoint=="R+7"&system=="Intrinsic age acceleration"))))

################## Figure Supplementary 3 ##################

cells <- c("CD4Tnv","CD4Tmem","CD8Tnv","CD8Tmem","Bnv","Bmem","NK","Mono","Eos","Baso","Neu","Treg")
dunedin <- data %>% dplyr::select(Patient.ID,Crew.ID,Age,Female,Timepoint) %>% left_join(data[,c(1,32,which(colnames(data)%in%cells))]) 
dunedin <- dunedin %>% dplyr::select(Crew.ID,Timepoint,DunedinPACE) %>% set_names("crew","timepoint","DunedinPACE")
dunedin <- reshape2::melt(dunedin)
dunedin$timepoint <- factor(dunedin$timepoint, levels = c("L-45","FD+4","FD+7","R+1","R+7"))

dunedin %>% filter(variable=="DunedinPACE"&timepoint=="L-45") %>% dplyr::select(crew,value) %>% set_names("crew","t1") %>% 
  left_join(dunedin %>% filter(variable=="DunedinPACE"&timepoint=="FD+4") %>% dplyr::select(crew,value) %>% set_names("crew","t2")) %>%
  mutate(t = t2-t1) %>% pull(t) %>% mean

dunedin %>% filter(variable=="DunedinPACE"&timepoint=="FD+4") %>% dplyr::select(crew,value) %>% set_names("crew","t1") %>% 
  left_join(dunedin %>% filter(variable=="DunedinPACE"&timepoint=="FD+7") %>% dplyr::select(crew,value) %>% set_names("crew","t2")) %>%
  mutate(t = t2-t1) %>% pull(t) %>% mean

dunedin %>% filter(variable=="DunedinPACE"&timepoint=="FD+4") %>% dplyr::select(crew,value) %>% set_names("crew","t1") %>% 
  left_join(dunedin %>% filter(variable=="DunedinPACE"&timepoint=="R+1") %>% dplyr::select(crew,value) %>% set_names("crew","t2")) %>%
  mutate(t = t2-t1) %>% pull(t) %>% mean

dunedin %>% filter(variable=="DunedinPACE"&timepoint=="L-45") %>% dplyr::select(crew,value) %>% set_names("crew","t1") %>% 
  left_join(dunedin %>% filter(variable=="DunedinPACE"&timepoint=="R+1") %>% dplyr::select(crew,value) %>% set_names("crew","t2")) %>%
  mutate(t = t2-t1) 

pdf(file = "./output/Supplementary_Figure3.pdf", width=4.5, height=4.58)
ggplot(dunedin %>% filter(variable=="DunedinPACE"), aes(timepoint, value, group = crew, fill = crew))+
  geom_rect(
    data = group_colors_ax,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = c('grey80','white','grey80'),
    inherit.aes = FALSE, alpha = 0.5
  ) +
  geom_line(color = "black")+
  geom_point(shape = 21, size = 3)+
  scale_fill_manual(values = c(colors$green[2],
                               colors$orange[2],
                               colors$blue[2],
                               colors$purple[2]))+
  guides(color = "none")+
  theme_pubr(border = TRUE)+
  labs(x = "Time point", y = "DunedinPACE", fill = "Crew")
dev.off()

################## Supplementary Table 4 ##################

cells <- c("CD4Tnv","CD4Tmem","CD8Tnv","CD8Tmem","Bnv","Bmem","NK","Mono","Eos","Baso","Neu","Treg")

cell_effect <- all_estimates %>% filter(system=="Age acceleration") %>% ungroup %>% dplyr::select(crew,timepoint,mean_ba) %>% set_names("crew","timepoint","mean_2") %>% 
  left_join(all_estimates %>% filter(system=="Intrinsic age acceleration") %>% ungroup %>% dplyr::select(crew,timepoint,mean_ba) %>% set_names("crew","timepoint","mean_3")) %>% 
  mutate(cell_effect = mean_3-mean_2) %>%
  dplyr::select(crew,timepoint,cell_effect)

temp_data <- data[,c("Crew.ID", "Timepoint", cells)] 
colnames(temp_data)[1] <- "crew"
colnames(temp_data)[2] <- "timepoint"
temp_data <- temp_data %>% left_join(cell_effect) %>% filter(timepoint!="R+7")
temp_data <- temp_data[,3:ncol(temp_data)]

model <- lm(cell_effect ~ ., data = temp_data)
da <- dominanceAnalysis(model)
da$contribution.average$r2 %>% enframe %>% arrange(desc(value)) %>% set_names("Cell type","R-squared")
writexl::write_xlsx(da$contribution.average$r2 %>% enframe %>% arrange(desc(value)) %>% set_names("Cell type","R-squared"), path = "./output/Supplementary_Table_4.xlsx")

plot_cell <- reshape2::melt(data[,c("Crew.ID","Timepoint",cells)]) %>% left_join(tibble(data[,c("Crew.ID","Timepoint")], mean = rowMeans(data[,clocks])) )
plot_cell$Timepoint <- factor(plot_cell$Timepoint, levels = c("L-45","FD+4","FD+7","R+1","R+7"))

pdf(file = "./output/Supplementary_Figure4.pdf", width=13.26, height=9.18)
ggplot(plot_cell, aes(Timepoint, value, color = Crew.ID, group = Crew.ID))+
  geom_point()+
  geom_line()+
  scale_color_manual(values = c(colors$green[2],
                               colors$orange[2],
                               colors$blue[2],
                               colors$purple[2]))+
  theme_pubr(border = TRUE)+
  facet_wrap(.~variable, scale = "free")+
  labs(x = "Timepoint", y = "Cell proportion")
dev.off()

################## Figure Supplementary 4 ##################

cor_method1_clocks <- cor(t(avg_method1))
cor_method2_clocks <- cor(t(avg_method2))
cor_method3_clocks <- cor(t(avg_method3))

cor_method1_clocks_grouped <- reshape2::acast(reshape2::melt(cor_method1_clocks) %>% mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>% mutate(group1 = sapply(Var1, get_group)) %>% mutate(group2 = sapply(Var2, get_group)) %>% 
                                                group_by(group1,group2) %>% summarise(mean = mean(value)), group1~group2, value.var = "mean")

cor_method2_clocks_grouped <- reshape2::acast(reshape2::melt(cor_method2_clocks) %>% mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>% mutate(group1 = sapply(Var1, get_group)) %>% mutate(group2 = sapply(Var2, get_group)) %>% 
                                                group_by(group1,group2) %>% summarise(mean = mean(value)), group1~group2, value.var = "mean")

cor_method3_clocks_grouped <- reshape2::acast(reshape2::melt(cor_method3_clocks) %>% mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>% mutate(group1 = sapply(Var1, get_group)) %>% mutate(group2 = sapply(Var2, get_group)) %>% 
                                                group_by(group1,group2) %>% summarise(mean = mean(value)), group1~group2, value.var = "mean")

col_fun <- colorRamp2(seq(-1,1,0.2), c(colors$blue[5:1],"white",colors$red[1:5]))

ht1_clock <- Heatmap(
  cor_method1_clocks_grouped,
  rect_gp = gpar(col = "black", lwd = 1, lty = 3),
  border = TRUE,
  name = "R",
  column_title = "Age difference",
  row_dend_side = "right",
  column_dend_side = "top",
  row_names_side = "left",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_heatmap_legend = FALSE,
  col = col_fun
)

ht2_clock <- Heatmap(
  cor_method2_clocks_grouped,
  rect_gp = gpar(col = "black", lwd = 1, lty = 3),
  border = TRUE,
  name = "R",
  column_title = "Age acceleration",
  row_dend_side = "right",
  column_dend_side = "top",
  row_names_side = "left",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_heatmap_legend = FALSE,
  col = col_fun
)

ht3_clock <- Heatmap(
  cor_method3_clocks_grouped,
  rect_gp = gpar(col = "black", lwd = 1, lty = 3),
  border = TRUE,
  name = "R",
  column_title = "Intrinsic age acceleration",
  row_dend_side = "right",
  column_dend_side = "top",
  row_names_side = "left",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_heatmap_legend = TRUE,
  col = col_fun
)

p1_clock <- grid.grabExpr({draw(ht1_clock)})
p2_clock <- grid.grabExpr({draw(ht2_clock)})
p3_clock <- grid.grabExpr({draw(ht3_clock)})

library(cowplot)
pdf(file = "./output/Supplementary_Figure5.pdf", width=24.42, height=8.12)
plot_grid(p1_clock,p2_clock,p3_clock, nrow = 1, rel_widths = c(1,1,1.1))
dev.off()

################## Figure Supplementary 6 ##################

all_estimates_temp$group <- sapply(all_estimates_temp$clock, get_group)
per_group <- all_estimates_temp %>% group_by(group, timepoint, system) %>% summarise(ba = mean(ba))

aa_group <- per_group %>% filter(system=="Age acceleration"&timepoint=="L-45") %>% ungroup %>% dplyr::select(group, system, timepoint, ba) %>% set_names("group","system","t1","pre") %>% 
left_join(per_group %>% filter(system=="Age acceleration"&timepoint!="L-45") %>% ungroup %>% dplyr::select(group, timepoint, ba) %>% set_names("group","t2","post")) %>%
  mutate(dif = post-pre) 

iaa_group <- per_group %>% filter(system=="Intrinsic age acceleration"&timepoint=="L-45") %>% ungroup %>% dplyr::select(group, system, timepoint, ba) %>% set_names("group","system","t1","pre") %>% 
  left_join(per_group %>% filter(system=="Intrinsic age acceleration"&timepoint!="L-45") %>% ungroup %>% dplyr::select(group, timepoint, ba) %>% set_names("group","t2","post")) %>%
  mutate(dif = post-pre) 

group_comparison <- rbind(aa_group,iaa_group)
writexl::write_xlsx(group_comparison, path = "./output/Supplementary_Table_5.xlsx")

pdf(file = "./output/Supplementary_Figure6.pdf", width=11.77, height=5.27)
ggplot(per_group, aes(timepoint, ba, group = group, fill = group))+
  geom_hline(yintercept = 0, color = "black", lty = 2)+
  geom_rect(
    data = group_colors_ax,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = c('grey80','white','grey80','grey80','white','grey80','grey80','white','grey80'),
    inherit.aes = FALSE, alpha = 0.5
  ) +
  geom_line(color = "black")+
  facet_grid(.~system)+
  geom_point(shape = 21, size = 3)+
  scale_fill_manual(values = c(colors$tone[3],colors$grey[3],colors$olive[3],colors$teal[3],colors$purple[3],colors$yellow[3],colors$skin[3],colors$blue[3]))+
  guides(color = "none")+
  theme_pubr(border = TRUE)+
  labs(x = "Time point", y = "Average value", fill = "Category")
dev.off()

