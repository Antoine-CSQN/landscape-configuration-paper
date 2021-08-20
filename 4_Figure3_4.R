# Script to realize the Figure of the paper
Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252")

## Figure 3
# Figure 3a: Optimization plans
opti_plan_df <- allCorr %>% filter(param %in% c("NO3","TP"),
                                   source.name == "sau_dont_bh_2010_2019",
                                   date == "2001-01-04")
opti_plan_df$label_param <- if_else(opti_plan_df$param == "NO3","NO[3]^{\"-\"}","TP")

max_by_param <- opti_plan_df %>% 
  arrange(-R) %>%
  group_by(param) %>%
  slice(1)

max10pc_by_param <- opti_plan_df %>%
  group_by(param) %>%
  arrange(-R) %>%
  slice(1:as.integer(n()/10))

fig3a <- ggplot(data = opti_plan_df,
                aes(x = coef.fls, y = coef.facc, fill = R))+
  geom_tile() + 
  scale_fill_viridis_b(breaks = seq(0.3,0.8,0.05), limits = c(0.2,0.85), name ="ρ") + 
  # Marquage des 10% des meilleures corrélations et de la meilleure corrélation
  geom_point(data = max10pc_by_param, pch = 16, color = "black", size = 0.5) + #0.25
  geom_point(data = max_by_param, pch = 0, color = "red", size = 2) + #1
  # One facet per parameter
  facet_wrap(.~label_param, labeller = label_parsed) +
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal() + 
  
  labs( x = "Coefficient (b) applied to inverse distance to stream/ditch (FLS)",
        y = "Coefficient (a) applied\nto flow accumulation (FAcc)") +
  
  theme_bw() + 
  theme(legend.position = "bottom") + 
  theme(legend.key.height =  unit(10, "pt")) + 
  theme(legend.key.width = unit(80, "pt")) +
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black"))
fig3a
ggsave("fig3a.tiff", plot = fig3a, width = 174, units = "mm", dpi = 600)  

# Figure 3b: Show the scatter plot for median TP with LCI(0,0) and LCI(1.4, 2.2)
paraM <- "TP"
sourceName <- "sau_dont_bh_2010_2019"
# Median TP concentration values
paraM.median <- geochem.median %>% dplyr :: select(subcat.name = p.sampling, !!paraM)
# Value of mean lumped value for meso scale catchment (= pre procressing area)
lumped.meso <- preprocess.results$mean.index[which(preprocess.results$source.name == sourceName &
                                                     preprocess.results$coef.facc == 0 &
                                                     preprocess.results$coef.fls == 0)]
# %of param in each subcatchment
lumped.sourceName <- res.df %>% filter(coef.facc == 0,
                                       coef.fls == 0,
                                       source.name == sourceName) %>%
  mutate(index.value = index.value) %>%
  dplyr :: select(subcat.name, LDI = index.value) %>%
  mutate(metric = "I(0,0)")

# Best coef for mixed, paraM, sourceName, median paraM (date = 2001-01-04)
paraM.best.coef <- bestCor %>% filter(date == '2001-01-04',
                                      source.name == !!sourceName,
                                      metric == "mixed",
                                      param == paraM)

paraM.lumped.coef <- bestCor %>% filter(date == '2001-01-04',
                                        source.name == !!sourceName,
                                        metric == "lumped",
                                        param == paraM)

mixed.sourceName <- res.df %>% filter(coef.facc == paraM.best.coef$coef.facc,
                                      coef.fls == paraM.best.coef$coef.fls,
                                      source.name == sourceName) %>%
  dplyr :: select(subcat.name,LDI = index.value)  %>%
  mutate(metric = "I(1.4,2.2)")

# Prepare df for plotting
paraM.optim <- rbind.data.frame(lumped.sourceName, mixed.sourceName) %>%
  inner_join(paraM.median, by = "subcat.name")
#Prepared ranked data
lumped.sourceName.ranked <- lumped.sourceName
lumped.sourceName.ranked$LDI <- rank(lumped.sourceName.ranked$LDI)
lumped.sourceName.ranked$metric <- "rank_I(0,0)"

mixed.sourceName.ranked <- mixed.sourceName
mixed.sourceName.ranked$LDI <- rank(mixed.sourceName.ranked$LDI)
mixed.sourceName.ranked$metric <- "rank_I(1.4,2.2)"

paraM.median.ranked <- paraM.median
paraM.median.ranked[,paraM] <-rank(paraM.median.ranked[,paraM])

# Helpr function to have exactly the name number of digits
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# join the data
paraM.optim.ranked <- rbind.data.frame(lumped.sourceName.ranked, mixed.sourceName.ranked) %>%
  inner_join(paraM.median.ranked, by = "subcat.name") 

paraM.optim <-  rbind.data.frame(paraM.optim, paraM.optim.ranked) 

paraM.optim.w <- paraM.optim %>% pivot_wider(names_from = metric, values_from =  c(!!paraM, LDI))
param.optim.labels <- paraM.optim.w %>% arrange(`TP_I(1.4,2.2)`)
param.optim.labels <- param.optim.labels[c(1:3,16:19),]

g.compo <- ggplot(paraM.optim.w, aes(`LDI_rank_I(0,0)`, `TP_rank_I(0,0)`)) + 
  geom_point(color = "#5bc0de", size = 2) + 
  geom_text_repel(data =  param.optim.labels,
                  aes(label = paste0(specify_decimal(`LDI_I(0,0)`,2)," : ", specify_decimal(`TP_I(0,0)`,3))),
                  size = 3.5,
                  nudge_x = c(3,3,3,-6,-6,-6,-6)) + 
  labs(title = paste0("ρ = ", specify_decimal(paraM.lumped.coef$R, 2),
                      ', p-val = ', specify_decimal(paraM.lumped.coef$p.value, 3)),
       x = "Rank of LCI(0, 0)\nLandscape Composition",
       y = "Rank of median [TP]") + 
  coord_equal() + 
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black"))
g.compo

g.distri <- ggplot(paraM.optim.w, aes(`LDI_rank_I(1.4,2.2)`, `TP_rank_I(1.4,2.2)`)) + 
  geom_point(color = "#337ab7", size = 2) + 
  geom_text_repel(data =  param.optim.labels,
                  aes(label = paste0(specify_decimal(`LDI_I(1.4,2.2)`,2)," : ", specify_decimal(`TP_I(1.4,2.2)`, 3))),
                  size = 3.5,
                  nudge_x = c(4,3,4,-6,-6,-6,-6)) + 
  labs(title = paste0("ρ = ", specify_decimal(paraM.best.coef$R, 2),
                      ', p-val = ', format(paraM.best.coef$p.value, digits = 3, scientific = TRUE)),
       x = "Rank of LCI(1.4, 2.2)\nLandscape Configuration",
       y = "Rank of median [TP]") + 
  coord_equal() + 
  theme_bw() + 
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black"))

fig3b <- g.compo + g.distri  #+ 
# plot_annotation(title = 'Rank correlation between LDI(0,0), LDI(1.4,2.2) and median [TP] ',
#                 #caption = "LDI value : median [TP] concentration (ug)",
#                 theme = theme(plot.title = element_text(size = 10))&
#                   theme(plot.caption = element_text(size = 7)))
fig3b
ggsave(filename = "fig3b.tiff", plot = fig3b, width = 174, units = "mm", dpi = 600)

fig3a/fig3b + plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave2(filename = "fig3.tiff", width = 174, units = "mm", dpi = 600)

## Figure 4
var.name <- "TP"
source.name.label <- "Agricultural surface including riparian buffer strips"
R.graph.data <- bestCor %>% filter (param == var.name) %>% 
  filter(metric %in% c("lumped", "mixed")) %>%
  filter(source.nut ==  "sau_dont_bh") %>%
  filter(grouped == FALSE) %>%
  dplyr :: select(date, metric, R) %>%
  mutate(metric = ifelse(metric == "lumped", "LCI(0,0):\nLandscape Composition", "Optimal LCI(a,b):\nLandscape Configuration"))

pval.graph.data <- bestCor %>% filter (param == var.name) %>% 
  filter(metric %in% c("lumped", "mixed")) %>%
  filter(source.nut ==  "sau_dont_bh") %>%
  filter(grouped == FALSE) %>%
  mutate(log10p.value = -log10(p.value)) %>%
  mutate(log10p.value = ifelse(log10p.value>6, 6, log10p.value)) %>%
  dplyr :: select(date, metric, log10p.value) %>%
  mutate(metric = ifelse(metric == "lumped", "LCI(0,0):\nLandscape Composition", "Optimal LCI(a,b):\nLandscape Configuration"))

coef.graph.data <- bestCor %>% filter (param == var.name) %>% 
  filter(metric == "mixed") %>%
  filter(source.nut ==  "sau_dont_bh") %>%
  filter(grouped == FALSE) %>%
  dplyr :: select(date, coef.facc, coef.fls) 
coef.graph.data <- melt(coef.graph.data, id.vars = "date") %>%
  rename(coefficient = variable) %>%
  mutate(coefficient = ifelse(coefficient == "coef.facc", "a", "b"))

weight.graph.data <- bestCor %>% filter (param == var.name) %>% 
  filter(metric %in% c("lumped","mixed")) %>%
  filter(source.nut ==  "sau_dont_bh") %>%
  filter(grouped == FALSE) %>%
  dplyr :: select(date, metric, pc.05) %>%
  mutate(metric = ifelse(metric == "lumped", "LCI(0,0):\nLandscape Composition", "Optimal LCI(a,b):\nLandscape Configuration"))

R.graph <- ggplot(data = R.graph.data,
                  aes (x = date, y = R, fill = metric, color = metric)) + 
  scale_color_manual(values = c("#5bc0de", "#337ab7")) + 
  geom_rect(xmin = as.Date("2018-07-17"), xmax = as.Date("2018-11-27"), ymin = -2, ymax = 2, 
            color = NA, fill = "grey90") + 
  geom_rect(xmin = as.Date("2019-07-02"), xmax = as.Date("2019-07-30"), ymin = -2, ymax = 2, 
            color = NA, fill = "grey90") +
  geom_line(size = 1) + geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  scale_x_date(breaks = '2 month', date_labels = "%b-%y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  scale_y_continuous(limits = c(-1,1)) +
  labs(title = NULL, x=NULL, y= "ρ") +
  theme(legend.position="bottom") 
R.graph 

median.a <- coef.graph.data %>% 
  filter(date %in% dates.flowing$date.sampling,
         coefficient == "a") %>%
  summarise(a = median(value))
median.b <- coef.graph.data %>% 
  filter(date %in% dates.flowing$date.sampling,
         coefficient == "b") %>%
  summarise(b = median(value))

coef.graph <- ggplot(data = coef.graph.data,
                     aes (x = date, y = value, fill = coefficient, color = coefficient)) +
  scale_color_manual(values = c("#33cc33", "#336600")) + 
  geom_rect(xmin = as.Date("2018-07-17"), xmax = as.Date("2018-11-27"), ymin = -1, ymax = 6, 
            color = NA, fill = "grey90") + 
  geom_rect(xmin = as.Date("2019-07-02"), xmax = as.Date("2019-07-30"), ymin = -1, ymax = 6, 
            color = NA, fill = "grey90") +
  geom_line(size = 1) + geom_point(size = 2) +
  scale_x_date(breaks = '2 month', date_labels = "%b-%y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  scale_y_continuous(limits = c(0,4)) +
  labs(x=NULL, y="Optimal (a, b)", title = NULL) +
  theme(legend.position="bottom")
coef.graph

pval.graph <- ggplot(data = pval.graph.data,
                     aes (x = date, y = log10p.value, fill = metric, color = metric)) + 
  scale_color_manual(values = c("#5bc0de", "#337ab7")) + 
  geom_rect(xmin = as.Date("2018-07-17"), xmax = as.Date("2018-11-27"), ymin = -1, ymax = 7, 
            color = NA, fill = "grey90") + 
  geom_rect(xmin = as.Date("2019-07-02"), xmax = as.Date("2019-07-30"), ymin = -1, ymax = 7, 
            color = NA, fill = "grey90") +
  geom_line(size = 1) + geom_point(size = 2) +
  scale_x_date(breaks = '2 month', date_labels = "%b-%y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  scale_y_continuous(limits = c(0,6)) +
  geom_hline(yintercept = -log10(5e-2), linetype = "dashed")+
  theme(legend.position="bottom") +
  labs(x=NULL, y="-log10(p-value)", title = NULL)
pval.graph

weight.graph <- ggplot(data = weight.graph.data,
                       aes (x = date, y = pc.05, fill = metric, color = metric)) + 
  scale_color_manual(values = c("#5bc0de", "#337ab7")) + 
  geom_rect(xmin = as.Date("2018-07-17"), xmax = as.Date("2018-11-27"), ymin = -1, ymax = 2, 
            color = NA, fill = "grey90") + 
  geom_rect(xmin = as.Date("2019-07-02"), xmax = as.Date("2019-07-30"), ymin = -1, ymax = 2, 
            color = NA, fill = "grey90") +
  geom_line(size = 1) + geom_point(size = 2) +
  scale_x_date(breaks = '2 month', date_labels = "%b-%y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.title = element_text(size = 10, face = "plain", colour = "black"),
        legend.text = element_text(size = 10, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 10, face = "plain", colour = "black"),    
        axis.title.y = element_text(size = 10, face = "plain", colour = "black"),    
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black"),
        strip.text.x = element_text(size = 10, face = "plain", colour = "black" ),
        strip.text.y = element_text(size = 10, face = "plain", colour = "black")) + 
  theme(plot.title = element_text(size=9)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed")+
  theme(legend.position="bottom") +
  labs(x=NULL, y="Weight of the 5% top weighted pixels")
weight.graph 


graph.pat <- (R.graph | pval.graph) / (coef.graph | weight.graph) + 
  plot_annotation(tag_levels = "a", tag_suffix = ")") + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Fig4.tiff", plot = graph.pat, width = 174, units = "mm", dpi = 600)
ggsave("Fig4.pdf", plot = graph.pat, width = 174, units = "mm", dpi = 600)

## Figure 5
load("F:/doctorat/R/IDW3/results/4_These.RData")
