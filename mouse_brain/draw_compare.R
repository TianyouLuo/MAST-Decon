library(tidyverse)
library(scales)

folder = "/pine/scr/t/i/tianyou/ST/mouse_brain_cell2location/comparison/"
f = list.files(folder)
dat = list()

for (i in 1:length(f)){
  comparison = read_csv(file.path(folder, f[i]))
  comparison_long = comparison %>%
    mutate(`MAST-Decon` = pmax(`low-smoothing`, `medium-smoothing`, `high-smoothing`)) %>%
    mutate(slice = str_sub(f[i],-13,-5)) %>%
    select(slice, nUMI, RCTD, `MAST-Decon`) %>%
    pivot_longer(cols = c("RCTD", "MAST-Decon"), names_to = "method", values_to = "ARI")
  dat[[i]] = comparison_long
}
dat = bind_rows(dat)

fig = dat %>%
  group_by(nUMI, method) %>%
  summarise(mean_ARI = mean(ARI), sd_ARI = sd(ARI)) %>%
  ggplot(aes(x = nUMI, y = mean_ARI)) +
  geom_point(aes(color = method), position = position_dodge(0.1)) +
  geom_line(aes(color = method), position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin = mean_ARI - sd_ARI, ymax = mean_ARI + sd_ARI, color = method), 
                width = 0.1, position = position_dodge(0.1)) +
  labs(x = "mean UMI count",
       y = "Adjusted Rand Index (ARI)") +
  theme(legend.position = c(1,0),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0.4)),
        text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave("comparison_max4.png", fig, width = 9, height = 12, units = "cm")


comparison %>%
  mutate(`MAST-Decon` = pmax(`low-smoothing`, `medium-smoothing`, `high-smoothing`)) %>%
  pivot_longer(cols = c("RCTD", "MAST-Decon", "Stereoscope"), names_to = "method", values_to = "ARI") %>%
  ggplot(aes(x = nUMI, y = ARI)) +
  geom_point(aes(color = method)) +
  geom_line(aes(color = method)) +
  labs(x = "mean UMI count",
       y = "Adjusted Rand Index (ARI)") +
  theme(legend.position = c(1,0),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0.4)),
        text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggsave("comparison_max4.png", width = 9, height = 12, units = "cm")

