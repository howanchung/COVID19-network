family_size <- as.data.table(table(COVID19$family_index2))
setnames(family_size, c('family_index', 'family_size'))
family_size$family_index <- as.numeric(family_size$family_index)

mean(family_size$family_size)

family_size_dist <- family_size %>% 
  group_by(family_size) %>%
  summarise(N = n()) %>%
  as.data.table()

family_size_dist[, prop := N / sum(N)][, family_size := factor(family_size)]

avg_family_size <- family_size %>%
  summarise(mean = round(mean(family_size), 3),
            se = round(sd(family_size) / sqrt(n()), 3))

ggplot(family_size_dist, aes(x = fct_inorder(family_size), y = prop)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 6, family = FONTS) +
    scale_y_continuous(breaks = seq(0, 0.8, 0.2),
                       labels = scales::percent_format(accuracy = 2),
                       limits = c(0, 0.8)) +
    labs(x = 'Family size', y = 'Proportion') +
    # annotation_custom(tableGrob(avg_family_size, rows = NULL, 
    #                             cols = c('平均家庭大小', '标准误'),
    #                             ttheme_default(base_size = 22,
    #                                            padding = unit(c(7, 5), 'mm'),
    #                                            base_family = FONTS)),
    #                   xmin = 6.5, ymin = 0.50) +
    theme_bw() +
    theme(text = element_text(size = 30, family = FONTS),
          legend.position = 'left',
          panel.border = element_rect(size = 1.1))

ggsave(filename = "figure/family_size_distribution.pdf", width = 15, height = 6, device = cairo_pdf)

COVID19 <- as.data.table(COVID19)
COVID19[, family_index := NULL]
setkey(COVID19, family_index2)

family_age_distribution <- list()
unique_size <- unique(family_size$family_size)
for(i in 1:length(unique_size)){
  family_age_distribution[[i]] <- COVID19[family_index2 %in% 
                                       family_size[family_size == unique_size[i], 
                                                   family_index], 
                                     c('age_year', 'family_index2')]
  family_age_distribution[[i]][, family_size := unique_size[i]]
}
family_age_distribution <- rbindlist(family_age_distribution)
  

# family_age_distribution[, age_group := cut(age_year, breaks = c(seq(0, 65, 5), 100), 
#                                            right = FALSE)]

family_age_distribution %>%
  group_by(family_size, family_index2) %>%
  summarise(concat = paste("(", paste(sort(age_year), collapse = ", "), ")", 
                           collapse="", sep = ""), 
            .groups="drop") %>%
  group_by(family_size) %>%
  summarise(concat = paste(concat, collapse = ", "), .groups="drop") %>%
  as.data.table() %>%
  write.csv("data/family_age_distribution_concat.csv")

family_age_distribution[, age_group := cut(age_year, breaks = c(seq(0, 65, 5), 100), 
                                           right = FALSE)]
write.csv(family_age_distribution, file = 'family_age_distribution.csv')
write.csv(family_size, file = 'family_size.csv')


a <- cut(COVID19$age_year, breaks = c(seq(0, 65, 5), 100), 
    right = FALSE, labels = FALSE)
b <- table(a)
b/sum(b)