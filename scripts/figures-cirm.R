# specify libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# load output file
outputs <- read.csv("results/results_all_contacts.csv", 
                    stringsAsFactors = FALSE)

# plot results over n_infectors
# x-axis n_infectors (5, 25, 50, 100)
# y-axis estimated median and SD of mean and true mean...we don't have logmean_mean...

# accuracy
# add logmean_mean to df
# violin/boxplot of logmean_mean across replicates 

# precision
# add lower_95CI and upper_95CI?
# plot distribution of width of CIs?
# violin/boxplot of logmean_sd across replicates? 


# make n_infectors and alpha_mean factors
plot_df <- outputs %>% 
  mutate(n_infectors = as.factor(n_infectors),
         alpha_mean = as.factor(alpha_mean))

p <- ggplot(data = plot_df,
            mapping = aes (x=n_infectors,
                           y = logmean_median,
                           ymin = logmean_median - logsd_median,
                           ymax = logmean_median + logsd_median,
                           colour = alpha_mean)) +
  geom_point(size = 4,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(size = 1,
                width = 0.25,
                position = position_dodge(width = 0.5)) + 
  geom_hline(aes(yintercept = true_logmean, 
                 linetype = "true mean"), color = "black") +
  theme_bw() + 
  theme(
    axis.line = element_line(colour = "grey50"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        plot.title = element_text(size = 22,
                                  hjust = 0.5)) +
  scale_linetype_manual(name = "", 
                        values = c("true mean" = "dashed")) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(name = 'median (log mean)',
                     limits = c(0, 3)) 
  #scale_x_discrete(name='') +
  #labs(title='')
  
 p 

pdf('outputs/test_plot1.pdf',
    width = 5,
    height = 4)
plot(p)
dev.off()

