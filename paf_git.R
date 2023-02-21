# Generating 95%CIs for PAFs using Monte Carlo simulation
## Load ----------------------------------------------------------------

# Set working directory
setwd("~/.")

# Load packages
pacman::p_load(dplyr, readr, forcats,
        ggplot2, ggpubr, ggbreak, ggimage,
        grid, gridExtra, stringr)

# Load RRs and 95%CIs
sim_input <- read_csv(file ="simulation_input_git1.csv")

# disable scientific notation
options(scipen = 999)

## Log  --------------------------------------------------------------------

log_input <- sim_input

# Log to get symmetrical confidence interval and standard error
log_input <- log_input %>%
  mutate(logrr = log(rr)) %>%
  mutate(loglci = log(rr_lci)) %>%
  mutate(loguci = log(rr_uci)) %>%
  mutate(se = (logrr-loglci)/qnorm(0.975))

## Simulating RRs ----------------------------------------------------------

# Create function to generate 10000 theoretical observations of the log RR for each rr_outcome
simulate_rr <- function(rownum) {
  set.seed(1)
  sim_rr <- exp(rnorm(10000, mean = log_input$logrr[rownum], sd = log_input$se[rownum]))
  sim_rr <- as.data.frame(sim_rr)
}

# Simulate RR for all outcomes
sim_rr_1 <- simulate_rr(1)
sim_rr_2 <- simulate_rr(2)
sim_rr_3 <- simulate_rr(3)
sim_rr_4 <- simulate_rr(4)
sim_rr_5 <- simulate_rr(5)
sim_rr_6 <- simulate_rr(6)
sim_rr_7 <- simulate_rr(7)
sim_rr_8 <- simulate_rr(8)
sim_rr_9 <- simulate_rr(9)
sim_rr_10 <- simulate_rr(10)
sim_rr_11 <- simulate_rr(11)
sim_rr_12 <- simulate_rr(12)
sim_rr_13 <- simulate_rr(13)
sim_rr_14 <- simulate_rr(14)
sim_rr_15 <- simulate_rr(15)
sim_rr_16 <- simulate_rr(16)
sim_rr_17 <- simulate_rr(17)
sim_rr_18 <- simulate_rr(18)
sim_rr_19 <- simulate_rr(19)
sim_rr_20 <- simulate_rr(20)
sim_rr_21 <- simulate_rr(21)
sim_rr_22 <- simulate_rr(22)
sim_rr_23 <- simulate_rr(23)
sim_rr_24 <- simulate_rr(24)
sim_rr_25 <- simulate_rr(25)

## Simulating population proportion ----------------------------------------

# Log to get symmetrical confidence interval and standard error
log_input <- log_input %>%
  mutate(logp = log(p)) %>%
  mutate(logplci = log(p_lci)) %>%
  mutate(logpuci = log(p_uci)) %>%
  mutate(pse = (logp-logplci)/qnorm(0.975))

# create standard error for multiple IG group's p, where margin of error (moe) is 10% of p, and where se = moe/1.96 (at the 95% level)
# keep this to later check the point estimate of the Monte Carlo vs PERT
log_input <- log_input %>% 
  mutate(pse = ifelse(is.na(pse), p*0.1/1.96, pse))

# Create function to generate 10000 theoretical observations of p for each rr_outcome
simulate_p_nolog <- function(rownum) {
  set.seed(1)
  sim_p <- rnorm(10000, mean = log_input$p[rownum], sd = log_input$pse[rownum])
  sim_p <- as.data.frame(sim_p)
}

simulate_p <- function(rownum) {
  set.seed(1)
  sim_p <- exp(rnorm(10000, mean = log_input$logp[rownum], sd = log_input$pse[rownum]))
  sim_p <- as.data.frame(sim_p)
}

# Simulate p for all outcomes
sim_p_1 <- simulate_p_nolog(1)
sim_p_2 <- simulate_p(2)
sim_p_3 <- simulate_p(3)
sim_p_4 <- simulate_p(4)
sim_p_5 <- simulate_p(5)
sim_p_6 <- simulate_p(6)
sim_p_7 <- simulate_p(7)
sim_p_8 <- simulate_p(8)
sim_p_9 <- simulate_p(9)
sim_p_10 <- simulate_p(10)
sim_p_11 <- simulate_p(11)
sim_p_12 <- simulate_p(12)
sim_p_13 <- simulate_p(13)
sim_p_14 <- simulate_p(14)
sim_p_15 <- simulate_p(15)
sim_p_16 <- simulate_p(16)
sim_p_17 <- simulate_p(17)
sim_p_18 <- simulate_p(18)
sim_p_19 <- simulate_p(19)
sim_p_20 <- simulate_p(20)
sim_p_21 <- simulate_p(21)
sim_p_22 <- simulate_p(22)
sim_p_23 <- simulate_p(23)
sim_p_24 <- simulate_p(24)
sim_p_25 <- simulate_p(25)

## Calculate PAFs ----------------------------------------------------------

# Create function to calculate PAFs
simulate_paf <- function(sim_rr, sim_p) {
  pafs <- cbind(sim_rr, sim_p) %>%
    mutate(paf = sim_p*(sim_rr-1) / (sim_p*(sim_rr-1)+1))
  quantitle <- as.data.frame(quantile(pafs$paf, probs = c(0.5, 0.025, 0.975)))
}

pafs_1 <- as.data.frame(t(simulate_paf(sim_rr_1, sim_p_1))) %>% mutate(rr_outcome = log_input$rr_outcome[1]) 
pafs_2 <- as.data.frame(t(simulate_paf(sim_rr_2, sim_p_2))) %>% mutate(rr_outcome = log_input$rr_outcome[2])
pafs_3 <- as.data.frame(t(simulate_paf(sim_rr_3, sim_p_3))) %>% mutate(rr_outcome = log_input$rr_outcome[3]) 
pafs_4 <- as.data.frame(t(simulate_paf(sim_rr_4, sim_p_4))) %>% mutate(rr_outcome = log_input$rr_outcome[4]) 
pafs_5 <- as.data.frame(t(simulate_paf(sim_rr_5, sim_p_5))) %>% mutate(rr_outcome = log_input$rr_outcome[5]) 
pafs_6 <- as.data.frame(t(simulate_paf(sim_rr_6, sim_p_6))) %>% mutate(rr_outcome = log_input$rr_outcome[6]) 
pafs_7 <- as.data.frame(t(simulate_paf(sim_rr_7, sim_p_7))) %>% mutate(rr_outcome = log_input$rr_outcome[7])
pafs_8 <- as.data.frame(t(simulate_paf(sim_rr_8, sim_p_8))) %>% mutate(rr_outcome = log_input$rr_outcome[8]) 
pafs_9 <- as.data.frame(t(simulate_paf(sim_rr_9, sim_p_9))) %>% mutate(rr_outcome = log_input$rr_outcome[9])
pafs_10 <- as.data.frame(t(simulate_paf(sim_rr_10, sim_p_10))) %>% mutate(rr_outcome = log_input$rr_outcome[10]) 
pafs_11 <- as.data.frame(t(simulate_paf(sim_rr_11, sim_p_11))) %>% mutate(rr_outcome = log_input$rr_outcome[11])
pafs_12 <- as.data.frame(t(simulate_paf(sim_rr_12, sim_p_12))) %>% mutate(rr_outcome = log_input$rr_outcome[12])
pafs_13 <- as.data.frame(t(simulate_paf(sim_rr_13, sim_p_13))) %>% mutate(rr_outcome = log_input$rr_outcome[13]) 
pafs_14 <- as.data.frame(t(simulate_paf(sim_rr_14, sim_p_14))) %>% mutate(rr_outcome = log_input$rr_outcome[14]) 
pafs_15 <- as.data.frame(t(simulate_paf(sim_rr_15, sim_p_15))) %>% mutate(rr_outcome = log_input$rr_outcome[15]) 
pafs_16 <- as.data.frame(t(simulate_paf(sim_rr_16, sim_p_16))) %>% mutate(rr_outcome = log_input$rr_outcome[16]) 
pafs_17 <- as.data.frame(t(simulate_paf(sim_rr_17, sim_p_17))) %>% mutate(rr_outcome = log_input$rr_outcome[17]) 
pafs_18 <- as.data.frame(t(simulate_paf(sim_rr_18, sim_p_18))) %>% mutate(rr_outcome = log_input$rr_outcome[18]) 
pafs_19 <- as.data.frame(t(simulate_paf(sim_rr_19, sim_p_19))) %>% mutate(rr_outcome = log_input$rr_outcome[19])
pafs_20 <- as.data.frame(t(simulate_paf(sim_rr_20, sim_p_20))) %>% mutate(rr_outcome = log_input$rr_outcome[20])
pafs_21 <- as.data.frame(t(simulate_paf(sim_rr_21, sim_p_21))) %>% mutate(rr_outcome = log_input$rr_outcome[21])
pafs_22 <- as.data.frame(t(simulate_paf(sim_rr_22, sim_p_22))) %>% mutate(rr_outcome = log_input$rr_outcome[22])
pafs_23 <- as.data.frame(t(simulate_paf(sim_rr_23, sim_p_23))) %>% mutate(rr_outcome = log_input$rr_outcome[23])
pafs_24 <- as.data.frame(t(simulate_paf(sim_rr_24, sim_p_24))) %>% mutate(rr_outcome = log_input$rr_outcome[24])
pafs_25 <- as.data.frame(t(simulate_paf(sim_rr_25, sim_p_25))) %>% mutate(rr_outcome = log_input$rr_outcome[25])

# Combine all outcomes
pafs_all <- rbind(pafs_1, pafs_2, pafs_3, pafs_4, pafs_5,
                  pafs_6, pafs_7, pafs_8, pafs_9, pafs_10,
                  pafs_11, pafs_12, pafs_13, pafs_14, pafs_15,
                  pafs_16, pafs_17, pafs_18, pafs_19, pafs_20,
                  pafs_21, pafs_22, pafs_23, pafs_24, pafs_25) %>%
  rename(paf = `50%`) %>%
  rename(paf_lci = `2.5%`) %>%
  rename(paf_uci = `97.5%`)


## Replace with PERT for all-cause estimate -------------------------

d <- read.csv("simulation_input_git2.csv")

# 10000 simulations
B <- 1e4

# rate ratios
d$logrr <- log(d$rr)
d$se <- (d$logrr - log(d$rr_lci)) / qnorm(0.975)
rrs <- mapply(rnorm, mean = d$logrr, sd = d$se, n = list(B), SIMPLIFY = T)
rrs <- t(exp(rrs))

# populations using modified PERT distribution
pops <- mapply(rpert, min = d$pop_lower, mode = d$pop_mean, max = d$pop_upper, n = list(B), SIMPLIFY = T)
pops <- t(pops)

# calculate PAFs
pafs <- pops * (rrs-1) / (pops * (rrs-1)+1)
paf_results <- t(apply(pafs, 1, quantile, probs = c(0.5, 0.025, 0.975))) 
paf_results <- as.data.frame(paf_results)

# combine with rest of pafs
pafs_all <- pafs_all %>% 
  mutate(paf = replace(paf, rr_outcome == "All-cause (inclusion health)", paf_results$`50%`)) %>% 
  mutate(paf_lci = replace(paf_lci, rr_outcome == "All-cause (inclusion health)", paf_results$`2.5%`)) %>% 
  mutate(paf_uci = replace(paf_uci, rr_outcome == "All-cause (inclusion health)", paf_results$`97.5%`))

## Calculate total number attributable --------------------------------------------------

# Select total number of cases from original input
cases <- log_input %>%
  select(rr_outcome, totalnum, pse)

# Calculate number attributable using PAF
pafs_all <- pafs_all %>%
  left_join(cases, by = c("rr_outcome")) %>%
  mutate(numattrib = totalnum*paf) %>%
  mutate(numattrib_lci = totalnum*paf_lci) %>%
  mutate(numattrib_uci = totalnum*paf_uci) %>%
  select(-totalnum)

# join other columns of info
pafs_all <- pafs_all %>%
  left_join(sim_input, by = c("rr_outcome"))

# turn PAFs into % for paper, and combine 95%CI columns for ease of reporting
pafs_all_percent <- pafs_all %>%
  mutate(numattrib = round(numattrib, 0)) %>%
  mutate(numattrib_lci = round(numattrib_lci, 0)) %>%
  mutate(numattrib_uci = round(numattrib_uci, 0)) %>%
  mutate(paf = round(paf*100, 2)) %>% 
  mutate(paf_lci = round(paf_lci*100, 2)) %>% 
  mutate(paf_uci = round(paf_uci*100, 2)) %>% 
  mutate(paf_ci = paste0(paf,'(',paf_lci,'-',paf_uci,')')) %>%
  mutate(numattrib_ci = paste0(numattrib,'(',numattrib_lci,'-',numattrib_uci,')')) %>%
  mutate(rr_ci = paste0(rr,'(',rr_lci,'-',rr_uci,')')) 
  
# Save
write.csv(pafs_all_percent, file = "simulation_output.csv", row.names = F)

## Prepare outputs for plots -------------------------------------------------------

stacked_totalnum <- pafs_all %>%
  select(rr_outcome, totalnum, numattrib_lci, numattrib_uci, paf, paf_lci, paf_uci) %>%
  rename(numattrib = totalnum) %>%
  mutate(group = 'All deaths') %>%
  mutate(numattrib = numattrib/1000) %>% # turn into thousands
  mutate(numattrib_lci = numattrib_lci/1000) %>%
  mutate(numattrib_uci = numattrib_uci/1000)

stacked_numattrib <- pafs_all %>%
  select(rr_outcome, numattrib, numattrib_lci, numattrib_uci, paf, paf_lci, paf_uci) %>%
  mutate(group = 'Deaths attributable to inclusion health group(s)') %>%
  mutate(numattrib = numattrib/1000) %>%
  mutate(numattrib_lci = numattrib_lci/1000) %>%
  mutate(numattrib_uci = numattrib_uci/1000)

stacked_input <- rbind(stacked_numattrib, stacked_totalnum) %>%
  mutate(paf = paf*100) %>%
  mutate(paf_lci = paf_lci*100) %>%
  mutate(paf_uci = paf_uci*100) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  mutate(paf_ci = paste0(paf,'(',paf_lci,'-',paf_uci,')')) %>%
  mutate(paf_ci = ifelse(group == 'All deaths', '', paf_ci))

## Shared axis plot --------------------------------------------------------

# reorder causes of death
stacked_input <- stacked_input %>%
  mutate(rr_outcome = as.factor(rr_outcome)) %>%
  mutate(rr_outcome = fct_relevel(rr_outcome,
                     "Suicide", "Accidents", "XX: EXTERNAL CAUSES", 
                     "Liver disease", "XI: DIGESTIVE",
                     "Influenza & pneumonia", "COPD", "X: RESPIRATORY",
                     "Cerebrovascular", "Ischaemic heart", "IX: CIRCULATORY",
                     "VI: NERVOUS SYSTEM",
                     "Respiratory", "Lymphoid & haematopoietic", "Female genital", "Digestive", "Breast", "II: CANCERS",
                     "Viral hepatitis", "HIV", "I: INFECTIONS",
                     "Non-drug-related", "Drug-related",
                     "All-cause (opioid users)", "All-cause (inclusion health)"))

g.mid<-ggplot(stacked_input,aes(y=rr_outcome, x=1)) +
  geom_text(aes(label=str_wrap(rr_outcome, width = 30)), size = 5)+
  ylab("NULL")+
  theme(axis.title.y=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        panel.margin = unit(c(0, 0, 0, 0), "mm"),       
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(4,0,10,0), "mm")) + ## cc edit margins
  labs(x = "") # add this to get the spacing the same as the other 2 plots
g.mid

## cc add
stacked_input <- stacked_input %>% 
  mutate(sub = ifelse(group=="Deaths attributable to inclusion health group(s)",
                      numattrib, NA)) %>% 
  arrange(rr_outcome) %>% 
  tidyr::fill(sub) %>% 
  mutate(numattrib = ifelse(numattrib==sub,numattrib,numattrib-sub))
#fin

g1 <- ggplot(data = stacked_input, aes(fill=group, x = rr_outcome, y = numattrib)) +
  geom_bar(position="stack", stat = "identity") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-10,1,1), "mm"),
        legend.position = c(0.3, 0.2),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.text=element_text(size=15),
        text=element_text(size=20)) +
  guides(fill=guide_legend(title = "")) +
  labs(y = "Number of deaths (thousands)") +
  scale_y_reverse() + coord_flip()
g1

# add axis break

g_break <- g1 

g2 <- g_break +
  scale_y_break(c(250, 625), ticklabels=c(600,650,700)) +
  theme(legend.position="none",
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())

g2
## extract legend from original plot
leg = cowplot::get_legend(g_break)

## redraw the figure
p3 <- ggplotify::as.ggplot(print(g2))

## place the legend
g1 <- p3 + ggimage::geom_subview(x=.4, y=.4, subview=leg)

# filter to get PAFs only for the IH groups 
stacked_input_bothsexes_paf <- stacked_input %>% 
  filter(group == "Deaths attributable to inclusion health group(s)")

g2 <- ggplot(data = stacked_input_bothsexes_paf, aes(x = rr_outcome, y = paf)) +xlab(NULL)+
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = paf_lci, ymax = paf_uci), width = .2, col = "black") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(4,4,4,-15), "mm"), ## cc edit margins
        text=element_text(size=20)) +
  labs(y = "PAF (%)") +
  coord_flip()

gg1 <- ggplot_gtable(ggplot_build(g1+theme(plot.margin = unit(c(-1,0,-1,1), "mm"))))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid+theme(plot.margin = unit(c(4,0,10,-17), "mm"))))

grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(3.75/9,1.5/9,3.75/9))
g <- arrangeGrob(gg1,gg.mid,gg2,ncol=3,widths=c(3.75/9,1.5/9,3.75/9)) 
ggsave(file="shared_axis_plot.jpg", g, width = 20, height = 10, dpi = 450)

