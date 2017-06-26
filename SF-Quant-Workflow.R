library(tidyverse)
setwd("~/Microscopy/20170310 mlck and rock2 on off rfpla pmlc ppmlc")

#LOAD DATA
SF <- read.delim(file="SF Combined Analysis.txt", sep="\t", stringsAsFactors=T)

#Average (Int/Area)
SF.summary <- SF %>%                                                     
  filter(Channel != "SF") %>%
  mutate(Channel = as.factor(ifelse(Channel == "p-mlc" | 
                            Channel == "p-MLC" |
                            Channel == "pMLC", "pMLC",
                          ifelse(Channel == "pp-mlc" |
                                   Channel == "ppMLC", "ppMLC", "oops"))),
         Compartment = as.factor(ifelse(Compartment == "c" 
                                        | Compartment == "C", "C",
                              ifelse(Compartment == "p" |
                                       Compartment == "P", "P", "oops"))),
         Kinase = as.factor(ifelse(Kinase == "ROCK" |
                                     Kinase == "ROCK2", "ROCK2",
                                   ifelse(Kinase == "MLCK", "MLCK",
                                          ifelse(Kinase == "Empty", "Empty", "Kinase Error"))))) %>%
  mutate(Dox=as.factor(Dox),
         Channel = as.factor(Channel)) %>%                            # format Dox column
  group_by(Date,Kinase, Dox, Compartment, Channel, Cell, Line) %>%     # group by conditional factors & per cell
  mutate(NormInt=IntDen/Area) %>%                       # calculate avg int/area per condition
  summarize(NormInt = mean(NormInt)) # calculate avg int/area per cell

SF.means.0 <- SF.summary %>%                                #calculate mean at each 0 dox condition 
  filter(Dox == 0) %>%                                      #for Kinase/channel/compartment groups
  group_by(Date,Kinase,Compartment, Channel, Line) %>%
  summarize(Mean0Int = mean(NormInt))

SF.norm <- SF.summary %>% 
  left_join(SF.means.0, by = c('Date', 'Kinase', 'Compartment', 'Channel', 'Line')) %>%
  mutate(NewNorm = NormInt/Mean0Int) %>%
  dplyr::select(Kinase, Dox, Compartment, Channel, NewNorm, Line)

SF.stats <- SF.norm %>%
  group_by(Kinase,Dox, Compartment,Channel, Line) %>%
  summarise(IntList=list(NewNorm)) %>%
  spread(key = Dox, value = IntList, sep = '-') %>%
  group_by(Kinase,Compartment,Channel, Line) %>%
  #mutate(SW0 = shapiro.test(unlist(`Dox-0`))$p.value, SW200 = shapiro.test(unlist(`Dox-200`))$p.value) %>%  #run if you want to check if normal
  #mutate(p_value = wilcox.test(unlist(`Dox-0`),unlist(`Dox-800`),var.equal=FALSE)$p.value) %>%  #run if U251 data, not U2OS data
  mutate(p_value = ifelse(is.null(unlist(`Dox-200`)), wilcox.test(unlist(`Dox-0`),unlist(`Dox-800`),var.equal=FALSE)$p.value,wilcox.test(unlist(`Dox-0`),unlist(`Dox-200`),var.equal=FALSE)$p.value)) %>%
  mutate(sig = ifelse(p_value < 0.001, '***',ifelse(p_value <0.01, '**',ifelse(p_value <0.05, '*',''))))

SF.asterisk.MLCK  <- SF.stats %>%
  filter(Kinase == "MLCK") %>%
  #mutate(xaxis = 2, yaxis = 1.1*(ifelse(is.null(max(unlist(`Dox-200`)))) %>% run for U251, not U2OS
  mutate(xaxis = 2, yaxis = 1.08*(ifelse(is.null(unlist(`Dox-200`)), max(unlist(`Dox-800`)), max(unlist(`Dox-200`))))) %>%
  dplyr::select(xaxis,yaxis, sig, Kinase, Compartment, Channel, Line) %>%
  mutate(CompartmentLabel = 
           ifelse(Compartment == "C", 
                  "Center", 
                  ifelse(Compartment == "P",
                         "Periphery",
                         "Oops")))
SF.asterisk.ROCK  <- SF.stats %>%
  filter(Kinase == "ROCK2") %>%
  #mutate(xaxis = 2, yaxis = 1.1*(ifelse(is.null(max(unlist(`Dox-200`)))) %>% run for U251, not U2OS
  mutate(xaxis = 2, yaxis = 1.01*(ifelse(is.null(unlist(`Dox-200`)), max(unlist(`Dox-800`)), max(unlist(`Dox-200`))))) %>%
  dplyr::select(xaxis,yaxis, sig, Kinase, Compartment, Channel, Line) %>%
  mutate(CompartmentLabel = 
           ifelse(Compartment == "C", 
                  "Center", 
                  ifelse(Compartment == "P",
                         "Periphery",
                         "Oops")))



#Graphing MLCK:  
SF.norm %>% 
  filter(Kinase == "MLCK") %>%
  mutate(CompartmentLabel = 
           ifelse(Compartment == "C", 
                  "Center", 
                  ifelse(Compartment == "P",
                         "Periphery",
                         "Oops"))) %>%
  ggplot() +                                                #plot
  aes(x=Dox, y=NewNorm) +
  stat_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, 
              aes(color = Compartment),
              size = 0.5,
              pch = 16,
              alpha = 1) +
  facet_grid(Channel ~ Line + CompartmentLabel, scales="free") + 
  scale_color_manual(values = c("C" = "gray42", "P" = "black"),
                     labels=c("Central SFs", "Peripheral SFs")) +
  scale_y_continuous(limits=c(0, NA)) +
  ylab("Normalized Fluorescence Intensity (AU)") + 
  xlab("Dox (ng/ml)") + 
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        panel.border = element_rect(fill = NA),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank()) +
  geom_text(data = SF.asterisk.MLCK, aes(x=xaxis, y=yaxis, label = sig), color = 'navyblue', size = 6)


# Graphing ROCK
SF.norm %>% 
  filter(Kinase == "ROCK2") %>%
  mutate(CompartmentLabel = 
           ifelse(Compartment == "C", 
                  "Center", 
                  ifelse(Compartment == "P",
                         "Periphery",
                         "Oops"))) %>%
  ggplot() +                                                #plot
  aes(x=Dox, y=NewNorm) +
  stat_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, 
              aes(color = Compartment),
              size = 1,
              pch = 16,
              alpha = 1) +
  facet_grid(Channel ~ Line + CompartmentLabel, scales="free") + 
  scale_color_manual(values = c("C" = "gray42", "P" = "black"),
                     labels=c("Central SFs", "Peripheral SFs")) +
  scale_y_continuous(limits=c(0, NA)) +
  ylab("Normalized Fluorescence Intensity (AU)") + 
  xlab("Dox (ng/ml)") + 
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        panel.border = element_rect(fill = NA),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank()) +
  geom_text(data = SF.asterisk.ROCK, aes(x=xaxis, y=yaxis, label = sig), color = 'navyblue', size = 6)


# Graphing U373
SF.asterisk.U373  <- SF.stats %>%
  filter(Line == "U373", Kinase != "Empty") %>%
  mutate(xaxis = 2, yaxis = 1.01*(ifelse(is.null(unlist(`Dox-200`)), max(unlist(`Dox-800`)), max(unlist(`Dox-200`))))) %>%
  dplyr::select(xaxis,yaxis, sig, Kinase, Compartment, Channel, Line) %>%
  mutate(CompartmentLabel = 
           ifelse(Compartment == "C", 
                  "Center", 
                  ifelse(Compartment == "P",
                         "Periphery",
                         "Oops"))) 
SF.norm %>% 
  filter(Line == "U373", Kinase != "Empty") %>%
  mutate(CompartmentLabel = 
           ifelse(Compartment == "C", 
                  "Center", 
                  ifelse(Compartment == "P",
                         "Periphery",
                         "Oops"))) %>%
  ggplot() +                                                #plot
  aes(x=Dox, y=NewNorm) +
  stat_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, 
              aes(color = Compartment),
              size = 0.5,
              pch = 16,
              alpha = 1) +
  facet_grid(Channel ~ Kinase + CompartmentLabel, scales="free") + 
  scale_color_manual(values = c("C" = "gray42", "P" = "black"),
                     labels=c("Central SFs", "Peripheral SFs")) +
  scale_y_continuous(limits=c(0, NA)) +
  ylab("Normalized Fluorescence Intensity (AU)") + 
  xlab("Dox (ng/ml)") + 
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        panel.border = element_rect(fill = NA),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 20),
        panel.grid = element_blank()) +
  geom_text(data = SF.asterisk.U373, aes(x=xaxis, y=yaxis, label = sig), color = 'midnightblue', size = 4)


#Add number of conditions with:
#   stat_summary(fun.data=function(x) c(y=0, label=length(x)), geom='text') +



#Plot To Look Like Elena's
######
SF.norm %>% 
  filter(Kinase == "MLCK", Channel == "pMLC") %>%
  ggplot() +                                                #plot
  aes(x=Dox, y=NewNorm) +
  stat_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, 
              aes(fill = Compartment),
              color = "black",
              pch = 21,
              size = 1) +
  facet_grid(. ~ Compartment) + 
  scale_fill_manual(values = c("C" = "white", "P" = "black"),
                    labels=c("Central SFs", "Peripheral SFs")) +
  scale_y_continuous(limits=c(0, NA)) +
  ylab("p-MLC Intensity (AU)") + 
  xlab("Dox (ng/ml)") + 
  labs(title="MLCK - pMLC") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(), 
        axis.line = element_line(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.justification = c(0,1),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),
        strip.text.x = element_blank()) 
SF.norm %>% 
  filter(Kinase == "MLCK", Channel == "ppMLC") %>%
  ggplot() +                                                #plot
  aes(x=Dox, y=NewNorm) +
  stat_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, 
              aes(fill = Compartment),
              color = "black",
              pch = 21,
              size = 1) +
  facet_grid(. ~ Compartment) + 
  scale_fill_manual(values = c("C" = "white", "P" = "black"),
                    labels=c("Central SFs", "Peripheral SFs")) +
  scale_y_continuous(limits=c(0, NA)) +
  ylab("pp-MLC Intensity (AU)") + 
  xlab("Dox (ng/ml)") + 
  labs(title="MLCK - ppMLC") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(), 
        axis.line = element_line(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.justification = c(0,1),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),
        strip.text.x = element_blank()) 
SF.norm %>% 
  filter(Kinase == "ROCK2", Channel == "pMLC") %>%
  ggplot() +                                                #plot
  aes(x=Dox, y=NewNorm) +
  stat_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, 
              aes(fill = Compartment),
              color = "black",
              pch = 21,
              size = 1) +
  facet_grid(. ~ Compartment) + 
  scale_fill_manual(values = c("c" = "white", "p" = "black"),
                    labels=c("Central SFs", "Peripheral SFs")) +
  scale_y_continuous(limits=c(0, NA)) +
  ylab("p-MLC Intensity (AU)") + 
  xlab("Dox (ng/ml)") + 
  labs(title="ROCK2 - pMLC") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(), 
        axis.line = element_line(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.justification = c(0,1),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),
        strip.text.x = element_blank()) 
SF.norm %>% 
  filter(Kinase == "ROCK2", Channel == "ppMLC") %>%
  ggplot() +                                                #plot
  aes(x=Dox, y=NewNorm) +
  stat_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, 
              aes(fill = Compartment),
              color = "black",
              pch = 21,
              size = 1) +
  facet_grid(. ~ Compartment) + 
  scale_fill_manual(values = c("c" = "white", "p" = "black"),
                    labels=c("Central SFs", "Peripheral SFs")) +
  scale_y_continuous(limits=c(0, NA)) +
  ylab("pp-MLC Intensity (AU)") + 
  xlab("Dox (ng/ml)") + 
  labs(title="ROCK2 - ppMLC") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(), 
        axis.line = element_line(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.justification = c(0,1),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),
        strip.text.x = element_blank()) 
