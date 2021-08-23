#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("hrbrthemes")
library(hrbrthemes)

# Load dataset from github
##export CheckSSCMax.txt file
checksscmax <- as.data.frame(checksscmax)
data<-checksscmax

# Plot
## Drug 1
d1_sscmax_ck<-data %>%
  ggplot( aes(x=Time, y=D1)) +
  geom_line( color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=0.5) +
  geom_line(aes(x=Time, y=SSCMax_D1), color="red")+
  #theme_ipsum() +
  ggtitle("Checking Steady State CMax")

  #print drug 1 sscmax graph
  d1_sscmax_ck

## Drug 2
d2_sscmax_ck<-data %>%
  ggplot( aes(x=Time, y=D2)) +
  geom_line( color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=0.5) +
  geom_line(aes(x=Time, y=SSCMax_D2), color="red")+
  #theme_ipsum() +
  ggtitle("Checking Steady State CMax")

  #print drug 2 sscmax graph
  d2_sscmax_ck

## Drug 3
d3_sscmax_ck<-data %>%
  ggplot( aes(x=Time, y=D3)) +
  geom_line( color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=0.5) +
  geom_line(aes(x=Time, y=SSCMax_D3), color="red")+
  #theme_ipsum() +
  ggtitle("Checking Steady State CMax")

  #print drug 3 sscmax graph
  d3_sscmax_ck
  
## Drug 4
d4_sscmax_ck <- data %>%
  ggplot( aes(x=Time, y=D4)) +
  geom_line( color="darkgrey") +
  #geom_point(shape=21, color="black", fill="#69b3a2", size=0.5) +
  geom_line(aes(x=Time, y=SSCMax_D4), color="red")+
  #theme_ipsum() +
  ggtitle("Checking Steady State CMax")

  #print drug 4 sscmax graph
  d4_sscmax_ck