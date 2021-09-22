library(reshape)
library(ggplot2)
library(tidyverse)
library(dbplyr)
library(dplyr)

#This is a growth rate calculator. The experiment is to test quality of FBS from 3 vendors.
#Which FBS is the best? 1e5 cells were seeded and counted every 24hr time points. 
#Calculate growth logistic and plot

#enter seed seeding density and K maximum capacity of the environment. 
#Typical cell culture confluency is 1e6cells/mL
seed = 1e5
K = 1e6

r = read.csv("fbs.csv")

J5A8C = r[1:2]
J5A8T = r[1:3]
J5A8T = J5A8T %>% select(-X5A8C)

par(mfrow = c(1:2))

Jurkat = read.csv("JURKAT.csv")
J5A8Cd = mutate(J5A8C, J5A8CdR= log(X5A8C/seed),
                J5A8CGrowthRate = J5A8CdR/Time,
                J5A8CGrowth = J5A8CGrowthRate*X5A8C*(X5A8C)/K)
J5A8Cd = mutate(J5A8Cd, J5A8CGrowth1 = X5A8C/max(X5A8C))
J5A8Cd%>%
  ggplot(aes(Time, J5A8CGrowth1)) +
  geom_point(col="blue") +
  geom_smooth( method = "glm",
               method.args = list(family="binomial"),
               data = J5A8Cd, se= F, col = "Blue")

J5A8Td = mutate(J5A8T, J5A8TdR= log(X5A8T/seed),
                J5A8TGrowthRate = J5A8TdR/Time,
                J5A8TGrowth = J5A8TGrowthRate*X5A8T*(X5A8T)/K)
J5A8Td = mutate(J5A8Td, J5A8TGrowth1 = X5A8T/max(X5A8T))
J5A8Td%>%
  ggplot(aes(Time, J5A8TGrowth1)) +
  geom_point(col="red") +
  geom_smooth( method = "glm",
               method.args = list(family="binomial"),
               data = J5A8Td, se= F, col = "red")

confluency = max(J10.6d$J10.6GrowthRate1)*0.8
J10.6 = read.csv("JLAT10.6.csv")
J10.6d = mutate(J10.6, J10.6dR= log(Count/seed),
               J10.6GrowthRate = J10.6dR/Time,
               J10.6Growth= J10.6GrowthRate*Count*(Count)/K)
J10.6d = mutate(J10.6d, J10.6GrowthRate1 = Count/max(Count))
J10.6d %>%
  ggplot(aes(Time, J10.6GrowthRate1)) +
  geom_point(col="red")+
  geom_smooth( method = "glm",
               method.args = list(family="binomial"(link="probit")),
               data = J10.6d, se= F, col = "red") +
  geom_segment(data=J10.6d,
               aes(x=Time, xend=Time, y=0, yend=J10.6GrowthRate1),
               colour="blue", linetype = "11") +
  geom_segment(data=J10.6d,
               aes(x=Time, xend=Time, y=J10.6GrowthRate1, yend=J10.6GrowthRate1),
               colour="blue", linetype = "11") +
  geom_text(data=J10.6d,
            aes(label=round(Time, 1), x=Time, y=-0.03),
            size=3, colour="blue")
  theme


J6.3 = read.csv("JLAT6.3.csv")
J6.3d = mutate(J6.3, J6.3dR= log(Count/seed),
              J6.3GrowthRate = J6.3dR/Time,
              J6.3Growth= J6.3GrowthRate*Count*(Count)/K)
J6.3d = mutate(J6.3d, J6.3GrowthRate1 = Count/max(Count))
J6.3d %>%
  ggplot(aes(Time, J6.3GrowthRate1)) +
  geom_point(col="blue") +
  geom_smooth( method = "glm",
               method.args = list(family="binomial"),
               data = J6.3d, se= F, col = "blue")
theme


J5A8CG = J5A8Cd[c("Time", "J5A8CGrowth1")]
J5A8TG = J5A8Td[c("Time", "J5A8TGrowth1")]
J6.3G = J6.3d[c("Time", "J6.3GrowthRate1")]
gather = left_join(J5A8CG, J5A8TG, by="Time")
gather2 = left_join(gather, J6.3G, by = "Time")

gather <- melt(gather, id = "Time")
#saveRDS(gather2, file = "gather.Rda") #save elements(e.g. data.frame) as Rfile
#loadRDS(file = "gather.Rda")

ggplot(data = gather, aes(x = Time, y = value, color = variable)) +
  geom_point()+
  geom_smooth(mapping = aes(Time, J5A8CGrowth1),method = "glm",
              method.args = list(family="binomial"),data = J5A8CG,
              se= F, col = "red")+
  geom_smooth(mapping = aes(Time, J5A8TGrowth1),method = "glm",
              method.args = list(family="binomial"), data = J5A8TG,
              se= F, col = "Green")+ ylab("Growth")+
  ggtitle("JLAT 5A8 GIBCO FBS Growth Rates:
          Current Lot:1891244 - Testing Lot:2156016")

  geom_smooth(mapping = aes(Time, J6.3GrowthRate1),method = "glm",
              method.args = list(family="binomial"),
              data = J6.3G, se= F, col = "Blue")
plt
theme

model <- glm(formula= J6.3GrowthRate1 ~ Time, data=J6.3G, family=binomial)
summary(model)
plt  + theme
theme <-  theme(
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot and panel background
  plot.background=element_rect(fill = "black"),
  panel.background = element_rect(fill = 'grey14'),
  # Change legend
  legend.position = c(0.6, 0.07),
  legend.direction = "horizontal",
  legend.background = element_rect(fill = "grey14", color = "grey14"),
  legend.key = element_rect(color = "gray", fill = "black"),
  legend.title = element_text(color = "white"),
  legend.text = element_text(color = "white")
)

Time=seq(0,33, by=0.1)
Seed = 3e6
K = 3e7

#Volume = 10*3

newdata = data.frame(Time = Time)
p = predict(model, newdata, type="response")

df = data.frame(Time)
df = df %>% mutate(Pred = p, PredPop = ((K*Pred)))
rate = df[c("Time", "PredPop")]

points <- tibble(
  age = rate$Time,
  prop = rate$PredPop)

colors <- list(
  data = "#41414550",
  # data = "grey80",
  fit = "#414145")
xs <- seq(0, 48, length.out = 80)

# Create the curve from the equation parameters
trend <- tibble(
  age = rate$Time,
  asymptote = .8,
  scale = .2,
  midpoint = 48,
  prop = asymptote / (1 + exp((midpoint - age) * scale)))

ggplot(points) +
  aes(x = age, y = prop) +
  geom_line(data = trend, color = colors$fit) +
  geom_point(size = 3.5, color = colors$data) +
  scale_x_continuous(
    name = "Age in months",
    limits = c(0, 40),
    breaks = scales::extended_breaks(Q = c(0, 6))) +
  scale_y_continuous(
    name = "Intelligibility",
    limits = c(0, NA),
    labels = scales::percent_format(accuracy = 1))





MaxTime = 20
rngTime <-filter(df,Time < MaxTime)
PopDiff = max(rngTime$PredPop)
PopDiff

MaxPop = 2e6
rng <-filter(df,PredPop < MaxPop)
TimeDiff = max(rng$Time)
TimeDiff

newD = data.frame(Time = 6)
pred = predict(model, newD, type="response")*K
pred

df%>%
  ggplot(aes(Time, PredPop)) +
  geom_point(col="red") +
  geom_smooth( method = "glm",
               method.args = list(family="binomial"),
               data = df, se= F, col = "red")

#ggplot(gather2, aes(Time, y = value, color = CellLine)) +
#geom_point(aes(y = JurkatG, col = "Jurkat")) +
#geom_point(aes(y = J10.6G, col = "JLAT10.6"))+
#geom_point(aes(y = J6.3G, col = "JLAT6.3"))+
#  geom_smooth(mapping = aes(Time, JurkatG),method = "glm", method.args = list(family="binomial"),
#  data = Jurkatd, se= F, col = "blue")+
#geom_smooth(mapping = aes(Time, J10.6G), method = "glm", method.args = list(family="binomial"),
#   data = J10.6d, se= F, col = "red") +
#geom_smooth(mapping = aes(Time, J6.3G), method = "glm", method.args = list(family="binomial"),
#data = J6.3d, se= F, col = "green")


