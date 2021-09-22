library(dbplyr)
library(dplyr)
library(corrplot)
library(tidyverse)
library(tidyselect)
library(tidyr)
#import dataset and set to dataframe variable
dataset =read.csv("corr2u.csv")
r <- dataset %>% select(Lectin,T.cell.Function,Coefficient)
r1 = data.frame(r)
#sets an index to variables, this removes duplicates
df = r1 %>% 
  group_by(T.cell.Function) %>% 
  mutate(grouped_id = row_number())

#After data has been group, variable data is spread against dependent variables
df1 = df %>% 
  spread(T.cell.Function, Coefficient) %>% 
  select(-grouped_id)
#data shaping and set to write csv
df2 = data.matrix(df1)
rownames(df2) = df1$Lectin
df2 = df2[,2:16]
#write csv
write.csv(df2,file = "corr2u2.csv", na = "NA")
# AFTER ROW UNITE FROM DUPLICATES, dataframe is prepped and shaped for corrplot
dataset1 =read.csv("corr2u2.csv")
data = data.matrix(dataset1)
rownames(data) = dataset1$X
colnames(data) = colnames(df1, do.NULL = TRUE)
data = data[,2:16]
res1 <- cor.mtest(data, conf.level = 0.95)
col3 <- colorRampPalette(c("dark blue", "white", "dark red")) 
corrplot(data, method="circle", is.corr=FALSE, tl.col = "Black", tl.srt = 45, col = col3(100))
