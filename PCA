library("FactoMineR")
library("factoextra")


data = read.csv("/Users/zeitgeist/Desktop/lectin/TimTissue/TIMCOLON-LECTIN.csv")
dataSample = data[,-c(2)]
rownames(dataSample) = dataSample[,1]
dataSample = dataSample[,-c(1)]

res.pca <- PCA(dataSample, graph = TRUE)
print(res.pca)
eig.val <- get_eigenvalue(res.pca)
eig.val
#~~~~~Scree Plot~~~~~~~~
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

var <- get_pca_var(res.pca)
var
# Coordinates
head(var$coord[,1:2])
# Cos2: quality on the factore map
head(var$cor[,1:2])
# Contributions to the principal components
head(var$contrib)

#~~~~~~BiPlot~~~~~~~~
library("corrplot")
corrplot(t(var$cor[,1:2]), is.corr=FALSE)
corrplot(var$contrib, is.corr=FALSE)
# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)
# Color by cos2 values: quality on the factor map
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$cor[,1:2], centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(res.pca, col.var = grp,
             legend.title = "Cluster")

res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)

# The variable Species (index = 5) is removed
# before PCA analysis
df = res.pca$call$X
category = decathlon2$Competition[1:23]
all = cbind(df,dataSample)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (not "text")
             ellipse.level=0.5,
             col.ind = data$GROUP, # color by groups
             palette = c("#999999", "#FF4D47", "#03D4FB"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")

#PCA Variable Corr heat map
col3 <- colorRampPalette(c("dark blue", "blue","white","red", "dark red"))
corrplot(t(var$cor[,1:2]), is.corr=FALSE, col = col3(100))

var$cor[,1:2]

df = iris[,-5]
df.pca = prcomp(df)

fviz_pca_ind(res.pca,
             # Fill individuals by groups
             geom.ind = "point",
             pointshape = 21,
             pointsize = 2.5,
             fill.ind = data$Group,
             col.ind = "black",
             # Color variable by groups
             #col.var = factor(c("sepal", "sepal", "petal", "petal")),
             
             legend.title = list(fill = "Species", color = "Clusters"),
             repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("jco")      # Variable colors

