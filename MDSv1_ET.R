# Script to create multi-dimensional scaling network 


rm(list=ls())  # Clear workspace 
graphics.off() # Close all graphics 

setwd("/Users/Eldin/Desktop/My_working_Bix/R")

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

mydata <- read.csv("output.csv")

d <- dist(mydata) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim

fit # view results
family <- as.factor(d) 

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric	MDS",	type="n", pch=19, col=family)
text(x, y, labels = colnames(mydata), cex=.7)

