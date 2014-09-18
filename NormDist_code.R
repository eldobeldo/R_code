#Determine whether a vector is normally distributed (e.g. 0.001, 0.01, 0.1, 0.0001, 0.02, etc)

#Input file is minorpopulation.csv; output histogram with normal distribution curve 

rm(list=ls())           #Clear workspace == ensure no leftover variables
graphics.off();         #Close all graphics windows

getwd()                 #Get current working directory

setwd("/Users/eldin/Desktop/R working dir/")  #Replace with your desired working directory

search()                #See current packages loaded


data = as.numeric(unlist(read.csv("minorpopulation.csv",stringsAsFactors=FALSE)))

hist(data)

#Log transform data and run mean, standard deviation and create histogram with normal dist curve

hist(log10(data),ylab="Frequency",xlab="Log of frequency of minor population",main="Distribution of frequencies of minor populations")
mean(log10(data))
sd(log10(data))
x = seq(from=-15, to =0, by = 0.01)
y = dnorm(x,mean(log10(data)),sd(log10(data)))
points(x,y*max(hist(log10(data),plot=FALSE)$counts)/max(y),typ="l",lwd=2,col="blue")
1-pnorm(log10(0.002),mean(log10(data)),sd(log10(data)))

hist((data))
mean((data))
sd((data))
x = seq(from=-15, to =0, by = 0.01)
y = dnorm(x,mean((data)),sd((data)))
points(x,y*max(hist((data))$counts)/max(y),typ="l",lwd=2,col="blue")
1-pnorm((0.001),mean((data)),sd((data)))
