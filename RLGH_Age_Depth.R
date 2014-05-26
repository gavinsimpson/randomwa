setwd("D:\\Text\\Teaching\\Bariloche\\Age-Depth\\Software\\WinBacon_2.2")

source("Bacon.R")

x <- read.csv("cores\\RLGH3\\RLGH3.csv")

Bacon("RLGH3", depths.file=TRUE)
y
