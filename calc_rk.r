library(spatstat)
library(maptools)

args = commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]
border_file <- args[3]
print(input_file)
print(output_file)

d = read.csv(input_file,sep=',')

epsg <- "+init=epsg:32736"
border <- readShapePoly(border_file, proj4string = CRS(epsg))

# assumes that the first row is the header
myppp <- ppp(d[2:dim(d)[1],1],d[2:dim(d)[1],2],as.owin(border))
k <- Kest(myppp)
write.csv(k,file=output_file,row.names=FALSE)
