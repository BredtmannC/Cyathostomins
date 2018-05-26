setwd("I:/WE13/Parasitologie/Christina Bredtmann/Method/2018.02.14-MALDI-R/Cylicostephanus longibursatus F/LON_F_01")
getwd()
library(MALDIquant)
library(MALDIquantForeign)
library(readMzXmlData)

LON_F_01_Directory <- ("I:/WE13/Parasitologie/Christina Bredtmann/Method/2018.02.14-MALDI-R/Cylicostephanus longibursatus F/LON_F_01")
LON_F_01<-importMzXml(LON_F_01_Directory,verbose = FALSE,centroided=FALSE)
names(LON_F_01) <- paste0("replicate",1:9)

table(sapply(LON_F_01, length))
typeof(LON_F_01)
class(LON_F_01)
summary(LON_F_01)

spec1 <- readMzXmlFile(file.path("I:/WE13/Parasitologie/Christina Bredtmann/Method/2018.02.14-MALDI-R/Cylicostephanus longibursatus F/LON_F_01/0_A1/1/1SLin/LON_f01_A1_1.mzXML"))
print(spec1$metaData)
plot(spec1$spectrum$mass, spec1$spectrum$intensity, type="l")



length(LON_F_01)
LON_F_01[1:9]

any(sapply(LON_F_01, isEmpty))

table(sapply(LON_F_01, length))

#not regular, don´t know why. Don´t know how to change it.
all(sapply(LON_F_01, isRegular))

#help by Sebastian Gibbs: "non regular" kann ignoriert werden, 
#wenn diese Plots dem Plot mit fiedler2009subset ähnlich sehen.
#Ja, tun sie!
mzDifflonf01 <- diff(mass(LON_F_01[[1]]))
plot(mzDifflonf01)
hist(mzDifflonf01)


plot(LON_F_01[[4]])
plot(LON_F_01[[2]])

spectra_LON_f01 <- transformIntensity(LON_F_01, method=("sqrt"))
spectra_LON_f01 <- smoothIntensity(spectra_LON_f01, method="SavitzkyGolay",
                           halfWindowSize=10)
baseline <- estimateBaseline(spectra_LON_f01[[6]], method="SNIP",
                             iterations=100)
plot(spectra_LON_f01[[6]])
lines(baseline, col="red", lwd=2)

spectra_LON_f01 <- removeBaseline(spectra_LON_f01, method="SNIP",
                          iterations=100)
plot(spectra_LON_f01[[1]])
plot(spectra_LON_f01[[4]])
plot(spectra_LON_f01[[7]])


spectra_LON_f01 <- calibrateIntensity(spectra_LON_f01, method="TIC")
spectra_LON_f01 <- alignSpectra(spectra_LON_f01,
                        halfWindowSize=20,
                        SNR=2,
                        tolerance=0.002,
                        warpingMethod="lowess")



avgLON_f01 <- averageMassSpectra(spectra_LON_f01,
                                 method="mean")

noise <- estimateNoise(avgLON_f01)
plot(avgLON_f01, xlim=c(4000, 5000), ylim=c(0, 0.002))
lines(noise, col="red")
lines(noise[,1], noise[, 2]*2, col="blue")
lines(noise[,1], noise[, 2]*3, col="green")

plot(avgLON_f01, xlim=c(3000, 6000), ylim=c(0, 0.002))
lines(noise, col="red")
lines(noise[,1], noise[, 2]*2, col="blue")
lines(noise[,1], noise[, 2]*3, col="green")

plot(avgLON_f01, xlim=c(2000, 20000), ylim=c(0, 0.002))
lines(noise, col="red")
lines(noise[,1], noise[, 2]*2, col="blue")
lines(noise[,1], noise[, 2]*3, col="green")

peaks <- detectPeaks(avgLON_f01, method="MAD",
                     halfWindowSize=20, SNR=2)
plot(avgLON_f01, xlim=c(4000, 5000), ylim=c(0, 0.002))
points(peaks, col="red", pch=4)

plot(avgLON_f01, xlim=c(2000, 20000), ylim=c(0, 0.002))
points(peaks, col="red", pch=4)

str(peaks)
list(peaks)
class(peaks)
length(peaks)

isMassPeaks(peaks)

#not working because only one list available!
# peaks<-binPeaks(peaks) 
# peaks<-filterPeaks(peaks, minFrequency = 0.25)
