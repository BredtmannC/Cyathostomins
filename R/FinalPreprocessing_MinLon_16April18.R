## Load libraries and data
library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)

# load directory + import data + save data
# MINLON_Directory <- ("../rawData/")
# MINLON <- importMzXml(MINLON_Directory,verbose = FALSE,centroided = FALSE)
# save(file = paste0(MINLON_Directory, "MINLON.RData"),
#      list="MINLON")

load(file = "../data/MINLON.RData")

## Preprocessing
sc.results <- screenSpectra(MINLON) #calculate atypicability score
summary(sc.results) # no atypical spectra! # changes to atypical (#142 after trimming)
plot(sc.results, labels=T) 
rm(sc.results)

# Some plots of raw data to be happy about ourselves

plot(MINLON[[1]])
plot(MINLON[[46]])
plot(MINLON[[84]])
plot(MINLON[[117]])
plot(MINLON[[142]])

# updating MINLON without faulty spectra - not necessary because no faulty spectra!
# MINLON <- sc.results$fspectra

# Create informations for our samples
TinaDF <- data.frame(spectrum_ID =  paste0("spectrum",1:162),
                     sample_ID = NA,
                     species = NA,
                     OTU = NA, #new
                     sex = NA,
                     species_sex = NA,
                     OTU_sex = NA, #new
                     ID_OTU = NA) #new

# Horrible nasty regular expression to extract correct names
changeName <- function(x){
  sub('_([^_]*)$', '', sub('_([^_]*)$', '', sub(".*\\\\([^.]*).*", "\\1", gsub("/", "\\\\",x))))
}

# Create labels for all our samples, in order to average the spectra
labels <- NA
for (i in 1:length(MINLON)){
  labels <- c(labels, changeName(MINLON[[i]]@metaData[["file"]]))
}

TinaDF$sample_ID <- na.omit(labels)

# Create species information
TinaDF$species[grep(pattern = "LON", TinaDF$sample_ID)] <- "Cylicostephanus_longibursatus"
TinaDF$species[grep(pattern = "MIN", TinaDF$sample_ID)] <- "Cylicostephanus_minutus"

# Create sex information
TinaDF$sex[grep(pattern = "_f", TinaDF$sample_ID)] <- "female"
TinaDF$sex[grep(pattern = "_m", TinaDF$sample_ID)] <- "male"

# Create species*sex information
TinaDF$species_sex <- paste0(TinaDF$species, "_", TinaDF$sex)

# Create OTU-information
TinaDF$OTU[grep(pattern = "MIN_f02", TinaDF$sample_ID)] <- "OTU_1"

TinaDF$OTU[grep(pattern = "MIN_f01|MIN_m02|MIN_m09|MIN_m03|MIN_m04", 
                TinaDF$sample_ID)] <- "OTU_3"

TinaDF$OTU[grep(pattern = 
                  "MIN_f04|MIN_f03|MIN_m05|MIN_m06|MIN_m07|MIN_m01|MIN_m10|MIN_f05|MIN_m08", 
                TinaDF$sample_ID)] <- "OTU_2"

TinaDF$OTU[grep(pattern = "LON", TinaDF$sample_ID)] <- "OTU_LON"

# Create OTU*sex information:
TinaDF$OTU_sex <- paste0(TinaDF$OTU, "_", TinaDF$sex)

# Create OTU*ID information:
TinaDF$ID_OTU <- paste0(TinaDF$sample_ID, "_", TinaDF$OTU)

# Create average Tina info HAHAHA
avgTina.info <- TinaDF[!duplicated(TinaDF$sample_ID),]
##################################################

# Name the spectra 
names(MINLON) <- TinaDF$spectrum_ID

rm(i, labels, changeName)

########
# Check data basically
typeof(MINLON)
# class(MINLON) same as typeof
length(MINLON) # n data 
table(sapply(MINLON, length)) # 27 values in range of 20613, 135 values in range of 20657
# summary(MINLON)-long table showing length of each spectrum 

# trim spectra to same length
MINLON <- trim(MINLON,c(2000,20000))

# check if some spectra are empty
any(sapply(MINLON, isEmpty))

# check if datapoints are evenly distributed
all(sapply(MINLON, isRegular)) #ignore, because datapoints are almost evenly distributed

#######
## Prepocessing:  #TODO find out the best method

# http://fmwww.bc.edu/repec/bocode/t/transint.html
# Square root
# The square root, x to x^(1/2) = sqrt(x), is a transformation with a
# moderate effect on distribution shape: it is weaker than the logarithm
# and the cube root. It is also used for reducing right skewness, and also
# has the advantage that it can be applied to zero values. Note that the
# square root of an area has the units of a length. It is commonly applied
# to counted data, especially if the values are mostly rather small.

## To have a look at the different possible transformations:
plot(MINLON[[1]]) # no transfo
plot(transformIntensity(MINLON, method=("sqrt"))[[1]])
plot(transformIntensity(MINLON, method=("log"))[[1]])

# We go for square root transformation
spectra <- transformIntensity(MINLON, method=("sqrt"))  
plot(spectra[[1]])

# estimation of best half window size for smoothing: 
# "Full width at half maximum of peaks"
plot(spectra[[1]], type="b", xlim=c(4650, 4750), ylim=c(80, 130))
abline(h=105, col=4, lty=2)   #count 10

spectra <- smoothIntensity(spectra, method="SavitzkyGolay",
                           halfWindowSize=10)
plot(spectra[[1]])

# estimate baseline and plot 

# Method= SNNIP: Define iteration steps: 25,50, 75, 100
iterations <- seq(from=25, to=100, by=25)
# define different colours for each step
col <- rainbow(length(iterations))
plot(spectra[[1]])
# draw different baseline estimates
for (i in seq(along=iterations)){
  baseline <- estimateBaseline (spectra[[1]], method="SNIP",
                                iterations=iterations[i])
  lines(baseline, col=col[i], lwd=2)}
legend("topright", legend=iterations, col=col,lwd=1)

#Method=TopHat
bTopHat1 <- estimateBaseline(spectra[[1]], method="TopHat", 
                             halfWindowSize=25)
bTopHat2 <- estimateBaseline(spectra[[1]], method="TopHat", 
                             halfWindowSize=75)
bTopHat3 <- estimateBaseline(spectra[[1]], method="TopHat",
                             halfWindowSize=150)
bTopHat4 <- estimateBaseline(spectra[[1]], method="TopHat",
                             halfWindowSize=200)

plot(spectra[[1]],main="TopHat Baseline")
lines(bTopHat1, col=2)
lines(bTopHat2, col=3)
lines(bTopHat3, col=4)
lines(bTopHat4,col=5)

legend("topright", lwd=1, legend = paste0("halfWindowSize=", c(25,75,150,200)),
       col=c(2,3,4,5))

# Choose which baseline correction is right! SNIP or TopHat?
rm(baseline)
rm(bTopHat1, bTopHat2, bTopHat3, bTopHat4)
rm(iterations)
rm(col)

# Substract baseline for all samples and plot to be happy
spectra <- removeBaseline(spectra, method="SNIP",
                          iterations=100)

# or: spectra <- removeBaseline(spectra, method="TopHat")

plot(spectra[[1]])
plot(spectra[[46]])
plot(spectra[[84]])
plot(spectra[[117]])

# Normalisation: put all spectra to the same scale
spectra <- calibrateIntensity(spectra, method="TIC")  ##PQN also possible
plot(spectra[[1]]) # intensity very low, too low?
plot(spectra[[46]])
plot(spectra[[84]])
plot(spectra[[117]])

# Align spectra: Reference: to which the samples should be aligned?
# If missing, reference peaks is used. 
# We don't have reference peaks, so we should give a reference. Maybe sample 1?
spectra <- alignSpectra(reference = spectra[[1]],
                        spectra,
                        halfWindowSize=10,
                        SNR=2,
                        tolerance=0.002,
                        warpingMethod="lowess")

# Average the spectra according to the sample_IDs
avgSpectra <- averageMassSpectra(l = spectra, 
                                 labels = TinaDF$sample_ID,
                                 method = "mean")

# Estimate the noise (SignalNoiseRatio) of the average spectra, for one example
# define SNRs steps: 1, 1.5, ..3
snrs <- seq(from=1, to=3, by=0.5)
## define different colors for each step
colSNR <- rainbow(length(snrs))
## estimate noise
noise <- estimateNoise(avgSpectra[[1]])
plot(avgSpectra[[1]], xlim=c(4000, 6000), ylim=c(0, 0.0015))
for (i in seq(along=snrs)){
  lines(noise[, "mass"],
        noise[, "intensity"]*snrs[i],
        col=colSNR[i], lwd=1)}
legend("topright", legend=snrs, col=colSNR, lwd=1)

# TODO try in different spectra, to check consistency
plot(avgSpectra[[1]], xlim=c(2000, 20000), ylim=c(0, 0.002))
for (i in seq(along=snrs)){
  lines(noise[, "mass"],
        noise[, "intensity"]*snrs[i],
        col=colSNR[i], lwd=1)}
legend("topright", legend=snrs, col=colSNR, lwd=1)

# noise <- estimateNoise(avgSpectra[[1]])
# plot(avgSpectra[[1]], xlim=c(3000, 6000), ylim=c(0, 0.002))
# lines(noise, col="red")
# lines(noise[,1], noise[, 2]*3, col="green")
# lines(noise[,1], noise[, 2]*2, col="blue")

rm(snrs, colSNR, noise, i)

# Assignment of peaks for all average spectra
peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=10, SNR=2)

# Plot one example to be happy
plot(avgSpectra[[1]], xlim=c(4000, 5000), ylim=c(0, 0.002))
points(peaks[[1]], col="red", pch=4)

plot(avgSpectra[[1]], xlim=c(2000, 20000), ylim=c(0, 0.002))
points(peaks[[1]], col="red", pch=4)

save(file = paste0("./","avgSpectra.RData"), list="avgSpectra")
save(file = paste0("./","avgTina.info.RData"), list="avgTina.info")
save(file = paste0("./","peaks.RData"), list="peaks")







