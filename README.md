# neuRosim
Simulate FMRI Data

---
title: "Codes for Generating fMRI Data - neuRosim package"
author: "Veronica"
date: "21-06-2016"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## To generate BOLD response

neuRosim uses a stimulus function

Example: to generate a stimulus function for a 20-second ON/OFF block design of 200s with a microtime resolution of 0.1s.

```{r, echo=T}
library(neuRosim)
setwd("C:\\Users\\Vera\\Desktop\\FMRI neuroimaging data")

library(pander)
panderOptions("digits", 2)

totaltime <- 200 #Total time of the design in seconds.
onsets <- seq(1, 200, 40) # onsets of the stimulus in seconds.
dur <- 20 # durations of the stimulus in seconds.
# accuracy is Microtime resolution in seconds
s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 0.1)
str(s)

plot(s, type="l", axes=F, xlab="", main="Block design")
axis(1,c(onsets*10,2010))
axis(2,seq(0,1,0.2))
box()
```

##Three different response functions 

To simulate the BOLD signal caused by the task, the stimulus function is convoluted with a haemodynamic response function. The BOLD signals based on the three convolution functions:

#Gamma HRF

The stimulus function is convoluted with a gamma-variate HRF as implemented in the function gammaHRF with a user-defined full width at half maximum (FWHM) value.

To modulate the strength of the activation in each condition, the argument effectsize should be specified.
                       
```{r}
gamma <- specifydesign(totaltime = 200, 
                       onsets = list(onsets),
                       durations = list(dur), 
                       effectsize = 1, TR = 2, conv = "gamma")

plot(gamma,type="l", axes=F,xlab="", main="gamma-variate HRF")
axis(1,seq(0,100,10))
axis(2,seq(0,1,0.2))
box()
```

#Double gamma HRF

The stimulus function is convoluted with a double-gamma HRF via canonicalHRF,
which models an initial dip and an undershoot of the BOLD signal.

```{r}
canonical <- specifydesign(totaltime = 200, 
                           onsets = list(onsets), 
                           durations = list(dur), 
                           effectsize = 1, TR = 2, conv = "double-gamma")
plot(canonical,type="l", axes=F,xlab="", main="double-gamma")
axis(1,seq(0,100,10))
axis(2,seq(-0.4,1,0.2))
box()
```


#Balloon model

Generates the BOLD signal based on the Balloon model of Buxton et al. (2004)

The solving of the differential equations in the model is based on the Runge-Kutta 
solver in the R package deSolve.

```{r}
balloon <- specifydesign( onsets = list(onsets), 
                          durations = list(dur), 
                          totaltime = 200, 
                          accuracy=0.1, effectsize = 1, TR = 2, conv = "Balloon")

#identical results : balloon(stim=s, totaltime=200, acc=0.1)

plot(balloon,type="l", axes=F, xlab="", main="Balloon model")
axis(1,seq(0,100,10))
axis(2,seq(0,1,0.2))
box()
```

##Spatial location of the activation

A region can be modeled in three ways, as a cube, as a sphere or manually.

#Sphere
For example, to define an activated sphere:

```{r}
a <- specifyregion(dim = c(62, 62), #Dimensions of the image space
                   coord = c(20, 20), #center of the region 
                   radius = 10, #distance in voxel from the center
                   form = "sphere", 
                   fading = 0.3) #0 no fading, 1 results in the fastest decay
#output = matrix 62x62
image(a, axes = T, main= "A region modeled as a sphere", col = grey(seq(1, 0, length = 256)))
box()
```

The fading option can be used to require that the region has the largest effect in the center and smaller activation towards the edges. This fading of the BOLD response is modeled as an exponential decay depending on the distance of the activated voxel to the center of the region.
the fading rate can vary between 0 and 1 with 0 meaning no decay and 1 indicating the strongest
decay.

#Activated voxels defined manually

```{r}
coord <- matrix(c(rep(20, 20), rep(26:30, each = 2), 
                  20:27, 20:27, 
                  rep(28, 6), 21:40, 30:21, 
                  rep(31, 8), rep(40, 8), 33:38), 
                ncol = 2, byrow = FALSE)
#output = matrix 52x2
#columns corresponding to the coordinates
b <- specifyregion(dim = c(64, 64), coord = coord, form = "manual")
image(b, axes = T, main= "Activated voxels defined manually", col = grey(seq(1, 0, length = 256)))
box()
```

#Noise

The noise present in fMRI data is caused by different sources.
The noise functions can be divided into four categories, namely (1) white noise, (2) colored noise, (3) temporal noise and (4) spatial noise.

##White noise

The white noise represents the system noise that is part of the fMRI data. 
Two types of system noise are considered: (1) system noise that is Rician distributed and (2) system noise that is Gaussian distributed.

```{r}
n.white <- systemnoise(dim = 1, 
                       nscan = 100,
                       sigma = 15, #user-defined or based on the desidered SNR
                       type = "rician") 
n.white2 <- systemnoise(dim = 1, 
                       nscan = 100,
                       sigma = 15, #user-defined or based on the desidered SNR
                       type = "gaussian")
plot(seq(1,100,length=100),n.white,type="l", ylim=c(-50,70), col="lightsalmon", main="White noise")
points(seq(1,100,length=100),n.white2,type="l", col="lightblue")
legend("bottomleft",legend = c("Rician", "Gaussian"), col = c("lightsalmon", "lightblue"), lty = 1, lwd = 1, box.col="white") 
box()
```

##Colored noise

Colored noise depends on either the signal, the timing or the location. 
neuRosim contains three types of signal-dependent noise, (1) low-frequency drift, (2) physiological noise and (3) task-related noise.

###Low-frequency drift

Low-frequency drift is a consequence of system noise that can be attributed to slow 
fluctuations in the scanner hardware. The drift is modeled as a basis of discrete cosine functions. The number of functions is determined by the frequency of the drift with a default value of 128 seconds.

Low-frequency drift by freq and Time Repetition

```{r, fig.width=7, fig.height=7}
par(mfrow=c(2,2))  #(1) low-frequency drift - Fluctuations in the scanner hardware   

n.low <- lowfreqdrift(dim = 1, nscan = 100, TR = 2, freq = 10) 
plot(seq(1,100,length=100),n.low,type="l", xlab="", main= "Low-frequency drift \n freq=10 TR=2")

n.low <- lowfreqdrift(dim = 1, nscan = 100, TR = 4, freq = 10) 
plot(seq(1,100,length=100),n.low,type="l", xlab="", main= "Low-frequency drift \n freq=10 TR=4")

n.low <- lowfreqdrift(dim = 1, nscan = 100, TR = 2, freq = 120) 
plot(seq(1,100,length=100),n.low,type="l", xlab="", main= "Low-frequency drift \n freq=120 TR=2")

n.low <- lowfreqdrift(dim = 1, nscan = 100, TR = 4, freq = 120) 
plot(seq(1,100,length=100),n.low,type="l", xlab="", main= "Low-frequency drift \n freq=120 TR=4")
```

###Physiological noise

Physiological noise is defined as possible cardiac and respiratory artefacts and as such accounts for the variability in the signal that is caused by the heart beat and respiratory rate. 

We choose to model the physiological noise separately because it is shown that the frequency of these artefacts is often higher than the scanner fluctuations. 

The physiological noise is modeled as sine and cosine functions with user-defined
frequencies. 

Physiological noise by sigma and Time Repetition

```{r, fig.width=7, fig.height=7}
par(mfrow=c(2,2)) #(2) physiological noise - variability caused by heart beat and respiratory rate

n.phys <- physnoise(dim = 1, nscan = 100, sigma = 15, TR = 2) 
plot(seq(1,100,length=100),n.phys,type="l", xlab="", ylim=c(-60,60), main= "Physiological noise \n sigma=15 TR=2")

n.phys <- physnoise(dim = 1, nscan = 100, sigma = 15, TR = 4) 
plot(seq(1,100,length=100),n.phys,type="l", xlab="", ylim=c(-60,60), main= "Physiological noise \n sigma=15 TR=4")

n.phys <- physnoise(dim = 1, nscan = 100, sigma = 30, TR = 2) 
plot(seq(1,100,length=100),n.phys,type="l", xlab="", ylim=c(-60,60), main= "Physiological noise \n sigma=30 TR=2")

n.phys <- physnoise(dim = 1, nscan = 100, sigma = 30, TR = 4) 
plot(seq(1,100,length=100),n.phys,type="l", xlab="", ylim=c(-60,60), main= "Physiological noise \n sigma=30 TR=4")
```


###Task-related noise

Task-related noise accounts for spontaneous neural activity due to the experimental
task and is operationalized by adding random noise only where and
when activation is present. 

The distribution of this noise can be either Gaussian or Rician. 
The task-related noise can be interpreted as residual noise from head motion that is not removed in the pre-processing stage. 

```{r, fig.width=7, fig.height=7}
par(mfrow=c(2,2))         #(3) task-related noise - residual noise from head motion

s <- stimfunction(totaltime = totaltime, onsets = onsets, durations = dur, accuracy = 0.1)
n.task <- tasknoise(act.image = s, sigma = 15,  type="rician" )
plot(n.task,type="l", ylim=c(-100,150), xlab="", main= "Rician \n sigma=15")
n.task <- tasknoise(act.image = s, sigma = 30, type="rician" ) 
plot(n.task,type="l", ylim=c(-100,150),  xlab="", main= "Rician \n sigma=30")

n.task <- tasknoise(act.image = s, sigma = 15,  type="gaussian" )
plot(n.task,type="l", ylim=c(-100,150),  xlab="", main= "Gaussian \n sigma=15")
n.task <- tasknoise(act.image = s, sigma = 30, type="gaussian" ) 
plot(n.task,type="l", ylim=c(-100,150),  xlab="", main= "Gaussian \n sigma=30")
```

#Temporal noise

Temporal noise accounts for the fact that fMRI data are repeated measurements. 
The function temporalnoise generates noise based on an auto-regressive model of order p (AR(p)).

```{r}
n.temp <- temporalnoise(dim = 1, sigma = 15, nscan = 100, rho = c(0.4, -0.2)) 
# The length of the vector determines the order of the autoregressive model.
# error normal distributed 
plot(seq(1,100,length=100),n.temp,type="l", xlab="", ylab= "Temporal noise", ylim=c(-80,80), 
     main="AR(2) sd=15")
```

#Spatial noise 

Spatial noise models the spatial dependencies in fMRI data.

The function spatialnoise incorporates three types of spatial noise models, namely (1) an autoregressive correlation structure, (2) a Gaussian random field and (3) a Gamma random field. 

The first structure correlates the voxels with each other based on random Gaussian or Rician noise. 
The strength of the correlation depends on the value of the auto-correlation coefficient (default value is rho = 0.75) and the distance between the voxels. 

For example, to generate spatially correlated noise for a 20x20 slice

```{r}
d <- c(20, 20)
n.corr <- spatialnoise(dim = d, sigma = 15, nscan = 100, method = "corr", rho = 0.7)
n.gaus <- spatialnoise(dim = d, sigma = 15, nscan = 100, method = "gaussRF", FWHM = 4)
n.gamma <- spatialnoise(dim = d, sigma = 15, nscan = 100, method = "gammaRF", FWHM = 4, gamma.shape = 3, gamma.rate = 2)

library(plot3D)
library(misc3d)
library(rgl)
colo <-  grey(seq(1, 0.5, length = 256))
#persp3d(n.corr[,,1], col=colo)  #per nscan==1
#persp3d(n.gaus[,,1], col=colo)  #per nscan==1
#persp3d(n.gamma[,,1], col=colo) #per nscan==1
```

#Simulates an fMRI time series

The simTSfmri function generates fMRI time series for a specified design matrix and with
an additive noise structure.
As an example, we will generate a time series for a block design with two conditions. 
The experiment lasts 100 scans with TR = 2 and the first condition has activation blocks of 20s, while the second condition had activation blocks of 7 seconds.

```{r, fig.width=6, fig.height=6}
TR <- 2
nscan <- 100
total <- TR  * nscan
os1 <- seq(1, total, 40)
os2 <- seq(15, total, 40)
dur <- list(20, 7) # block of 20s and block of 7s
os <- list(os1, os2)
effect <- list(3, 10)
w <- c(0.3, 0.3, 0.01, 0.09, 0.3)
# Rician noise (white noise) weight 0.30
# temporal noise weight 0.30
# low-frequency drift 0.01
# Physiological noise 0.09
# task-related noise 0.30
design <- simprepTemporal(totaltime = total, #duration experiment in seconds (200)
                          onsets = os, #onset of each conditions
                          durations = dur, #duration in each condition
                          effectsize = effect, #size1= 3 size2=10
                          TR = TR, #repetion time 2s
                          hrf = "double-gamma") #form of HRF

ts <- simTSfmri(design = design, base = 0, 
                SNR = 2, #variance of the noise
                noise = "mixture", #onr of noise specified before or a mixture
                type = "rician", 
                weights = w, #weights of different noises
                verbose = FALSE)
ts1 <- specifydesign(totaltime = total, 
                     onsets = os1, durations = 20, effectsize = 3, TR = TR, 
                     conv = "double-gamma")

##condition2 without noise

ts2 <- specifydesign(totaltime = total, 
                     onsets = os2, durations = 7, effectsize = 10, TR = TR,
                     conv = "double-gamma")

plot(seq(1,100,length=100),ts,type="l", lty=1,  xlab="", 
     main= "Generated time series \n based on an experiment with two conditions", 
     ylab="Simulated fMRI signal", col=1)
points(seq(1,100,length=100),ts1, type="l", lty=2, col=2) #red first condition
points(seq(1,100,length=100),ts2, type="l", lty=2, col=4) #blue second condition
```

The figure presents the resulting activation from this design in dashed lines. 
It is possible to change the form of the HRF (either gamma, double-gamma or balloon). 
The different noise components are weighted with a vector of weights specified by the user.

#Generating fMRI volumes

For the function simVOLfmri, not only a design matrix, defining when the activation occurs, has to be specified, but also a region, defining where the activation takes place, should be provided.

Suppose that we wish to simulate 2 activated regions that are part of a small network.
Further, we will generate the activation in both regions following the same design matrix as for the generation of the time series.

We will add a mixture of noise with the additional possibility that we can add spatially correlated noise.

```{r}
r1 <- specifyregion(dim = c(64, 64, 64), 
                    coord = c(10, 5, 24), 
                    radius = 10, 
                    form = "sphere")
r2 <- specifyregion(dim = c(64, 64, 64), 
                    coord = c(53, 29, 24), 
                    radius = 5, 
                    form = "sphere")

#image3d(r1, axes = T, col = grey(seq(1, 0, length = 256)))
#image3d(r2, axes = T, col = grey(seq(1, 0, length = 256)))
```

Similarly as for the design matrix, there is a preparation function (simprepSpatial) to specify the regions.

```{r}
regions <- simprepSpatial(regions = 2, 
                          coord = list(c(10, 5, 24), c(53, 29, 24)), 
                          radius = c(10, 5), 
                          form = "sphere")
onset <- list(os, os)
duration <- list(dur, dur)
effect1 <- list(2, 9)
effect2 <- list(14, 8)

design2 <- simprepTemporal(regions = 2, 
                           onsets = onset,
                           durations = duration,
                           TR = TR, 
                           hrf = "double-gamma", #form of HRF
                           effectsize = list(effect1, effect2), 
                           totaltime = total)
w <- c(0.3, 0.3, 0.01, 0.09, 0.1, 0.2)

# data <- simVOLfmri(dim = c(64, 64, 64), 
#                    base = 100, 
#                    design = design2, 
#                    image = regions, # regions = simrepSpatial
#                    SNR = 10, 
#                    noise = "mixture", 
#                    type = "rician", 
#                    weights = w, 
#                    verbose = FALSE)

```

##Generation fMRI volumes without noise

```{r}
data2 <- simVOLfmri(dim = c(64, 64, 64), 
                    base = 100, 
                    design = design2, #generated by simprepTemporal specifying the design
                    image = regions, #specifying the activated regions 
                    noise = "none", 
                    verbose = FALSE)

colo = grey(seq(1, 0, length = 256))
#image3d(data2[,,,1], col=colo)
#image3d(data2[,,,25], col = colo)
#image3d(data2[,,,50], col = colo)
#image3d(data2[,,,75], col = colo)
#image3d(data2[,,,100], col = colo)
#baseline2 <- apply(data2, 1:3, mean) #array 3d 64 x 64 x 64 
#image3d(baseline2, col=colo) #mean given the time
```

#Henson example

##Simulating and analyzing a 4D fMRI dataset

Consider the data from a repetition priming experiment performed using event-related fMRI.

The data are the result of a 2x2 factorial study with factors fame and repetition where famous and non-famous faces were presented twice against a checkerboard.

First we define some parameters like the dimension of the image space, the number of scans and TR. Then, since we simulate an event-related design, we also assign the onsets for each condition.

```{r}
dim <- c(53, 63, 46)
nscan <- 351
TR <- 2
total.time <- nscan * TR
onsets.N1 <- c(6.75, 15.75, 18, 27, 29.25, 31.5, 36, 42.75, 65.25,
               + 74.25, 92.25, 112.5, 119.25, 123.75, 126, 137.25, 141.75,
               + 144, 146.25, 155.25, 159.75, 162, 164.25, 204.75, 238.5) * TR
onsets.N2 <- c(13.5, 40.5, 47.25, 56.25, 90, 94.5, 96.75, 135,
               + 148.5, 184.5, 191.25, 202.5, 216, 234, 236.25, 256.5, 261,
               + 281.25, 290.25, 303.75, 310.5, 319.5, 339.75, 342) * TR
onsets.F1 <- c(0, 2.25, 9, 11.25, 22.5, 45, 51.75, 60.75, 63,
               + 76.5, 78.75, 85.5, 99, 101.25, 103.5, 117, 130.5, 150.75,
               + 171, 189, 227.25, 265.5, 283.5, 285.75, 288, 344.25) * TR
onsets.F2 <- c(33.75, 49.5, 105.75, 153, 157.5, 168.75, 177.75,
               + 180, 182.25, 198, 222.75, 240.75, 254.25, 267.75, 270, 274.4,
               + 294.75, 299.25, 301.5, 315, 317.25, 326.25, 333, 335.25,
               + 337.5, 346.5)
```

## Prepare spatial parameters

We will consider five regions.
The first three are general regions that activate when faces are presented, the fourth region is only activated if famous faces are shown, while in the last region adaptation to the repetition of faces is modeled.

```{r}
region.1A.center <- c(13, 13, 11)
region.1A.radius <- 4
region.1B.center <- c(40, 18, 9)
region.1B.radius <- 6
region.1C.center <- c(10, 45, 24)
region.1C.radius <- 3
region.2.center <- c(15, 16, 31)
region.2.radius <- 5
region.3.center <- c(12, 16, 13)
region.3.radius <- 5
```

##Formatting parameters

In each region, the same design matrix will be considered. However, the effect size in each
condition will vary over conditions.

```{r}
onsets <- list(onsets.N1, onsets.N2, onsets.F1, onsets.F2) #trasform list
onsets.regions <- list(onsets, onsets, onsets, onsets, onsets)
dur <- list(0, 0, 0, 0)
dur.regions <- list(dur, dur, dur, dur, dur)
region.1a.d <- list(160.46, 140.19, 200.16, 160.69)
region.1b.d <- list(140.51, 120.71, 160.55, 120.44)
region.1c.d <- list(120.53, 120.74, 140.02, 100.48)
region.2.d <- list(-0.24, 10.29, 80.18, 160.24)
region.3.d <- list(200.81, 50.04, 240.6, 50.83)
effect <- list(region.1a.d, region.1b.d, region.1c.d, region.2.d,region.3.d)
```

##Create baseline data

We will consider a baseline image. The baseline value for each voxel is determined
as the mean value of the measured time series of that voxel. Non-brain voxels are defined as
voxels with an average measured value less than 250 and are fixed to 0 in the baseline image.
The anatomical structure of the brain will be incorporated in the simulated data.

```{r}
library(oro.nifti)
Hensondata <- readNIfTI("preprocessed_face.nii.gz") #nifti object (4D) 56 x 63 x 46 x 351
str(Hensondata@.Data)

baseline <- apply(Hensondata@.Data, 1:3, mean) #array 53 x 64 x 46
#calculate average using data of nifti object

baseline.bin <- ifelse(baseline > 250, 1, 0) #if average > 250 --> 1
# I have an array 53x63x46 of 0/1

ix <- which(baseline == 1) #positions which average > 250
baseline[-ix] <- 0 #I assign 0 in the positions with average < 250
```

##Preparation functions

We can use the functions simprepTemporal and simprepSpatial to prepare the temporal and spatial structure of our simulated 4D fMRI data.

```{r}
design <- simprepTemporal(regions = 5, onsets = onsets.regions, 
                          durations = dur.regions, 
                          hrf = "double-gamma", 
                          TR = TR,
                          totaltime = total.time,
                          effectsize = effect)

spatial <- simprepSpatial(regions=5, 
                          coord=list(region.1A.center,region.1B.center,
                                     region.1C.center,region.2.center,region.3.center), 
                          radius=c(region.1A.radius,region.1B.radius,region.1C.radius,
                                   region.2.radius,region.3.radius),
                          form="sphere",
                          fading=0.01)

```

##Generate 4D fmri data

We can generate the dataset:
  
  ```{r}
# sim.data <- simVOLfmri(design=design, image=spatial,base=baseline, SNR=3.87,noise="mixture",type="rician", rho.temp=c(0.142, 0.108, 0.084),rho.spat=0.4, w=c(0.05,0.1,0.01,0.09,0.05,0.7), dim=c(53,63,46),nscan=351, vee=0, template=baseline.bin,spat="gaussRF")
sim.data<-readNIfTI("simulated data")
```

We plot the orthographic view of fMRI data for an event-related repetition priming study. Before the data measured by Henson et al. (2002) and after the data simulated by neuRosim.

```{r, fig.width=6, fig.height=6}
orthographic(Hensondata)
orthographic(sim.data)
```


We visualize the real and simulated data.

```{r, fig.width=20, fig.height=12}
image(Hensondata)
sim.data<-readNIfTI("simulated data")
image(sim.data)
```
