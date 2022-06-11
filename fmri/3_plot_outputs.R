# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("oro.nifti")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})
source('misc_functions.R')

# read in the image data
BoldData = extractData(read.NIFTI("bold.nii")) # fMRI image
Boldp1 = extractData(read.NIFTI("test_Select_scaled_spatial_pvalue1.nii")) # p-values

## plots
plotfunz = function(i){
  image(BoldData[,,i,1], col=grey.colors(12))
  image(Boldp1[,,i,1],add=T,col=colorRampPalette(c("red","black"))(8))
}

for(i in 1:33){
  plotfunz(i)
}
# z = 12

## plots
plotfunx = function(i){
  i = 65-i
  image(BoldData[i,,,1], col=grey.colors(12))
  image(Boldp1[i,,,1],add=T,col=colorRampPalette(c("red","black"))(8))
}

for(i in 1:64){
  plotfunx(i)
}
# x = 48

## plots
plotfuny = function(i){
  image(BoldData[,i,,1], col=grey.colors(12))
  image(Boldp1[,i,,1],add=T,col=colorRampPalette(c("red","black"))(8))
}

for(i in 1:64){
  plotfuny(i)
}
# y = 7

pdf('xslice.pdf',4.5,4); plotfunx(48); dev.off()
pdf('yslice.pdf',4.5,4); plotfuny(7); dev.off()
pdf('zslice.pdf',4,4.5); plotfunz(12); dev.off()
