# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("fmri","fda.usc","oro.nifti","parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})
source('misc_functions.R')

# read in the image data
ImageFileName = "bold.nii";
BOLD_sub007_run001 = read.NIFTI(ImageFileName); ##### reads in the fMRI image
BoldData = extractData(BOLD_sub007_run001); #### this is the 4-dim array
Mask = BOLD_sub007_run001$mask; #### Mask = TRUE are the place where there is actually 
								#### the active parts of a head
                
## after model selection: spatial model
# set non-selected indices NA for each voxel
coef.mat = matrix(NA, nrow=135168, ncol=30)
for(nchunk in 1:12){
  load(paste0('all_list_',nchunk,'.Rda'))
  for(pos in 1:length(all.list)){
    po = all.list[[pos]]$position
    a = (po[1]-1)*Dim[2]*Dim[3] + (po[2]-1)*Dim[3] + po[3]
    # if(a > 13e4){
    #   break
    # }
    best.a = all.list[[pos]]$best.index
    beta.vec = all.list[[pos]]$beta
    coef.mat[a,] = beta.vec
    coef.mat[a,-best.a] = NA
  }
  cat(paste("Chunk",nchunk,"done!\n"))
}

Dim = dim(BoldData)
BoldBeta = array(0, dim=c(Dim[-4],30))
for (i in 1: Dim[1]){
  for (j in 1: Dim[2]){
    for (k in 1: Dim[3]){
      a = (i - 1)*Dim[2]*Dim[3] + (j -1)*Dim[3] + k;
      BoldBeta[i,j,k,] = coef.mat[a,]
    }
  }
}

niftiObj = nifti(BoldBeta, datatype=16, pixdim=BOLD_sub007_run001$header$pixdim)
writeNIfTI(niftiObj, file="test_Select_scaled_spatial")

# voxelwise count of the number of non-zero effect neighbors
BoldEffect = array(NA, dim=Dim)
BoldF = array(NA, dim=c(Dim[-4],1))
Boldp = array(NA, dim=c(Dim[-4],1))
for (i in 1: Dim[1]){
  for (j in 1: Dim[2]){
    for (k in 1: Dim[3]){
      a = (i - 1)*Dim[2]*Dim[3] + (j -1)*Dim[3] + k;
      
      # get predictor matrix and record effects
      if(i==1 | j==1 | k==1 | i==64 | j==64 | k==33){
        neighbor.y = matrix(0, nrow=nrow(YData), ncol=26)
      } else{
        grid = expand.grid((i-1):(i+1),(j-1):(j+1),(k-1):(k+1))
        aa = as.numeric((grid[,1] - 1)*Dim[2]*Dim[3] + (grid[,2] -1)*Dim[3] + grid[,3])
        neighbor.y = YData[,aa]
      }
      ## control for all-zero responses
      if(min(diag(var(neighbor.y)))>0){
        total.x = cbind(neighbor.y[,-14],BoldX[,-3])
        scale.x = scale(total.x)
        
        # spatial effects across time
        best.index = which(!is.na(coef.mat[a,1:26]))
        p1 = length(best.index)
        if(p1>0){
          BoldEffect[i,j,k,] = matrix(scale.x[,best.index], ncol=p1, byrow=F) %*%
            matrix(coef.mat[a,best.index], ncol=1)
          SSreg1 = sum(BoldEffect[i,j,k,]^2)
          # scale back to original scale
          # BoldEffect[i,j,k,] = BoldEffect[i,j,k,]*sd(neighbor.y[,14]) + mean(neighbor.y[,14])
        }
        
        # F-statistic
        best.index.full = which(!is.na(coef.mat[a,]))
        p0 = length(best.index.full)
        if(p0>0){
          SSreg0 = sum((matrix(scale.x[,best.index.full], ncol=p0, byrow=F) %*%
                         matrix(coef.mat[a,best.index.full], ncol=1))^2)
          SSres = var(neighbor.y[,14]) + SSreg0 - SSreg1
          BoldF[i,j,k,1] = (SSreg1/SSres) * ((Dim[4]-p0+p1)/p1)
          Boldp[i,j,k,1] = pf(BoldF[i,j,k,1],p1,Dim[4]-p0+p1)
        }
      }
    }
  }
}

BoldEffect0 = BoldEffect
for(l in 1:dim(BoldEffect0)[4]){
  tempimg = BoldEffect0[,,,l]
  tempimg[which(Mask==FALSE,arr.ind=T)] = NA
  BoldEffect0[,,,l] = tempimg
}
BoldEffect1 = BoldEffect0
for(i in 1:Dim[1]){
  BoldEffect1[i,,,] = BoldEffect0[65-i,,,]
}

niftiObj = nifti(BoldEffect1, datatype=16, pixdim=BOLD_sub007_run001$header$pixdim)
writeNIfTI(niftiObj, file="test_Select_scaled_spatial_effect")

# # Mask
# tempimg = BoldEffect1[,,,1]
# tempimg[which(Mask==FALSE,arr.ind=T)] = NA
# image(tempimg[,,10])

# frames after face stimulus
Face.frames = ceiling(Face.Stimulus.Onsets[,1]/2)
niftiObj = nifti(BoldEffect1[,,,Face.frames],
                 datatype=16, pixdim=BOLD_sub007_run001$header$pixdim)
writeNIfTI(niftiObj, file="test_Select_scaled_spatial_effect_face")

FaceMean = array(0,c(Dim(BoldEffect)[-4],1))
FaceMean[,,,1] = apply(BoldEffect1[,,,Face.frames],1:3,mean)
niftiObj = nifti(FaceMean, datatype=16, pixdim=BOLD_sub007_run001$header$pixdim)
writeNIfTI(niftiObj, file="test_Select_scaled_spatial_effect_facemean")

# frames after face stimulus
Scramble.frames = ceiling(Scramble.Stimulus.Onsets[,1]/2)
niftiObj = nifti(BoldEffect1[,,,Scramble.frames],
                 datatype=16, pixdim=BOLD_sub007_run001$header$pixdim)
writeNIfTI(niftiObj, file="test_Select_scaled_spatial_effect_scramble")

ScrambleMean = array(0,c(Dim(BoldEffect)[-4],1))
ScrambleMean[,,,1] = apply(BoldEffect1[,,,Scramble.frames],1:3,mean)
niftiObj = nifti(ScrambleMean, datatype=16, pixdim=BOLD_sub007_run001$header$pixdim)
writeNIfTI(niftiObj, file="test_Select_scaled_spatial_effect_scramblemean")


# F statistic
BoldF[which(BoldF<0,arr.ind=T)] = NA
BoldF[which(BoldF>30,arr.ind=T)] = 30
tempimg = BoldF[,,,1]
tempimg[which(Mask==FALSE,arr.ind=T)] = NA
BoldF[,,,1] = tempimg
BoldF1 = BoldF
for(i in 1:Dim[1]){
  BoldF1[i,,,] = BoldF[65-i,,,]
}

niftiObj = nifti(BoldF1, datatype=16, pixdim=BOLD_sub007_run001$header$pixdim)
writeNIfTI(niftiObj, file="test_Select_scaled_spatial_Fstat")

# p-value
tempimg = Boldp[,,,1]
tempimg[which(Mask==FALSE,arr.ind=T)] = NA
Boldp[,,,1] = tempimg
Boldp1 = Boldp
for(i in 1:Dim[1]){
  Boldp1[i,,,] = Boldp[65-i,,,]
}
Boldp1[which(Boldp<.95,arr.ind=T)] = NA

niftiObj = nifti(Boldp1, datatype=16, pixdim=BOLD_sub007_run001$header$pixdim)
writeNIfTI(niftiObj, file="test_Select_scaled_spatial_pvalue1")

