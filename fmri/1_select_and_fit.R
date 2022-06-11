# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("fmri","fda.usc","parallel")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})
source('misc_functions.R')

# read in the image data
ImageFileName = "bold.nii";
BOLD_sub007_run001 = read.NIFTI(ImageFileName); ##### reads in the fMRI image
BoldData = extractData(BOLD_sub007_run001); #### this is the 4-dim array
Mask = BOLD_sub007_run001$mask; #### Mask = TRUE are the place where there is actually 
								#### the active parts of a head

# stimulus outputs
Face.Stimulus.Onsets = read.table(file= "cond001.txt", header = FALSE);
Durations.Face = rep(0, nrow(Face.Stimulus.Onsets));
HRF.Face = fmri.stimulus(scans = BOLD_sub007_run001$dim0[4], onsets = 0.5*Face.Stimulus.Onsets[,1],
durations = Durations.Face);

Scramble.Stimulus.Onsets = read.table(file= "cond002.txt", header = FALSE);
Durations.Scramble = rep(0, nrow(Scramble.Stimulus.Onsets));
HRF.Scramble = fmri.stimulus(scans = BOLD_sub007_run001$dim0[4], onsets = 0.5*Scramble.Stimulus.Onsets[,1],
durations = Durations.Scramble);

BoldX = fmri.design(cbind(HRF.Face, HRF.Scramble), order = 2);
dim(BoldX)

# flattening out the data for future use.
# comment out this part is this is not needed
Dim = dim(BoldData)
SpaceDim = Dim[1]*Dim[2]*Dim[3];
TimeDim = Dim[4];
YData = matrix(nrow = TimeDim, ncol = SpaceDim);
DataMap = matrix(nrow = 4, ncol = SpaceDim);


for (i in 1: Dim[1]){
	for (j in 1: Dim[2]){
		for (k in 1: Dim[3]){
			a = (i - 1)*Dim[2]*Dim[3] + (j -1)*Dim[3] + k;
			DataMap[ , a] = c(Mask[i , j, k], i, j, k);
			YData[ , a] = BoldData [i, j, k, ];
		}
	}
}

## do selection over all voxels
all.grid = expand.grid(2:(Dim[1]-1), 2:(Dim[2]-1), 2:(Dim[3]-1))
jfun = function(j) lm.select.spatial(X=BoldX, coords=as.numeric(all.grid[j,]),
                                     sd.vec = seq(.05, .5, by=.02), threshold=.8)
all.list = list()
chunks = seq(1,1001,by=100)-1
nc = 1
for(nchunk in 1:(length(chunks)-1)){
  all.list = mclapply((chunks[nchunk]+1):chunks[nchunk+1], jfun, mc.cores=16)
  cat(paste("Iteration",nc,"done\n"))
  save(all.list,file=paste0('all_list_',nc,'.Rda'))
  nc = nc+1
}
