# fMRI data analysis

Following are the description and functionality of each file.

- `bold.nii` loads the 4-dimensional array containing fMRI images across a number of time points for our test subjects.
- `cond00*.txt` contain our test subject's responses to the experimental stimulii.
- `1_select_and_fit.R` contains code for processing the datasets and selecting the best subset of predictors for each regression model.
- `2_process_outputs.R` processes outputs from the model selection procedure into voxelwise F-statistics and p-values.
- `3_plot_outputs.R` ingests the p-values and produces plots (corresponding to Fig. E.1 in the main paper)
- `test_Select_scaled_spatial_pvalue1.nii` contain a 3D array with voxelwise p-values.

The smoothed image in the main paper (Fig. 7.1c) is produced by visualizing a cross-section of the p-values with [FMRIB Software Library](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki).