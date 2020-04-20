# PRoNTo M/EEG atlas creation for MKL
Creating NIfTI atlas for multi-kernel learning in M/EEG data for [PRoNTo](http://www.mlnl.cs.ucl.ac.uk/pronto/).

## Introduction

"Multi-kernel learning" (MKL) needs an atlas to split the data in different (non-overlapping) regions. 

For Brain images, this consists in one of the standard  anatomical atlases, for example the AAL among many others, where each voxel is tagged with an index. All voxels with the same index belong to the same parcel of the brain. When applying the MKL approach, one kernel is built for each such brain parcel and the list of kernels are then linearly mixed by the algorithm. MKL includes a sparsity constraints such that only as few kernels, i.e. brain parcels, are retained in the trained model.

With M/EEG data, there is *stricto sensu* no image but rather time series of signal acquired from a set of channels, distributes on or around the head. And therefore there is no "atlas" as for brain images.

## Issues

With M/EEG, there are thus 2 issues if one wants the MKL machinery developed for brain images: 1/ turning the M/EEG signal into images, and 2/ building a meaningful atlas.

1. SPM provides tools to process M/EEG data and turn them into images for further statistical analysis.  The `spm_eeg_convert2images.m` function does the job for a series of formats

  ````
%   mode       - type of images to generate one of:
%                'scalp x time'
%                'scalp x frequency' (average over time)
%                'scalp' (average over time and frequency)
%                'source' (average over time and frequency)
%                'time x frequency' (average over channels)
%                'time' (1D average over channels, frequency)
%                'frequency' (1D average over channels, time)
%                'average' (average over all dimensions to get a single
%                           number)
  ````
  but not all possibilities are available. Scalp projection (or interpolation) is convenient for scalp-by-time analysis, i.e. without frequency component, or scalp-by-frequency analysis, i.e. after averaging over time.  There is no option though if one is interested in channel-by-time-by-frequency analysis.

2. If M/EEG data are written out as 3D images with dimensions `channels x time x frequency` then one should be able to specify an "atlas" that would be parcellating these 3D volumes into meaningful sets of channels, e.g. depending on their broad location over the scalp, time windows, e.g. early, medium, and late response, and/or or frequency bands, e.g. alpha, beta, and sigma bands.

## Solution

This repository contains a few function to address the aforementioned issues!

1. Two SPM functions, [`spm_eeg_convert2images.m`](spm_eeg_convert2images.m) and [`spm_cfg_eeg_convert2images.m`](spm_cfg_eeg_convert2images.m) were modified to allow the generation of 3D  `channels x time x frequency`  images. The `_cfg_`  function provides the new option in the batch-GUI, while the main function includes the required option:

   ````
   %                'channel x time x frequency' (no averaging)
   ````

   These 2 functions should be copy-pasted in the main SPM folder, overwriting the original function.
   Make sure you do that every time you update to a new SPM version, since they will be replaced by the original ones.

2. The new PRoNTo function [`prt_create_meegAtlas.m`](prt_create_meegAtlas.m)  can create the required "M/EEG atlas" by specifying as input the sets of channels, time windows, and frequency bands to consider. See help for further details on how to use it

  ````
% FORMAT
% [Vout,fn_labels] = prt_create_meegAtlas(Vin,fn_out, blocks,labels,dorder)
% 
% INPUT
% Vin       : header information (see 'spm_vol') of one typical data image
% fn_out    : base name for the generated files, atlas & labels
% blocks    : structure defining the way data are blocked per channels,
%             time-window and frequency-band.
%   .chans  : [Ncg x 1] cell array of channels grouping, each cell contains
%             the list of channels to group together.
%             Shortcuts, use 0 for all channels together, i.e. one group,
%             and 1 for each channel separately, i.e. Nch groups.
%   .twind  : [Ntw x 2] array with time-window boundaries
%   .fband  : [Nfb x 2] array with frequency-band boundaries
% labels    : structure with the corresponding labels for channels,
%             time-windows, and frequency-bands. If not provided, then
%             labels will be built based on the index of channel group,
%             time-window, freq-band.
%   .chans  \
%   .twind  | cell array with labels.
%   .fband  /
% dorder    : ordering of the data as a string made of 'c' -> channels,
%             'f' -> frequency, 't' -> time. Default 'cft'.
  ````

Finally the small script [`cp_scr_meegAtlas.m`](cp_scr_meegAtlas.m) just demonstrates how the function can be used to generate the M/EEG atlas on example data.
