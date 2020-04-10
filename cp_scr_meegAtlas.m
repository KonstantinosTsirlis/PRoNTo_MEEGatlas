%% SCRIPT to test/use the function prt_create_meegAtlas

% get the basics sorted
pth_root = 'C:\Dox\1_Code\PRoNTo_Hackathon2020';

% define the input parameters for 1chan/2D case
fn_data = fullfile(pth_root,'condition_Famous.nii');
Vin = spm_vol(fn_data);
blocks_def = struct( ...
    'chans', 1, ...
    'twind', [0 50 ; 51 250 ; 251 800], ...
    'fband', [6 9 ; 12 20]);
labels_def.chans = {'chanXYZ'};
labels_def.twind = cellstr(char('early','medium','late'));
labels_def.fband =  cellstr(char('beta','alpha')) ; 
dorder = 'ftc'; % order of axis in current data: frequency, time, channels
fn_out = fullfile(pth_root,'test_MEEGatlas.nii');
% call to prt_create_meegAtlas
[Vout,fn_labels] = prt_create_meegAtlas(Vin,fn_out,blocks_def,labels_def,dorder);


% Test image creation
fn_meeg = 'rtf_EEG_mapMcbdspmeeg_run_01_sss.mat';
S.D = spm_eeg_load(fullfile(pth_root,fn_meeg));
S.mode = 'channel x time x frequency';
S.channels = 'all'
[images, outroot] = spm_eeg_convert2images(S)
S.mode = 'time x frequency';
[images, outroot] = spm_eeg_convert2images(S)


% define the input parameters for N_chan/3D case
fn_data = fullfile(pth_root,'rtf_EEG_mapMcbdspmeeg_run_01_sss', ...
    'condition_Famous.nii');
Vin = spm_vol(fn_data);
blocks_def = struct( ...
    'chans', 1, ... % Each channel individually
    'twind', [0 50 ; 51 250 ; 251 500], ...
    'fband', [6 9 ; 12 20]);
labels_def = struct( ...
    'chans', {{''}}, ... % Name based on channel index will be created
	'twind', {{'early';'medium';'late'}}, ...
    'fband', {{'beta';'alpha'}} ) ; 
dorder = 'cft'; % order of axis in current data: channels, frequency, time 
fn_out = fullfile(pth_root,'test3D_MEEGatlas.nii');
% call to prt_create_meegAtlas
[Vout,fn_labels] = prt_create_meegAtlas(Vin,fn_out,blocks_def,labels_def,dorder);

blocks_def = struct( ...
    'chans', {{1:35 , 36:70}}, ... % 2 groups of channels (1->35 & 36->70)
    'twind', [0 50 ; 51 250 ; 251 500], ...
    'fband', [6 9 ; 12 20]);
labels_def = struct( ...
    'chans', {{'ChSet1' ; 'ChSet2'}}, ...
	'twind', {{'early' ; 'medium' ; 'late'}}, ...
    'fband', {{'beta' ; 'alpha'}} ) ; 
% call to prt_create_meegAtlas
[Vout,fn_labels] = prt_create_meegAtlas(Vin,fn_out,blocks_def,labels_def,dorder);


% Test voxel list creation
M_dim = [5 7 8];
l_1 = 2:4;
l_2 = 3:6;
l_3 = 1:2;




