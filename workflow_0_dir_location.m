addpath linecurvature_version1b/
addpath sc_map/

N = 100;
resolutioN = 200;

isHawaiian = 0;
isMutant = 0;
isDiet = 0;
isInbred = 1;

%%
if isHawaiian
    
    dirLocation = ...
        '/Volumes/vasyl_alba/datasets/Edwards-Hawaiian Drosophila wings-complete set/';
    dirSource  = strcat(dirLocation,'img_01_source/');
    %dirSourceRotated  = strcat(dirLocation,'img_05_source_rotated/');
    
elseif isMutant
    
    %there are 5 mutants egfr, mam, samw, star, tkv
    
    %gray-scale images
    %dirLocation = strcat('/Volumes/vasyl_alba/datasets/data_mutants/',...
    %    'Microscope_1/40X_magnification_Olympus/');
    
    %color images 
    dirLocation = strcat('/Volumes/vasyl_alba/datasets/data_mutants/',...
        'Microscope_2/40X_magnification_Leica/');
    
    dirSource  = strcat(dirLocation,'img_00_originals/');
    %dirSourceRotated  = strcat(dirLocation,'img_05_originals_rotated/');
elseif isDiet
    
    dirLocation = '/Volumes/vasyl_alba/datasets/jamie_temp_diet/';
    dirSource = dirLocation;
    %dirSourceRotated  = strcat(dirLocation,'img_05_originals_rotated/');
elseif isInbred
     
    dirLocation = strcat('/Volumes/vasyl_alba/datasets/inbred_lines',...
         '/');
    %dirLocation = '/Volumes/home/wings_backup/'
    
    dirSource  = strcat(dirLocation,'img_00_originals/');
end





dirLayer = strcat(dirLocation,'img_01_layers/');
dirMask = strcat(dirLocation,'img_02_mask/');
dirMaskCropped = strcat(dirLocation,'img_03_mask_cropped/');
dirBoundary = strcat(dirLocation,'img_04_boundary/');
dirSourceRotated  = strcat(dirLocation,'img_05_originals_rotated/');
%dirCropped = strcat(dirLocation,'img_cropped/');
dirMap = strcat(dirLocation,'img_06_map/');

dirImage = sprintf('%simg_07_mapped/%d/',dirLocation,resolutioN);

dirAngularAligned = sprintf('%simg_08_angular_aligned/%d/',...
    dirLocation,resolutioN);

dirInitiallyAligned = sprintf('%simg_09_initially_aligned/%d/',...
    dirLocation,resolutioN);

dirVeinNormal = sprintf('%simg_10_masks_veins_normalized/',...
    dirLocation);

dirGauge = sprintf('%simg_11_gauge_transformation/%d/',...
    dirLocation,resolutioN);

dirVeinMapped = sprintf('%simg_14_veins_mapped/%d/',...
    dirLocation,resolutioN);

dirVeinAligned = sprintf('%simg_12_veins_aligned/%d/',...
    dirLocation,resolutioN);

dirPCA = sprintf('%simg_13_pca_wings/',dirLocation);



%dirAligned = sprintf('%simg_aligned/%d/',dirLocation,resolutioN);

isLazy = 0;

badItems = [2,13,18,22,24];
tic;