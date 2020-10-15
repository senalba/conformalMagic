workflow_0_dir_location


%   fileList = dir(strcat(dirCropped,'*.tif'));
%fileList = dir(strcat(dirSource,'*.tif'));
fileListLayers = dir(strcat(dirMask,'*.tif'));
%fileListVein = dir(strcat(dirSource,'*_veins.tif'));
L = length(fileListLayers);
fprintf('There are %d files\n',L)

if isMutant
    
    name_str = fileListLayers(1).name(1:end-11);
elseif isHawaiian
    
    name_str = fileListLayers(1).name(1:end-9);
end

f = uint8(imread(strcat(dirSource,name_str,'.tif')));
nX = size(f,1);
nY = size(f,2);
fWing = uint8(zeros(nX,nY));
fVein = uint8(zeros(nX,nY));
fModel = uint8(zeros(nX,nY));


seDisk1 = strel('disk',1);
seDisk2 = strel('disk',2);
seDisk3 = strel('disk',3);
seDisk5 = strel('disk',5);
seDisk7 = strel('disk',7);
seDisk8 = strel('disk',8);
seDisk10 = strel('disk',10);
seDisk12 = strel('disk',12);
seDisk20 = strel('disk',20);
seDisk40 = strel('disk',40);
seDisk50 = strel('disk',50);
seDisk100 = strel('disk',100);
seDiamond = strel('diamond',1);
se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);

tic;
for indK = 1:1%[1:15,480:495,1000:1015,1500:1515,2000:2015]
    
    tStart = toc;
    
    if isMutant
        name_str = fileListLayers(indK).name(1:end-11);
    elseif isHawaiian
        name_str = fileListLayers(indK).name(1:end - 9);
    end
    %name_str  = 'egfr_M_R_lei_4X_94';
    
    fprintf('%3d/%3d || %s',indK,L, name_str)
    
    f = im2uint8(imread(strcat(dirSource,name_str,'.tif')));
    
    if isMutant
        fLayers = uint8(imread(strcat(dirSource,name_str,'_layers.tif')));
    elseif isHawaiian
        fLayers = uint8(imread(strcat(dirMask,name_str,'_mask.tif')));
    end
    
    
    indWing = find(fLayers == 255); % wing's surface
    %  indVein = find(fLayers == 170); % vein's
    indBackground = find(fLayers == 85); % backgorund
    
    if isMutant
        indBackground = find(fLayers == 85); % backgorund
    elseif isHawaiian
        indBackground = find(fLayers == 0); % backgorund
    end
    
    fWing(:) = uint8(0);
    
    if ~isHawaiian
        fVein(:) = uint8(0);
        
        fVein(indVein) = uint8(255);
        
        fVein = imerode(fVein,seDisk3);
        fVein = imfill(imdilate(fVein,seDisk8));
        fVein = imerode(fVein,seDisk12);
        
        fVein = uint8(255.*bwareaopen(fVein,50000));
        
        indBackgroundClean = find(fVein ~= 255);
        
        indWingExtra = intersect(indBackgroundClean,indBackground);
        
        indBackground = indBackgroundClean;
        
        indWing = union(indWingExtra,indWing);
        
        indWing = setdiff(indWing,indBackgroundClean);
        fWing(indWing) = uint8(255);
        
        
        fWing = uint8(255.*bwareaopen(imfill(fWing),1000));
        
        
        fWingOpen = imclose(fWing,seDisk40);
        
        fWingDilate = imdilate(fWingOpen,seDisk1);
        
        fVeinErode = imerode(imfill(fVein),seDisk2);
        
        fFinal = imopen(fWingDilate  + fVeinErode ,seDisk5);
    else
        fLayers = imfill(fLayers);
        
        fFinal = imerode(fLayers,seDisk5);
        fFinal = imclose(fLayers,seDisk10);
        bndr = bwperim(fLayers);
    end
    
    
    
    [L1, N] = bwlabel(fFinal);
    D = regionprops(L1, 'Area');
    
    if N > 1
        area_values = [D.Area];
        %        idx = find(max(area_values));
        % we remove all small objects
        fFinal = bwareaopen(fFinal,max(area_values));
        %fprintf('More than one region left. Extra images deleted.\n')
    end
    
    
    indF = find(fFinal == 0);
    
    bndr = bwperim(fFinal);
    indBndr = find(bndr == 1);
    
    f1 = f(:,:,1);
    f2 = f(:,:,2);
    %f3 = f(:,:,3);
    
    f1(indF) = 0;
    f2(indF) = 0;
    %f3(indF) = 0;
    
    f(:,:,1) = f1;
    f(:,:,2) = f2;
    %f(:,:,3) = f3;
    
    fWingsmaller = imerode(fWing,seDisk2);
    
    indWing = find(fWingsmaller == 255);
    indVein = setdiff(setdiff(1:nX*nY,indBackgroundClean),indWing);
    
    fModel(:) = 0;
    fModel(indWing) = 255;%85;
    fModel(indVein) = 170;%255;
    fModel(indBackground) = 85;%170;
    fModel(indBndr) = 0;
    %size(indBackground)
    %size(indBackgroundClean)
    
    %imshowpair(f,fModel,'montage');
    
    save(strcat(dirMask,name_str,'.mat'),...
        'bndr','indF','indWing','indVein','indBackground','indBndr',...
        'name_str')
    imwrite(fModel,strcat(dirMask,name_str,'.tif'))
    tStop = toc;
    
    fprintf('\t|| %1.2f sec\n',tStop-tStart)
end
