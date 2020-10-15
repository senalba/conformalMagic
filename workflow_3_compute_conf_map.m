workflow_0_dir_location

if isLazy
    fileList = dir(strcat(dirBoundary,'*_lazy.mat'));
else
    fileList = dir(strcat(dirBoundary,'*.mat'));
end

L = length(fileList);
fprintf('There are %d files\n',L)

for indK = setdiff(276:L,[3,26,235,275])
    
    tic;
    tStartK = toc;

    if ~isLazy && ~sum(fileList(indK).name(end-7:end-4) == 'lazy')
        
    %load(strcat(dirCropped,name_str(1:end-4),'.mat'),'centerXY')
    name_str = fileList(indK).name(1:end-4);

    fprintf('For map #%3d it takes ',indK)
    
    load(strcat(dirBoundary,fileList(indK).name),'bndrRefined')
    load(strcat(dirCropped,name_str,'.mat'),'centerXY')
    
    %if isLazy
    %    SegOut = imread(strcat(dirCropped,name_str,'_rotated_lazy.tif'));
    %else
    SegOut = imread(strcat(dirCropped,name_str,'_rotated.tif'));
    %end
    
    pC = polygon([bndrRefined]);
    
    tStart = toc;
    f = diskmap(pC);
    f = center(f,0);
    fi = inv(f);
    
    
    save(strcat(dirMap,name_str,'.mat'),'f','fi','name_str')
    
    tStop = toc;
    fprintf('%3d sec.\t %s\n',round(tStop-tStart),name_str)
    
    end
    
end