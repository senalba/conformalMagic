workflow_0_dir_location

fileListDone = dir(strcat(dirMaskCropped,'*.mat'));

if isDiet
    
    fileListNew = dir(strcat(dirMask,'G_*.tif'));
else
    
    fileListNew = dir(strcat(dirMask,'*.tif'));
end

fileList = find_new_items(fileListNew,fileListDone);
fileList = fileListNew;

lFileList = length(fileList);
fprintf('There are %d files\n',lFileList)

seDisk1 = strel('disk',1);
seDisk2 = strel('disk',2);


%%
for indK = 1:lFileList%min(10,lFileList)
    
    
    name_str = fileList(indK).name(1:end-4);
    fprintf('%3d/%3d || %s\n',indK,lFileList, name_str)
    
    load(strcat(dirMask,name_str,'.mat'),'bndr')
    
    %f = rgb2gray(im2uint8(imread(strcat(dirSource,name_str,'.tif'))));
    f = im2uint8(imread(strcat(dirSource,name_str,'.tif')));
    
    if isMutant
         
        name_cell = strsplit(name_str,'_');
    
    elseif isHawaiian
         
        name_cell = strsplit(name_str,{' '});
         
    elseif isDiet
        
        name_cell = strsplit(name_str,'_');
    elseif isInbred
        
        name_cell = strsplit(name_str,'_');
    end
    
    
    fMask = imread(strcat(dirMask,name_str,'.tif'));
    [n1, n2] = size(f);
    
    bndrImage = imfill(255.*uint8(bndr));
    
%    imshow(f.*uint8(imcomplement(bndr)))
   
%    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
     if isMutant
         title(sprintf('indK = %d/%d, %s %s %s %s',...
        indK,lFileList, name_cell{1}, name_cell{2},...
        name_cell{3}, name_cell{6}))
    elseif isHawaiian
         title(sprintf('indK = %d/%d, %s %s %s %s',...
        indK,lFileList, name_cell{1}, name_cell{2},...
        name_cell{3}, name_cell{4}))
    elseif isDiet
        
        title(sprintf('indK = %d/%d, %s %s %s %s %s',...
        indK,lFileList, name_cell{2}, name_cell{3},...
        name_cell{4}, name_cell{5}, name_cell{7}))
    elseif isInbred
        
%        title(sprintf('indK = %d/%d, %s %s %s %s',...
%        indK,L, name_cell{2}, name_cell{3},...
%        name_cell{4}, name_cell{5}))
	end
    
    
%    hold on
    
    load(   strcat(dirMaskCropped,name_str,'.mat')  ,...
            'centerXY','hbXY','alXY','cXY'          );

    %hbXY = ginput(1); % humeral break
%    plot(hbXY(1) + 1i* hbXY(2),'*')
%    plot(round(hbXY(1)) + 1i* round(hbXY(2)),'b*')
    
    %alXY = ginput(1); % alula notch
%    plot(alXY(1) + 1i* alXY(2),'r*')
    
    
    %cXY = ginput(1); %extra point for a smooth cut
%    plot(cXY(1) + 1i* cXY(2),'g*')
    
    
    index1 = lineBetweenPoints(n1,n2,hbXY,alXY);
    index2 = lineBetweenPoints(n1,n2,cXY,alXY);
    
    bndrImage(index1) = 0;  % Set the line points to black
    bndrImage(index2) = 0;  % Set the line points to black
    
    bndrImage = imerode(bndrImage,seDisk1);
    
    
    [L1, N] = bwlabel(bndrImage);
    D = regionprops(L1, 'Area');
    
    if N > 1
        area_values = [D.Area];
        idx = find(max(area_values));
        %we remove all small objects
        bndrImage = 255.*uint8(bwareaopen(bndrImage,max(area_values)));
        %fprintf('More than one region left. Extra images deleted.\n')
    end
    
    bndrImage = imdilate(bndrImage,seDisk1);
    
    bndrImage(index1) = 255;  % Set the line points to white
    bndrImage(index2) = 255; % Set the line points to white
    
    bndrImage = imfill(bndrImage);
    
    bndrImage = imopen(bndrImage,seDisk2);
    
    bndrImage = imdilate(bndrImage,seDisk1);
    bndr = bwperim(bndrImage);
    indBndr = find(bndr == 1);
    
%    imshow(f.*uint8(imcomplement(bndr)))
    
    %centerXY = round(ginput(1));
%    plot(centerXY(1) + 1i* centerXY(2),'g*')
    
%    hold off
    %bndrC = curve2polygon(bndr, centerXY, hbXY, alXY);
    
    save(strcat(dirMaskCropped,name_str,'.mat'), 'bndr', 'indBndr',...
        'centerXY', 'name_str','hbXY','alXY','cXY')
    
    %'bndrC', 'centerXY', 'name_str','hbXY','alXY','cXY')
    
    f = f.*uint8(imcomplement(bndr));
    imwrite(f,strcat(dirMaskCropped,name_str,'.tif'))
    
end

close all