workflow_0_dir_location

%fileList = dir(strcat(dirImage,'*.tif'));

fileListNew = dir(strcat(dirImage,'*.tif'));
fileListDone = dir(strcat(dirAngularAligned,'*.tif'));

fileList = find_new_items(fileListNew,fileListDone);

%fileList = fileListNew;

L = length(fileList);
fprintf('There are %d files\n',L)

angle = [-170:0.5:170];

maxAngle = zeros(L,1);


if isDiet

    imgRef = imread(strcat(dirImage,'G_18C_high_F_L_1_21.tif'));
elseif isInbred
    
    imgRef = imread(strcat(dirImage,'G_25203_25C_01_L.tif'));
else
    
    
    imgRef = imread(strcat(dirImage,'samw_F_L_lei_4X_110.tif'));
end
%imgRef = imread(strcat(dirImage,'cognata M1a z16 s103.tif'));
indWhite = find(imgRef == 255);

angleBest = zeros(1,L);
ccBest = zeros(1,L);
tic;

for indK = 1:L
    
    tStart = toc;
    name_str = fileList(indK).name(1:end-4);
    img = imread(strcat(dirImage,fileList(indK).name));
    %img= imrotate255(imgRef,7.5,'bicubic','crop');
    %img(indWhite) = 255;
    %img = imgBest;
    [imgBest,angleBest(indK), ccBest(indK)] =...
        perform_angular_alignment(imgRef,img,indWhite,angle);
    
    %imshow(imabsdiff(imgRef,imgBest))
    
    imwrite(imgBest,strcat(dirAngularAligned,name_str,'.tif'))
   
    
    tStop = toc;
    fprintf('#%4d, angle = %4.2f || %s\t|| %2.2f sec\n', ...
        indK, angleBest(indK),name_str,tStop-tStart)
    
    
end

save(strcat(dirAngularAligned,'best_angles','.mat'),...
    'angleBest', 'ccBest','fileList')
