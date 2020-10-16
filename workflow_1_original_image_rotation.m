workflow_0_dir_location

if isDiet
    
    fileList = dir(strcat(dirSource,'G_*.tif'));
else
    
    fileList = dir(strcat(dirSource,'*.tif'));
end

%fileList = dir(dirSource);

L = length(fileList);
fprintf('There are %d files\n',L)


for k = 1:1
    
    name_str = fileList(k).name;
    fprintf('k = %3d\t %s\n',k,name_str(1:end-4))
    
    C1 = strsplit(fileList(k).name,'_');
    
    if isLazy
        [bndr, SegOut, centroidXY] =...
            edge_detection_mutant_lazy(dirSource,name_str);
    else
        [bndr, SegOut, centroidXY, veinX] =...
            edge_detection_mutant(dirSource,name_str);
    end
    
    
    [NX, NY, ~] = size(SegOut);
    
    fig = figure;
    imshow(SegOut)
    title(sprintf('Choose center for the conformal map %s %s %s %s %s %s',...
        C1{1},C1{2},C1{3},C1{4},C1{5},C1{6}(1:end-4)))
    
    hold on
    %plot(centroidXY(1),centroidXY(2),'y*')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    centerXY = round(ginput(1));
    
    
    hold off
    close(fig)
    
    bndrC = curve2polygon(bndr,centerXY);
    
    
    if isLazy
        imwrite(SegOut,strcat(dirCropped,name_str(1:end-4),...
            '_lazy_rotated.tif'))
        save(strcat(dirCropped,name_str(1:end-4),'_lazy.mat'),...
            'bndr','bndrC','centerXY','centroidXY','name_str')
    else
        
        imwrite(SegOut,strcat(dirCropped,name_str(1:end-4),'_rotated.tif'))
        save(strcat(dirCropped,name_str(1:end-4),'.mat'),'bndr','bndrC',...
            'centerXY','centroidXY','name_str','veinX')
    end
    
end
