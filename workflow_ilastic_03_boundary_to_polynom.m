workflow_0_dir_location

fileListDone = dir(strcat(dirBoundary,'*.mat'));
fileListNew = dir(strcat(dirMaskCropped,'*.mat'));

fileList = find_new_items(fileListNew,fileListDone);

fileList  = fileListNew; 

lFileList = length(fileList);
fprintf('There are %d files\n',lFileList)
tic;

%%
for indK =  253:lFileList%min(10,lFileList)
    
    tStart = toc;
    
    name_str = fileList(indK).name(1:end-4);
    
    load(strcat(dirMaskCropped,name_str,'.mat'), 'bndr','centerXY',...
        'hbXY', 'alXY')
    
    if isDiet
        
        zz = hbXY;
        hbXY = alXY;
        alXY = zz;
    elseif isInbred
        
        zz = hbXY;
        hbXY = alXY;
        alXY = zz;
    end
    
    % A stupid way around some technical difficulties.
    % It is necessary to remove 4 parts and restore two.
    
    bndrFill = imfill(255.*uint8(bndr));
    [Label, ~] = bwlabel(bndrFill);
    D = regionprops(Label, 'Area','Orientation');
    
    angle = D.Orientation + 30;
    Rot = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];
    
    bndrFillRotated = imrotate(bndrFill,-angle,'loose', 'bilinear');
    bndrRotated = bwperim(bndrFillRotated);
    
    [lY, lX] = size(bndr);
    [lY1, lX1] = size(bndrRotated);
    
    
    centerXYrotated = round(Rot*[(centerXY(1) - lX/2),...
        (centerXY(2) -lY/2)]' +[lX1/2;lY1/2])';
    
    hbXYrotated = round(Rot*[(hbXY(1) - lX/2),...
        (hbXY(2) -lY/2)]' +[ lX1/2;lY1/2]);
    
    alXYrotated = round(Rot*[(alXY(1) - lX/2),...
        (alXY(2) -lY/2)]' +[lX1/2; lY1/2]);
    
    %bndrC = curve2polygon(bndr, centerXY, hbXY, alXY);
    bndrC = curve2polygon(bndrRotated, centerXYrotated,...
        hbXYrotated, alXYrotated);
    
    bndrX = real(bndrC);
    bndrY = imag(bndrC);
    
    xMin = min(bndrX);
    nMin = find(bndrX == xMin,1);
    
    xMax1 = bndrX(1);
    nMax1 = 1;
    
    xMax2 = bndrX(end);
    nMax2 = length(bndrX);
    
    
    % We would like to smooth boundary. We use cubic splines.
    % There are at most four part of the boundary. The first interval is
    % from humeral break to xMax1( doesn't always exist). The second
    % interval exists always, it starts at xMax and ends at xMin. The third
    % interval exists always, it starts at xMin and ends at xMax2. The
    % fourth interval doesn't always exist. It starts at xMax2 and ends at
    % alula notch.
    
    pSmooth = 0.001;
    
    %first part
    if nMax1 > 2
        
        curvSmooth1 = csaps(bndrX(1:nMax1),bndrY(1:nMax1),pSmooth);
    end
    
    %second part
    curvSmooth2 = csaps(bndrX(nMax1:nMin),bndrY(nMax1:nMin),pSmooth);
    
    %third part
    %curvSmooth3 = csaps(bndrX(nMin:nMax2),bndrY(nMin:nMax2),pSmooth);
    if nMax2 - nMin > 2
        curvSmooth3 = csaps(bndrX(nMin:nMax2),bndrY(nMin:nMax2),pSmooth);
    end
    
    %forth part
    if nMax2 < length(bndrC) - 1
        
        curvSmooth4 = csaps(bndrX(nMax2+1:end),bndrY(nMax2+1:end),pSmooth);
    end
    
    
    %fnplt(curvSmooth1,'r')
    %plot([c3(1,:)+centerXY(1)],[c3(2,:)+(size(f,1)-centerXY(2))])
    
    
    delta = 30;
    % the smaller pSmoot the more smooth the curve is.
    pSmoothCurvature = 0.2;
    
    
    %first part
    if nMax1 > 2*delta + 1
        
        interval1 = bndrX(1):delta:xMax1;
    elseif nMax1 > 2 && nMax1 <= 2*delta + 1
        
        interval1 = [bndrX(1):xMax1];
    else
        
        interval1  = [];
    end
    
    
    if nMax2 - nMin > 2
        
        interval3 = unique([xMin:delta:xMax2,xMax2]);
    else
        interval3 = [];
    end
    
    
    if ~isempty(interval1)
        
        curV1 = LineCurvature2D([interval1;fnval(curvSmooth1,interval1)]');
        curvSpline1 = csaps(interval1,abs(curV1),pSmoothCurvature);
        curvF1 = @(x) fnval(curvSpline1,x);
        %normConst = integral(curvF1,interval1(1),interval1(end));
    end
    
    
    %second part
    interval2 = unique([xMin:delta:xMax1,xMax1]);
    curV2 = LineCurvature2D([interval2;fnval(curvSmooth2,interval2)]');
    curvSpline2 = csaps(interval2,abs(curV2),pSmoothCurvature);
    curvF2 = @(x) fnval(curvSpline2,x);
    %normConst2 = integral(curvF2,interval2(1),interval2(end));
    
    
    %third part
    %interval3 = unique([xMin:delta:xMax2,xMax2]);

    if ~isempty(interval3)
        
        curV3 = LineCurvature2D([interval3;fnval(curvSmooth3,interval3)]');
        curvSpline3 = csaps(interval3,abs(curV3),pSmoothCurvature);
        curvF3 = @(x) fnval(curvSpline3,x);
    %normConst3 = integral(curvF3,interval3(1),interval3(end));
    end
    
    %fourth part
    if nMax2 < length(bndrC)- 2*delta
        
        interval4 = bndrX(end):delta:xMax2;
    elseif nMax2 >= length(bndrC)- 2*delta && nMax2 < length(bndrC)-1
        
        interval4 = [bndrX(end):xMax2];
    else
        interval4 = [];
    end
    
    if ~isempty(interval4)
        
        curV4 = LineCurvature2D([interval4;fnval(curvSmooth4,interval4)]');
        curvSpline4 = csaps(interval4,abs(curV4),pSmoothCurvature);
        curvF4 = @(x) fnval(curvSpline4,x);
        %normConst4 = integral(curvF4,interval4(1),interval4(end));
    end
    
    %figure
    %hold on
    %fnplt(curvSpline)
    %fnplt(curvSpline2)
    
    %legend('upper part','lower part')
    
    if nMax1 > 2
        tau1 = bndrX(1):xMax1;
    else
        tau1 = [];
    end
    tau2 = ceil(xMin):floor(xMax1);
    
    if nMax2 - nMin > 2
        tau3 = ceil(bndrX(nMin+1)):floor(xMax2);
    else
        tau3 = [];
    end
    
    
    if nMax2 < length(bndrC) - 1
        
        tau4 = bndrX(end):xMax2;
    else
        tau4 = [];
    end
    
    
    n = 1/((N-1));
    
    curvInt1 = zeros(1,length(tau1));
    curvInt2 = zeros(1,length(tau2));
    curvInt3 = zeros(1,length(tau3));
    curvInt4 = zeros(1,length(tau4));
    %comulative curvature
    
    for r = 1:length(tau1)
        
        curvInt1(r) = integral(curvF1,interval1(1),tau1(r));
    end
    
    for r = 1:length(tau2)
        curvInt2(r) = integral(curvF2,interval2(1),tau2(r));
    end
    
    for r = 1:length(tau3)
        curvInt3(r) = integral(curvF3,interval3(1),tau3(r));
    end
    
    for r = 1:length(tau4)
        curvInt4(r) = integral(curvF4,interval4(1),tau4(r));
    end
    
    %%
    normConstFull = 0;
    
    if nMax1 > 2
        normConstFull = normConstFull + curvInt1(end);
    end
    
    normConstFull = normConstFull + curvInt2(end) + curvInt3(end);
    
    if nMax2 < length(bndrC) - 1
        normConstFull = normConstFull + curvInt4(end);
        nPoints4 = ceil(N*curvInt4(end)/normConstFull);
        if nPoints4 < 2
            nPoints4 = 2;
        end
        
        curvInt4 = curvInt4./curvInt4(end);
    else
        
        nPoints4 = 0;
    end
    
    if nMax1 > 2
        
        nPoints1 = ceil(N*curvInt1(end)/normConstFull);
        if nPoints1 < 2
            nPoints1 = 2;
        end
        curvInt1 = curvInt1./curvInt1(end);
    else
        
        nPoints1 = 0;
    end
    
    nPoints2 = ceil((N-nPoints1-nPoints4)*curvInt2(end)/normConstFull);
    nPoints3 = N - nPoints1 - nPoints2 - nPoints4 + 1;
    
    curvInt2 = curvInt2./curvInt2(end);
    curvInt3 = curvInt3./curvInt3(end);

    
    r1 = zeros(1,nPoints1);
    r2 = zeros(1,nPoints2);
    r3 = zeros(1,nPoints3);
    r4 = zeros(1,nPoints4);
    
    if nMax1 > 2
        r1(1) = tau1(1);
    else
        r1 = [];
    end
    
    r2(1) = tau2(1);
    r3(1) = tau3(1);
    
    if nMax2 < length(bndrC) - 1
        r4(1) = tau4(1);
    else
        r4 = [];
    end
    
    n1 = 1;
    n2 = 1;
    n3 = 1;
    n4 = 1;
    
    for r = 1:length(tau1)
        
        if n1/(nPoints1-1) <= curvInt1(r)
            
            n1 = n1 + 1;
            r1(n1) = tau1(r);
        end
    end
    
    for r = 1:length(tau2)
        
        if n2/(nPoints2-1) <= curvInt2(r)
            
            n2 = n2 + 1;
            r2(n2) = tau2(r);
        end
    end
    
    for r = 1:length(tau3)
        
        if n3/(nPoints3) <= curvInt3(r)
            
            n3 = n3 + 1;
            r3(n3) = tau3(r);
        end
    end
    
    for r = 1:length(tau4)
        
        if n4/(nPoints4-1) <= curvInt4(r)
            
            n4 = n4 + 1;
            r4(n4) = tau4(r);
        end
    end
    
    %r1(end+1) = tau(r);
    %r2(end+1) = tau(r);
    
    r1 = unique(r1,'stable');
    r2 = unique(fliplr(r2),'stable');
    r3 = unique(r3,'stable');
    r4 = unique(fliplr(r4),'stable');
    
    
    %figure(1)
    %hold on
    
    %plot(r1,fnval(curvSmooth,r1),'*')
    %plot(r2,fnval(curvSmooth2,r2),'*')
    if nMax1 > 2
        bndr1 = r1 + 1i*fnval(curvSmooth1,r1);
    else
        bndr1 =[];
    end
    
    bndr2 = r2 + 1i*fnval(curvSmooth2,r2);
    bndr3 = r3 + 1i*fnval(curvSmooth3,r3);
    
    if nMax2 < length(bndrC) - 1
        bndr4 = r4 + 1i*fnval(curvSmooth4,r4);
    else
        bndr4 =[];
    end
    
    
    bndrRefined = ...
        [bndr2(2:end-1),(bndr2(end) + bndr3(1))/2,bndr3(2:end-1)];
    
    if nMax1 > 2
        bndrRefined =...
            [bndr1(1:end-1),(bndr1(end) + bndr2(1))/2, bndrRefined];
    else
        bndrRefined =[bndr2(1),bndrRefined];
    end
    
    if nMax2 < length(bndrC) - 1
        bndrRefined =[bndrRefined,(bndr3(end)+bndr4(1))/2,bndr4(2:end)];
    else
        bndrRefined =[bndrRefined,bndr3(end)];
    end
    
    %figure
    %plot(bndrRefined,'*-')
    %axis equal
    
    if isDiet || isInbred
        
        f = im2uint8(imread(strcat(dirSource,name_str,'.tif')));
    else
        
        f = rgb2gray(im2uint8(imread(strcat(dirSource,name_str,'.tif'))));
    end
    name_cell = strsplit(name_str,'_');
    
    if isMutant
        if strcmp(name_cell{3},'R')
            
            f = fliplr(f);
        end
    end
    
    B = imrotate(f, -angle, 'loose');
    
%    plot(	[real(bndrRefined)+centerXY(1)],...
%           [-imag(bndrRefined) + centerXY(2)],...
%            'LineWidth', 2                                  );
    areaC = round(polyarea(real(bndrC),imag(bndrC)));
    
    imwrite(B,strcat(dirSourceRotated,name_str,'.tif'))
    
    save(strcat(dirBoundary,name_str,'.mat'),'bndrC', 'bndrRefined',...
        'name_str', 'centerXYrotated', 'hbXYrotated', 'alXYrotated',...
        'areaC', 'angle')
    
    tStop = toc;
    
    fprintf('k = %3d, %d vertices| %s | Area = %7d || %2.2f sec\n',...
        indK, length(bndrRefined),name_str,areaC,tStop-tStart)
    
    
    clear('bndrC','name_str')
end