workflow_0_dir_location

if isLazy
    fileList = dir(strcat(dirCropped,'*_lazy.mat'));
else
    fileList = dir(strcat(dirCropped,'*.mat'));
end

L = length(fileList);
fprintf('There are %d files\n',L)

for k =  1:L
    
    load(strcat(dirCropped,fileList(k).name),'bndrC','name_str')
    
    tic;
    tStart = toc;
    
    bndrX = real(bndrC);
    bndrY = imag(bndrC);
    
    xMin = min(bndrX);
    xMax = max(bndrX);
    
    nMax = find(bndrX == xMax);
    nMin = find(bndrX == xMin);
    
    pSmooth = 0.01;
    
    % upper part
    %curvSmooth = csaps(bndrX(nMax:end),bndrY(nMax:end),pSmooth);
    curvSmooth = csaps(bndrX(1:nMin),bndrY(1:nMin),pSmooth);
    
    % lower part
    %curvSmooth2 = csaps(bndrX(1:nMax),bndrY(1:nMax),pSmooth);
    curvSmooth2 = csaps(bndrX(nMin:end),bndrY(nMin:end),pSmooth);
    
    %figure
    
    %plot(centerXY(1) + i * centerXY(2),'r*')
    %hold on
    %plot(centerXY2(1) + i * centerXY2(2),'g*')
    
    %plot(bndrX,bndrY,'g.')
    %hold on
    %fnplt(curvSmooth,'r')
    %fnplt(curvSmooth2,'b')
    
    %axis equal
    %xlim([1.1*xMin 1.1*xMax])
    
    %plot(skeleton,'.')
    %hold off
    
    
    delta = 50;
    curV =...
        LineCurvature2D([xMin:delta:xMax;fnval(curvSmooth,xMin:delta:xMax)]');
    curvSpline = csaps(xMin:delta:xMax,abs(curV),0.001);
    curvF = @(x) fnval(curvSpline,x);
    normConst = integral(curvF,xMin,xMax);
    
    curV2 = ...
        LineCurvature2D([xMin:delta:xMax;fnval(curvSmooth2,xMin:delta:xMax)]');
    curvSpline2 = csaps(xMin:delta:xMax,abs(curV2),0.001);
    curvF2 = @(x) fnval(curvSpline2,x);
    normConst2 = integral(curvF2,xMin,xMax);
    
    
    
    %figure
    %hold on
    %fnplt(curvSpline)
    %fnplt(curvSpline2)
    
    %legend('upper part','lower part')
    
    tau = ceil(xMin):floor(xMax);
    
    
    n = 1/((N-1));
    
    curvInt = zeros(1,length(tau));
    curvInt2 = zeros(1,length(tau));
    
    for r = 1:length(tau)
        
        curvInt(r) = integral(curvF,xMin,tau(r))/normConst;
        curvInt2(r) = integral(curvF2,xMin,tau(r))/normConst2;
    end
    
    
    %figure
    
    %plot(tau,curvInt)
    %hold on
    %plot(tau,curvInt2)
    
    %xlim([1.1*xMin 1.1*xMax])
    
    r1 = [tau(1)];
    r2 = [tau(1)];
    
    n1 = 1;
    n2 = 1;
    
    for r = 1:length(tau)
        if n1 * n <= curvInt(r)
            r1(end+1) = tau(r);
            n1 = n1 + 1;
        end
        
        if n2 * n <= curvInt2(r)
            r2(end+1) = tau(r);
            n2 = n2 + 1;
        end
    end
    
    r1(end+1) = tau(r);
    r2(end+1) = tau(r);
    
    r1 = unique(r1);
    
    r1 = r1(end:-1:1);
    r2 = unique(r2);
    
    %figure(1)
    %hold on
    
    %plot(r1,fnval(curvSmooth,r1),'*')
    %plot(r2,fnval(curvSmooth2,r2),'*')
    
    bndrUP = r1 + 1i*fnval(curvSmooth,r1);
    bndrLOWER = r2 + 1i*fnval(curvSmooth2,r2);
    
    
    bndrRefined = [bndrUP(1:end-1),(bndrLOWER(1) + bndrUP(end))/2,bndrLOWER(2:end-1)];
    
    
    %figure
    %plot(bndrRefined,'*-')
    %axis equal
    
    %pC = polygon([bndrRefined]);
    %pC = polygon([bndrRefined(end:-1:1)]);
    fprintf('k = %d, length of polygon = %d vertices. %s\n',...
        k, length(bndrRefined),fileList(k).name(1:end-4))
    
    if isLazy
        
        save(strcat(dirBoundary,name_str(1:end-4),'_lazy.mat'),...
            'bndrRefined','name_str')
    else
        
        save(strcat(dirBoundary,name_str(1:end-4),'.mat'),...
            'bndrRefined','name_str')
    end
    clear('bndrC','name_str')
end