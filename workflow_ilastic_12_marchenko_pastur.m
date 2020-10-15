% INIT
workflow_0_dir_location

if 1 == 1
    
    %all wings
    dirIn = dirInitiallyAligned;
    %fileList = dir(strcat(dirInitiallyAligned,sprintf('G_*.tif')));
    fileList = dir(strcat(dirInitiallyAligned,sprintf('*.tif')));
    
elseif 1 == 0
    
    %angular
    dirIn = dirGauge;
    fileList = dir(strcat(dirGauge,'*_a_*.tif'));
elseif 1 == 0
    
    %linear
    dirIn = dirGauge;
    fileListNew = dir(strcat(dirGauge,'*.tif'));
    fileListAngular = dir(strcat(dirGauge,'*_a_*.tif'));
    
    fileList = find_new_items(fileListNew,fileListAngular);
end


%load('diet_radon_distance_all_indices.mat')
%load('radon_distance_all_indices.mat')

lFileList = length(fileList);

fprintf('There are %d files\n',lFileList)

img = imread(strcat(dirIn,fileList(1).name));
%indEmpty = find(img == 255);
%img(indEmpty) = 0;
%imagesc(img)

R = radon(img);
%imagesc(R)

%% 0 mutants

iW = 3;
idxW = intersect(idxMutants{iW},idxSex{1});


images = zeros(resolutioN, resolutioN,length(idxW));
radonImages = zeros(size(R,1), size(R,2),length(idxW));


for k = 1:length(idxW)
    
    img = imread(strcat(dirIn,fileList(idxW(k)).name));
    img(indEmpty) = 0;
    img = adapthisteq(img);
    
    images(:,:,k) = img;
    radonImages(:,:,k) = radon(img);
end

qRadon = reshape(radonImages,[size(R,1)*size(R,2) length(idxW)])';
qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));


[~,~,eigVl] = pca(qRadon);

qRadon = reshape(radonImages,[size(R,1)*size(R,2) length(idxW)])';
qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));

%cumSumEigVal = cumsum(eigVl)./sum(eigVl);
%plot(cumSumEigVal ,'o-')

histogram(eigVl,'Normalization','pdf')
hold on

clear('qRadon','eigValue')
nIttr = 30;

for itrr = nIttr:-1:1
    
    fprintf('itteration #%d\n',itrr)
    for k = 1:resolutioN
        for l = 1:resolutioN
            
            idx = randperm(length(idxW));
            
            images(k,l,:) = images(k,l,idx);
        end
    end
    
    
    for k = 1:length(idxW)
        
        radonImages(:,:,k) = radon(images(:,:,k));
    end
    
    
    qRadon = reshape(radonImages,[size(R,1)*size(R,2) length(idxW)])';
    qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
    
    
    [~,~,eigValue(itrr,:)] = pca(qRadon);
    
    eigValue(itrr,:) = sort(eigValue(itrr,:));
    
    
end

figure(2)

histogram(eigValue,'Normalization','pdf')
hold on

%% 1  mutants

iW = 3;
idxW = intersect(idxMutants{iW},idxSex{1});

for nIttr = [200]

%mutants
for iM = 1:-1:1
    
    clear('qRadon','images','radonImages')
     
    idx = intersect(idxMutants{iM},idxSex{1});
    lFileListP = length(idx);
    
    fprintf('There are %d images for %s %s\n',...
        lFileListP,sex{1},mutants{iM})

    
    images = zeros(resolutioN, resolutioN,lFileListP);
    radonImages = zeros(size(R,1), size(R,2),lFileListP);
    
    %fprintf('%dx%dx%d\n',size(radonImages))
    
    for k = 1:lFileListP
        
        img = imread(strcat(dirIn,fileList(idx(k)).name));
        %fprintf('\t \t %s\n',fileList(idx(k)).name)
        
        img(indEmpty) = 0;
        img = adapthisteq(img);
        
        images(:,:,k) = img;
        radonImages(:,:,k) = radon(img);
    end
    
    %fprintf('%dx%dx%d\n',size(radonImages))
    
    qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
%    qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
    
    [~, ~, eigVlM{iM}] = pca(qRadon);
%    figure(iM + 10)
%    histogram(eigVlM{iM},'Normalization','pdf')
%    set(gca,'XScale','log')
    
    figure(iM + nIttr)
    histogram(eigVlM{iM},'Normalization','pdf',...
        'DisplayName',mutants{iM})
    hold on
    
    clear('qRadon')
    
    for itrr = nIttr:-1:1
        
        fprintf('\titteration #%d / %d\n',nIttr-itrr+1,nIttr)
        
        for k = 1:resolutioN
            for l = 1:resolutioN
                
                idxP = randperm(lFileListP);
                
                images(k,l,:) = images(k,l,idxP);
            end
        end
        
        
        for k = 1:lFileListP
            
            radonImages(:,:,k) = radon(images(:,:,k));
        end
        
        
        qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
%        qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
        
        
        [~,~,eigValueM{iM}(itrr,:)] = pca(qRadon);
        
        eigValueM{iM}(itrr,:) = sort(eigValueM{iM}(itrr,:));
    end
    
    
    meanEigVal = mean(eigValueM{iM});
    aM(iM) = meanEigVal(1);
    bM(iM) = meanEigVal(end);
    
    idxSM{iM} = find(eigVlM{iM} > bM(iM));
    lM(iM) = length(idxSM{iM});
    
    fprintf('\t There are %d out of %d eigenwings for %s\n',...
        lM(iM),length(eigValueM{iM}(1,:)), mutants{iM})
    
    figure(iM + nIttr)
    hold on
    histogram(eigValueM{iM},'Normalization','pdf',...
        'DisplayName',sprintf('%d permutations',nIttr))
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    
    
    
    meanEigValD = meanEigVal;
    aeM = aM(iM);
    beM = bM(iM);


    sqrtC = (sqrt(beM) - sqrt(aeM))/(sqrt(beM) + sqrt(aeM));

    %sM = sqrt(bM)/(1+sqrtL);
    sM = (sqrt(beM) + sqrt(aeM))/2;
    %sM = (sqrt(beM) - sqrt(aeM))/2;
    cM = sqrtC^2;


    deltaAm = (beM-aeM)/1000;
    lambdaM = (0.9*aeM:deltaAm:1.01*beM);

    ft = @(lambda,a,b,c,s) ...  
        sqrt((b-lambda).*(lambda-a))./(2*pi*lambda*c*s^2);


    Fm = ft(lambdaM,aeM,beM,cM,sM);

%F = F/sum(F);
    Fm(isnan(Fm)) = 0;
    figure(iM + nIttr)
    hold on
%    plot(lambdaM,Fm,'r--','LineWidth',2,'DisplayName','MP Distribution')

    %title('Eigen Value Distribution')
    title(sprintf('Eigen Value Distribution. %d EigenVectors are well sampled',lM(iM)))
    legend show
    
end

end

%save('data/eig_mutants_100_perm.mat','eigValueM','eigVlM','lM','aM','bM')

%% 1 for Presentation. Diet  Female. High

nIttr = 100;

%temperature
for kT = 3:-1:1
    
    figure(kT)
    hold on
    
    idx = intersect(intersect(idxTemp{kT},idxF),idxHigh);
    %idx = 1:690;
    lFileListP = length(idx);
    
    clear('qRadon','images','radonImages')
    
    
    images = zeros(resolutioN, resolutioN,lFileListP);
    radonImages = zeros(size(R,1), size(R,2),lFileListP);
    
    for k = 1:lFileListP
        
        img = imread(strcat(dirIn,fileList(idx(k)).name));
        img(indEmpty) = 0;
        img = adapthisteq(img);
        
        images(:,:,k) = img;
        radonImages(:,:,k) = radon(img);
    end
    
    qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
    %qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
    
    [~,~,eigVlT{kT}] = pca(qRadon);
    histogram(eigVlT{kT},'Normalization','pdf'  ,...
                        'DisplayName','Female 25C. High diet') 
	hold on
    
    clear('qRadon','eigValueT')
    
    for itrr = nIttr:-1:1
        
        fprintf('itteration #%d\n',itrr)
        
        for k = 1:resolutioN
            for l = 1:resolutioN
                
                idxP = randperm(lFileListP);
                
                images(k,l,:) = images(k,l,idxP);
            end
        end
        
        
        for k = 1:lFileListP
            
            radonImages(:,:,k) = radon(images(:,:,k));
        end
        
        
        qRadon = ...
            reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
        %qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
        
        
        [~,~,eigValueT{kT}(itrr,:)] = pca(qRadon);
        
        eigValueT{kT}(itrr,:) = sort(eigValueT{kT}(itrr,:));
    end
    
    
    meanEigVal = mean(eigValueT{kT});
    aT(kT) = meanEigVal(1);
    bT(kT) = meanEigVal(end);
    
    idxST{kT} = find(eigVlT{kT} > bT(kT));
    lT(kT) = length(idxST{kT});
    
	histogram(eigValueT{kT},'Normalization','pdf',...
                            'DisplayName',sprintf('%d permutations',nIttr))
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    %figure(2)
    %histogram(eigValueT{kT},'Normalization','pdf')
    title(sprintf(...
        'Eigen Value Distribution. There are %d well sampled EigenVectors for %s',...
        lT(kT),temp{kT}))
    legend show
    
    %savefig(sprintf('plots/diet_MP_distribution_female_%s_high.fig',...
    %    temp{kT}))
    %print(sprintf('plots/diet_MP_distribution_female_%s_high.pdf',...
    %    temp{kT}),'-dpdf')
    
end


fprintf('\nThere are %d well sampled eigen vectors for %dC \n',...
        [lT; [18,25,29]])
    
%save(   'data/eig_diet_Female_High_100_perm.mat',...
%        'eigValueT','eigVlT','lT','aT','bT')

%% 1 for Presentation.  Diet  All

nIttr = 100;


%idx = intersect(intersect(idxTemp{kT},idxF),idxHigh);
idx = 1:lFileList;
lFileListP = lFileList;

clear('qRadon','images','radonImages')


images = zeros(resolutioN, resolutioN,lFileListP);
radonImages = zeros(size(R,1), size(R,2),lFileListP);

for k = 1:lFileListP
    
    img = imread(strcat(dirIn,fileList(idx(k)).name));
    img(indEmpty) = 0;
    img = adapthisteq(img);
    
    images(:,:,k) = img;
    radonImages(:,:,k) = radon(img);
end

qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
%qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));

[~,~,eigVlA] = pca(qRadon);
histogram(eigVlA,'Normalization','pdf', 'DisplayName', 'All Conditions')
hold on

clear('qRadon','eigValueT')

for itrr = nIttr:-1:1
    
    fprintf('itteration #%d\n',itrr)
    
    for k = 1:resolutioN
        for l = 1:resolutioN
            
            idxP = randperm(lFileListP);
            
            images(k,l,:) = images(k,l,idxP);
        end
    end
    
    
    for k = 1:lFileListP
        
        radonImages(:,:,k) = radon(images(:,:,k));
    end
    
    
    qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
    %qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
    
    
    [~,~,eigValueA(itrr,:)] = pca(qRadon);
    
    eigValueA(itrr,:) = sort(eigValueA(itrr,:));
end


meanEigVal = mean(eigValueA);
aA = meanEigVal(1);
bA = meanEigVal(end);

idxSA = find(eigVlA > bA);
lA = length(idxSA);

histogram(eigValueA,'Normalization','pdf',...
    'DisplayName',sprintf('%d permutations',nIttr))
set(gca,'YScale','log')
set(gca,'XScale','log')
%figure(2)
%histogram(eigValueT{kT},'Normalization','pdf')
title(sprintf(...
    'Eigen Value Distribution. %d EigenVectors are well sampled', lA));
legend show


%save(   'data/eig_diet_all_100_perm.mat',...
%        'eigValueA','eigVlA','lA','aA','bA')


%% 1

nIttr = 30;

%temperature
for kT = 3:-1:1
    
    idx = idxTemp{kT};
    lFileListP = length(idx);
    
    clear('qRadon','images','radonImages')
    
    
    images = zeros(resolutioN, resolutioN,lFileListP);
    radonImages = zeros(size(R,1), size(R,2),lFileListP);
    
    for k = 1:lFileListP
        
        img = imread(strcat(dirIn,fileList(idx(k)).name));
        img(indEmpty) = 0;
        img = adapthisteq(img);
        
        images(:,:,k) = img;
        radonImages(:,:,k) = radon(img);
    end
    
    qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
    qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
    
    [~,~,eigVlT{kT}] = pca(qRadon);
    %histogram(eigVlT,'Normalization','pdf'), hold on
    
    clear('qRadon','eigValueT')
    
    for itrr = nIttr:-1:1
        
        fprintf('itteration #%d\n',itrr)
        
        for k = 1:resolutioN
            for l = 1:resolutioN
                
                idxP = randperm(lFileListP);
                
                images(k,l,:) = images(k,l,idxP);
            end
        end
        
        
        for k = 1:lFileListP
            
            radonImages(:,:,k) = radon(images(:,:,k));
        end
        
        
        qRadon = ...
            reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
        qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
        
        
        [~,~,eigValueT{kT}(itrr,:)] = pca(qRadon);
        
        eigValueT{kT}(itrr,:) = sort(eigValueT{kT}(itrr,:));
    end
    
    
    meanEigVal = mean(eigValueT{kT});
    aT(kT) = meanEigVal(1);
    bT(kT) = meanEigVal(end);
    
    idxST{kT} = find(eigVlM{kT} > bT(kT));
    lT(kT) = length(idxST{kT});
    %figure(2)
    %histogram(eigValue,'Normalization','pdf')
    
end

save('eig_1_temp.mat','eigValueT','eigVlT')

%Sex
for jS = 2:-1:1
    
    idx = idxSex{jS};
    lFileListP = length(idx);
    
    clear('qRadon','images','radonImages')
    
    
    images = zeros(resolutioN, resolutioN,lFileListP);
    radonImages = zeros(size(R,1), size(R,2),lFileListP);
    
    for k = 1:lFileListP
        
        img = imread(strcat(dirIn,fileList(idx(k)).name));
        img(indEmpty) = 0;
        img = adapthisteq(img);
        
        images(:,:,k) = img;
        radonImages(:,:,k) = radon(img);
    end
    
    qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
    qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
    
    [~,~,eigVlS{jS}] = pca(qRadon);
    %histogram(eigVl,'Normalization','pdf'), hold on
    
    clear('qRadon')
    
    for itrr = nIttr:-1:1
        
        fprintf('itteration #%d\n',itrr)
        
        for k = 1:resolutioN
            for l = 1:resolutioN
                
                idxP = randperm(lFileListP);
                
                images(k,l,:) = images(k,l,idxP);
            end
        end
        
        
        for k = 1:lFileListP
            
            radonImages(:,:,k) = radon(images(:,:,k));
        end
        
        
        qRadon = ...
            reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
        qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
        
        
        [~,~,eigValueS{jS}(itrr,:)] = pca(qRadon);
        
        eigValueS{jS}(itrr,:) = sort(eigValueS{jS}(itrr,:));
    end
    
    
    meanEigVal = mean(eigValueS{jS});
    aS(jS) = meanEigVal(1);
    bS(jS) = meanEigVal(end);
    
    idxSS{jS} = find(eigVlS{jS} > bS(jS));
    lS(jS) = length(idxSS{jS});
    %figure(2)
    %histogram(eigValue,'Normalization','pdf')
    
end

save('eig_1_sex.mat','eigValueS','eigVlS')
%diet we don't plot

%% 2

%Temperature and Sex
for kT = 3:-1:1
    for jS = 2:-1:1
        
        idx = intersect(idxTemp{kT},idxSex{jS});
        lFileListP = length(idx);
        
        fprintf('There are %d wings\n',lFileListP)
        
        clear('qRadon','eigValue','images','radonImages')
        
        
        images = zeros(resolutioN, resolutioN,lFileListP);
        radonImages = zeros(size(R,1), size(R,2),lFileListP);
        
        for k = 1:lFileListP
            
            img = imread(strcat(dirIn,fileList(idx(k)).name));
            img(indEmpty) = 0;
            img = adapthisteq(img);
            
            images(:,:,k) = img;
            radonImages(:,:,k) = radon(img);
        end
        
        qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
        qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
        
        [~,~,eigVlTS{kT,jS}] = pca(qRadon);
        %histogram(eigVl,'Normalization','pdf'), hold on
        
        clear('qRadon','eigValueTS')
        
        for itrr = nIttr:-1:1
            
            fprintf('itteration #%d\n',itrr)
            
            for k = 1:resolutioN
                for l = 1:resolutioN
                    
                    idxP = randperm(lFileListP);
                    
                    images(k,l,:) = images(k,l,idxP);
                end
            end
            
            
            for k = 1:lFileListP
                
                radonImages(:,:,k) = radon(images(:,:,k));
            end
            
            
            qRadon = ...
                reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
            qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
            
            
            [~,~,eigValueTS{kT,jS}(itrr,:)] = pca(qRadon);
            
            eigValueTS{kT,jS}(itrr,:) = sort(eigValueTS{kT,jS}(itrr,:));
        end
        
        
        meanEigVal = mean(eigValueTS{kT,jS});
        aTS(kT,jS) = meanEigVal(1);
        bTS(kT,jS) = meanEigVal(end);
        
        idxSTS{kT,jS} = find(eigVlTS{kT,jS} > bTS(kT,jS));
        lTS(kT,jS) = length(idxSTS{kT,jS});
        %figure(2)
        %histogram(eigValue,'Normalization','pdf')
        
    end
end


save('eig_2_temp_sex.mat','eigValueTS','eigVlTS')
% we don't include diet


%% 3

%Temperature and Sex

for kT = 3:-1:1
    for jS = 2:-1:1
        for jD = 2:-1:1
            
            idx = intersect(intersect(idxTemp{kT},idxSex{jS}),...
                idxDiet{jD});
            lFileListP = length(idx);
            
            if lFileListP == 0
                
                fprintf('skip. T = %s, Sex = %s, Diet = %s\n',...
                    temp{kT},sex{jS}, diet{jD})
                break
            else
                
                fprintf('T = %s, Sex = %s, Diet = %s\n',...
                    temp{kT},sex{jS}, diet{jD})
            end
            
            fprintf('There are %d wings\n',lFileListP)
            
            
            
            clear('qRadon','eigValue','images','radonImages')
            
            
            images = zeros(resolutioN, resolutioN,lFileListP);
            radonImages = zeros(size(R,1), size(R,2),lFileListP);
            
            for k = 1:lFileListP
                
                img = imread(strcat(dirIn,fileList(idx(k)).name));
                img(indEmpty) = 0;
                img = adapthisteq(img);
                
                images(:,:,k) = img;
                radonImages(:,:,k) = radon(img);
            end
            
            qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
            qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
            
            [~,~,eigVlTSD{kT,jS,jD}] = pca(qRadon);
            %histogram(eigVl,'Normalization','pdf'), hold on
            
            clear('qRadon')
            
            for itrr = nIttr:-1:1
                
                fprintf('itteration #%d\n',itrr)
                
                for k = 1:resolutioN
                    for l = 1:resolutioN
                        
                        idxP = randperm(lFileListP);
                        
                        images(k,l,:) = images(k,l,idxP);
                    end
                end
                
                
                for k = 1:lFileListP
                    
                    radonImages(:,:,k) = radon(images(:,:,k));
                end
                
                
                qRadon = ...
                    reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
                qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
                
                
                [~,~,eigValueTSD{kT,jS,jD}(itrr,:)] = pca(qRadon);
                
                eigValueTSD{kT,jS,jD}(itrr,:) =...
                    sort(eigValueTSD{kT,jS,jD}(itrr,:));
            end
            
            
            meanEigVal = mean(eigValueTSD{kT,jS,jD});
            aTSD(kT,jS,jD) = meanEigVal(1);
            bTSD(kT,jS,jD) = meanEigVal(end);
            
            idxSTSD{kT,jS,jD} = find(eigVlTSD{kT,jS,jD} > bTSD(kT,jS,jD));
            lTSD(kT,jS,jD) = length(idxSTSD{kT,jS,jD});
            %figure(2)
            %histogram(eigValue,'Normalization','pdf')
            
        end
    end
end


save('eig_3_temp_sex_diet.mat','eigValueTSD','eigVlTSD')



%% kernel 3
addpath('kPCA_v3.1/code');

%Temperature and Sex

for kT = 3:-1:1
    for jS = 2:-1:1
        for jD = 2:-1:1
            
            idx = intersect(intersect(idxTemp{kT},idxSex{jS}),...
                idxDiet{jD});
            lFileListP = length(idx);
            
            if lFileListP == 0
                
                fprintf('skip. T = %s, Sex = %s, Diet = %s\n',...
                    temp{kT},sex{jS}, diet{jD})
                break
            else
                
                fprintf('T = %s, Sex = %s, Diet = %s\n',...
                    temp{kT},sex{jS}, diet{jD})
            end
            
            fprintf('There are %d wings\n',lFileListP)
            
            
            
            clear('qRadon','eigValue','images','radonImages')
            
            
            images = zeros(resolutioN, resolutioN,lFileListP);
            radonImages = zeros(size(R,1), size(R,2),lFileListP);
            
            for k = 1:lFileListP
                
                img = imread(strcat(dirIn,fileList(idx(k)).name));
                img(indEmpty) = 0;
                img = adapthisteq(img);
                
                images(:,:,k) = img;
                radonImages(:,:,k) = radon(img);
            end
            
            qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
            qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
            
            %[~,~,eigVlTSD{kT,jS,jD}] = pca(qRadon);
            para = 5;
            [~,~,kEigVlTSD{kT,jS,jD}] = kPCA(qRadon,2,'poly',para);
            [~,~,kEigVlTSDG{kT,jS,jD}] = kPCA(qRadon,2,'gaussian',para);
            %histogram(kEigVlTSD{kT,jS,jD},'Normalization','pdf'), hold on
            
            %cumSumEigVal = ...
            %    cumsum(kEigVlTSD{kT,jS,jD})./sum(kEigVlTSD{kT,jS,jD});
            %plot(cumSumEigVal ,'o-')
            
           % cumSumEigValG = ...
            %    cumsum(kEigVlTSDG{kT,jS,jD})./sum(kEigVlTSDG{kT,jS,jD});
            %plot(cumSumEigValG ,'o-')
            
            clear('qRadon')
            
            for itrr = nIttr:-1:1
                
                fprintf('itteration #%d\n',itrr)
                
                for k = 1:resolutioN
                    for l = 1:resolutioN
                        
                        idxP = randperm(lFileListP);
                        
                        images(k,l,:) = images(k,l,idxP);
                    end
                end
                
                
                for k = 1:lFileListP
                    
                    radonImages(:,:,k) = radon(images(:,:,k));
                end
                
                
                qRadon = ...
                    reshape(radonImages,[size(R,1)*size(R,2) lFileListP])';
                qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
                
                
                %[~,~,eigValueTSD{kT,jS,jD}(itrr,:)] = pca(qRadon);
                [~,~,kEigValueTSD{kT,jS,jD}(itrr,:)] = ...
                    kPCA(qRadon,2,'poly',para);
                
                kEigValueTSD{kT,jS,jD}(itrr,:) =...
                    sort(kEigValueTSD{kT,jS,jD}(itrr,:));
            end
            
            
            meanEigVal = mean(kEigValueTSD{kT,jS,jD});
            aTSDk(kT,jS,jD) = meanEigVal(1);
            bTSDk(kT,jS,jD) = meanEigVal(end);
            
            idxSTSDk{kT,jS,jD} = find(kEigVlTSD{kT,jS,jD} > bTSDk(kT,jS,jD));
            lTSDk(kT,jS,jD) = length(idxSTSDk{kT,jS,jD});
            %figure(2)
            %histogram(eigValue,'Normalization','pdf')
            
        end
    end
end

save('kernel_eig_3_temp_sex_diet.mat','kEigValueTSD','kEigVlTSD')

%%


meanEigVal = mean(eigValue);
aM = meanEigVal(1);
bM = meanEigVal(end);


sqrtC = (sqrt(bM) - sqrt(aM))/(sqrt(bM) + sqrt(aM));

%sM = sqrt(bM)/(1+sqrtL);
sM = (sqrt(bM) + sqrt(aM))/2;

cM = sqrtC^2;

a = aM;
b = bM;
c = cM;
s = sM;


deltaA = (b-a)/1000;
lambda = (0.9*a:deltaA:1.05*b);

ft = @(lambda,a,b,c,s) ...
    sqrt((b-lambda).*(lambda-a))./(2*pi*lambda*c*s^2);

F = ft(lambda,a,b,c,s);

%F = F/sum(F);
F(isnan(F)) = 0;



plot(lambda,F,'y','LineWidth',2)



%% Marchenko-Pastur works well
n = 100;
m = 200;

c = n/m;

sigma = 1;
%sigma = 0.5;
%R = sqrt(sigma);

nIttr = 100;
clear('D','aI','bI','sI')

for itrr = nIttr:-1:1
    
    %A = (rand(n,m)-0.5)*R;
    A = randn(n,m)*sigma;
    
    sI(itrr) = std(A(:));
    %s = R;
    
    
    D(itrr,:) = eig(A*(A')/m);
    
    D(itrr,:) = sort(D(itrr,:));
    
    aI(itrr) = (sI(itrr)^2)*(1-sqrt(c))^2;
    bI(itrr) = (sI(itrr)^2)*(1+sqrt(c))^2;
    
    
end

%a = mean(aI);
%b = mean(bI);
%s = mean(sI);

a = sigma^2*(1 - sqrt(c))^2;
b = sigma^2*(1 + sqrt(c))^2;
s = sigma;

deltaA = (b-a)/1000;
lambda = (0.9*a:deltaA:1.05*b);

ft = @(lambda,a,b,c,s) ...
    sqrt((b-lambda).*(lambda-a))./(2*pi*lambda*c*s^2);

F = ft(lambda,a,b,c,s);

%F = F/sum(F);
F(isnan(F)) = 0;
figure(4321)
histogram(D(:),'Normalization','pdf',...
    'DisplayName','Eigen Value Distribution')

hold on

plot(lambda,F,'y','LineWidth',2,'DisplayName','Theoretical Curve')


meanEigValD = mean(D);
aM = meanEigValD(1);
bM = meanEigValD(end);


sqrtC = (sqrt(bM) - sqrt(aM))/(sqrt(bM) + sqrt(aM));

%sM = sqrt(bM)/(1+sqrtL);
sM = (sqrt(bM) + sqrt(aM))/2;

cM = sqrtC^2;


deltaAm = (bM-aM)/1000;
lambdaM = (0.9*aM:deltaAm:1.01*bM);


Fm = ft(lambdaM,aM,bM,cM,sM);

%F = F/sum(F);
Fm(isnan(Fm)) = 0;
figure(4321)
hold on
plot(lambdaM,Fm,'r--','LineWidth',2,'DisplayName','Experimental Curve')



legend show

%% Experiment with Block-diagonal matrices 


n1 = 100;
m1 = 900;

c1 = n1/m1;

sigma1 = 3;

n2= 490;
m2 = 500;

c2 = n2/m2;

sigma2 = 1;

nIttr = 100;
clear('D','aI','bI','sI')

for itrr = nIttr:-1:1
    
    %A = (rand(n,m)-0.5)*R;
    A = randn(n1,m1)*sigma1;
    
    B = randn(n2,m2)*(sigma2);
    
    sI(itrr) = std(A(:));
    %s = R;
    
    Mat = blkdiag(A,B);
    
    mMat = size(Mat,2);
    
    D(itrr,:) = eig(Mat*(Mat')/mMat);
    
    D(itrr,:) = sort(D(itrr,:));
    
    aI(itrr) = (sI(itrr)^2)*(1-sqrt(c))^2;
    bI(itrr) = (sI(itrr)^2)*(1+sqrt(c))^2;
    
    
end

%a = mean(aI);
%b = mean(bI);
%s = mean(sI);

a = sigma^2*(1 - sqrt(c))^2;
b = sigma^2*(1 + sqrt(c))^2;
s = sigma;

deltaA = (b-a)/1000;
lambda = (0.9*a:deltaA:1.05*b);

ft = @(lambda,a,b,c,s) ...
    sqrt((b-lambda).*(lambda-a))./(2*pi*lambda*c*s^2);

F = ft(lambda,a,b,c,s);

%F = F/sum(F);
F(isnan(F)) = 0;
figure(432)
histogram(D(:),'Normalization','pdf',...
    'DisplayName','Eigen Value Distribution')

hold on

%plot(lambda,F,'y','LineWidth',2,'DisplayName','Theoretical Curve')


meanEigValD = mean(D);
aM = meanEigValD(1);
bM = meanEigValD(end);


sqrtC = (sqrt(bM) - sqrt(aM))/(sqrt(bM) + sqrt(aM));

%sM = sqrt(bM)/(1+sqrtL);
sM = (sqrt(bM) + sqrt(aM))/2;

cM = sqrtC^2;


deltaAm = (bM-aM)/1000;
lambdaM = (0.9*aM:deltaAm:1.01*bM);


Fm = ft(lambdaM,aM,bM,cM,sM);

%F = F/sum(F);
Fm(isnan(Fm)) = 0;
figure(432)
hold on
%plot(lambdaM,Fm,'r--','LineWidth',2,'DisplayName','Experimental Curve')



legend show

%% Marchenko Pastur Distribution 

% Marchenko Pastur Distribution
% In Random Matrix Theory, MP law gives the probability density function
% of singular values of large rectangular random matrices;
% when the dimensions of matrix tend to infinity.
% This contribution illustrates the PDF of matrix Y(N,N)= X*X^T/M,
% where X is random matrix whose entries X_i,j are independent
% and identically distributed random variables with zero mean
% and variance s^2. The program is applicable for both uniform and random
% distributions.
% Ref :
% Marchenko,V. A., Pastur, L. A. (1967) "Distribution of eigenvalues
% for some sets of random matrices", Mat. Sb. (N.S.), 72(114):4, 507?536

clear;
close all;
N = 200;
M = 400;
% Ratio of matrix dimensions
c = N/M;
% Sample
sigma = 1;
x = randn(N,M)*sigma; % Normal distribution with sigma
%x=rand(N,T);  % Uniform distribution
s = std(x(:));
% spectral matrix
r = x*x'/M;
%eigenvalues
l = eig(r);

% Probability Density Function

% number of points for measurement.

n = 30; % number of bins

% Boundaries
a = (s^2)*(1-sqrt(c))^2;
b = (s^2)*(1+sqrt(c))^2;
[f,lambda] = hist(l,linspace(a,b,n));
% Normalization
f = f/sum(f);
% Theoretical pdf
ft = @(lambda,a,b,c) ...
    (1./(2*pi*lambda*c*s^(2))).*sqrt((b-lambda).*(lambda-a));
F = ft(lambda,a,b,c);
% Processing numerical pdf
F = F/sum(F);
F(isnan(F))=0;
% Results
figure(1234)
h = bar(lambda,f);
set(h,'FaceColor',[.75 .75 .8]);
set(h,'LineWidth',0.25);
xlabel('Eigenvalue \lambda');
ylabel(' Probability Density Function f(\lambda)');
title(' Marchenko-Pastur distribution');
lmin=min(l);
lmax=max(l);
%axis([-1 2*lmax 0 max(f)+max(f)/4]);
hold on;
plot(lambda,F,'g','LineWidth',2);
hold off;
