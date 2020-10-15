%% Init
workflow_0_dir_location

fileList = dir(strcat(dirInitiallyAligned,'*.tif'));
lFileList = length(fileList);


labels = cell(lFileList,1);
%L = zeros(lFileList,1);
nameWT = 'samw';
mutants = {'egfr','mam','samw','star','tkv'};
sex = {'F','M'};
chirality = {'L','R'};

img = imread(strcat(dirVeinAligned,fileList(1).name));
indEmpty = find(img == 0);
R = radon(img);

imgRef = imread(strcat(dirVeinAligned,'samw_F_L_lei_4X_110.tif'));
imgRef(indEmpty) = 0;
%histRef = imhist(imgRef);

images = zeros(size(img,1), size(img,2),lFileList);
pImages = zeros(size(img,1), size(img,2),15);
qImages = zeros(size(img,1), size(img,2),15);
%pImg = zeros(size(img,1)*size(img,2),15);
%pqImg = zeros(size(img,1)*size(img,2),size(img,1)*size(img,2),15);

load('radon_distance_all_indices.mat')

for k = 1:lFileList
    
    name_str{k} = fileList(k).name;
    
    img = imread(strcat(dirVeinAligned,name_str{k}));
    img(indEmpty) = 0;
    %img = histeq(img, histRef);
    %img = adapthisteq(img);
    %img  = imgradient(img);
    
    images(:,:,k)  = img;
end

%% alternative Radon

for k = lFileList:-1:1
    
    img = 255-images(:,:,k);
    img(indEmpty) = 0 ;
    imagesInv(:,:,k) = img;
end


for k = lFileList:-1:1
    radonImages(:,:,k) = radon(imagesInv(:,:,k));
end

qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileList])';
%qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));


%% RADON


radonImages = zeros(size(R,1), size(R,2),lFileList);

for k = 1:lFileList
    
    radonImages(:,:,k) = radon(images(:,:,k));
end


qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileList])';
qRadon1 = qRadon;
%qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
%qRadon = bsxfun(@rdivide,qRadon,sum(qRadon,2));

Z = tsne(qRadon,'Algorithm','exact','Distance','euclidean');


for k = 1:lFileList
    
    %labels{k} = strcat(name_cell{k}{1},'_',name_cell{k}{2},'_',...
    %    name_cell{k}{3});
    %labels{k} = strcat(name_cell{k}{2},'_',name_cell{k}{3});
    labels{k} = strcat(name_cell{k}{1});
end

L(strcmp(labels,'egfr')) = 1;
L(strcmp(labels,'mam')) = 2;
L(strcmp(labels,'samw')) = 3;
L(strcmp(labels,'star')) = 4;
L(strcmp(labels,'tkv')) = 5;


figure
scatter(Z(:,1),Z(:,2),[],L,'filled')

%% playing with centroids 

idx1 = intersect(idxMutants{1},idxSex{2});
idx2 = intersect(idxMutants{2},idxSex{2});
idx3 = intersect(idxMutants{3},idxSex{2});
idx4 = intersect(idxMutants{4},idxSex{2});
idx5 = intersect(idxMutants{5},idxSex{2});

l1 = length(idx1);
l2 = length(idx2);
l3 = length(idx3);
l4 = length(idx4);
l5 = length(idx5);

idxS = idxSex{2};


rmImg1 = mean(qRadon(idx1,:),1);
rmImg2 = mean(qRadon(idx2,:),1);
rmImg3 = mean(qRadon(idx3,:),1);
rmImg4 = mean(qRadon(idx4,:),1);
rmImg5 = mean(qRadon(idx5,:),1);


qRadon1 = qRadon(idx1,:) - rmImg1;
qRadon2 = qRadon(idx2,:) - rmImg2;
qRadon3 = qRadon(idx3,:) - rmImg3;
qRadon4 = qRadon(idx4,:) - rmImg4;
qRadon5 = qRadon(idx5,:) - rmImg5;

dImg1 = rmImg1 - rmImg3;
dImg2 = rmImg2 - rmImg3;
dImg4 = rmImg4 - rmImg3;
dImg5 = rmImg5 - rmImg3;


eImg1 = dImg1./sqrt(dot(dImg1,dImg1));
eImg2 = dImg2./sqrt(dot(dImg2,dImg2));
eImg4 = dImg4./sqrt(dot(dImg4,dImg4));
eImg5 = dImg5./sqrt(dot(dImg5,dImg5));


eImg = [eImg1; eImg2; eImg4; eImg5]';

Q = orth([rmImg1;rmImg2;rmImg3;rmImg4;rmImg5]');

Qr= orth([rmImg1-rmImg3;rmImg2-rmImg3;rmImg4-rmImg3;rmImg5-rmImg3]');

%imagesc(iradon(reshape(rmImg1,[size(R,1) size(R,2)]),0:179))
%imagesc(reshape(Q(:,2)',[size(R,1) size(R,2)]))

%sum((Q'*Qr).^2,2)
%sum(sum((Q'*Qr).^2,2).^2)/4


%imagesc(dotProdMat(eImg,eImg))
imagesc(abs(dotProdMat(eImg,eImg)))
colorbar
title('Pairwise scalar product of mutant directions')
%% PCA prjection

idx3 = idxMutants{3};

lIdx2 = length(idx2);


[a,b,c] = pca(qRadon);

[a3,b3,c3] = pca(qRadon(idx3,:));



%%

%imagesc(dotProdMat(a3(:,1:20),eImg))
imagesc(abs(dotProdMat(a3(:,1:10),eImg)))
colorbar
title('Pairwise scalar product of PCA and mutant directions')



%%

qRadon1Proj = dotProdMat(qRadon1', eImg);
qRadon2Proj = dotProdMat(qRadon2', eImg);
qRadon3Proj = dotProdMat(qRadon3', eImg);
qRadon4Proj = dotProdMat(qRadon4', eImg);
qRadon5Proj = dotProdMat(qRadon5', eImg);

imagesc(qRadon3Proj)

stdqRadon1Proj = std(qRadon1Proj);
stdqRadon2Proj = std(qRadon2Proj);
stdqRadon3Proj = std(qRadon3Proj);
stdqRadon4Proj = std(qRadon4Proj);
stdqRadon5Proj = std(qRadon5Proj);
%std(qRadon1Proj,0,1)


stdqRadon3 = std(qRadon3,0,1);
stdqRadon1 = std(qRadon1,0,1);
stdqRadon2 = std(qRadon2,0,1);
stdqRadon4 = std(qRadon4,0,1);
stdqRadon5 = std(qRadon5,0,1);


alpha3 = stdqRadon3Proj/sum(stdqRadon3);


fprintf('sigma/sigma * Dim = %03.3f\n',size(qRadon,2)*alpha3)

%% random planes 

tic;
tStart = toc;

nIteration = 1200;
projAmount = zeros(1,nIteration);

for iR = 1:nIteration 
    
    idxS = idxSex{2};
    
    rIdx3 = idx3;
    idxS = setdiff(idxS,rIdx3);
    
    rIdx1 = datasample(idxS,l1,'Replace',false);
    idxS = setdiff(idxS,rIdx1);
    
    rIdx2 = datasample(idxS,l2,'Replace',false);
    idxS = setdiff(idxS,rIdx2);
    
    %rIdx3 = datasample(idxS,l3,'Replace',false);
    %idxS = setdiff(idxS,rIdx3);
    
    rIdx4 = datasample(idxS,l4,'Replace',false);
    idxS = setdiff(idxS,rIdx4);
    
    rIdx5 = datasample(idxS,l5,'Replace',false);
    idxS = setdiff(idxS,rIdx5);
    
    
    rmImg1t = mean(qRadon(rIdx1,:),1);
    rmImg2t = mean(qRadon(rIdx2,:),1);
    rmImg3t = mean(qRadon(rIdx3,:),1);
    rmImg4t = mean(qRadon(rIdx4,:),1);
    rmImg5t = mean(qRadon(rIdx5,:),1);
    
    
    qRadon1t = qRadon(rIdx1,:) - rmImg1t;
    qRadon2t = qRadon(rIdx2,:) - rmImg2t;
    qRadon3t = qRadon(rIdx3,:) - rmImg3t;
    qRadon4t = qRadon(rIdx4,:) - rmImg4t;
    qRadon5t = qRadon(rIdx5,:) - rmImg5t;
    
    
    Qt = orth([rmImg1t; rmImg2t; rmImg3t; rmImg4t; rmImg5t]');
    
    Qrt= orth([rmImg1t-rmImg3t;rmImg2t-rmImg3t;...
        rmImg4t-rmImg3t;rmImg5t-rmImg3t]');
    
    scalProd = (Qrt'*Qr).^2;
    projAmount(iR) = sum(sum((Qrt'*Qr).^2,2).^2)/4;
    
end


fprintf('projection  = %0.2f ± %0.2f \n',...
    mean(projAmount),std(projAmount))

figure(1)
plot(sort(projAmount))
hold on
plot([1 nIteration],mean(projAmount)*[1 1])
hold off

figure(2)
histogram(projAmount)

tStop = toc;

fprintf('time to compute = %5d sec\n',round(tStop-tStart))

%% PCA


idx1 = union(idxMutants{1},idxMutants{3});
idx2 = union(idxMutants{2},idxMutants{3});
idx4 = union(idxMutants{4},idxMutants{3});
idx5 = union(idxMutants{5},idxMutants{3});

idx3 = idxMutants{3};

lIdx1 = length(idx1);
lIdx2 = length(idx2);
lIdx3 = length(idx3);
lIdx4 = length(idx4);
lIdx5 = length(idx5);

[a,b,c] = pca(qRadon);

[a1,b1,c1] = pca(qRadon(idx1,:));
[a2,b2,c2] = pca(qRadon(idx2,:));
[a3,b3,c3] = pca(qRadon(idx3,:));
[a4,b4,c4] = pca(qRadon(idx4,:));
[a5,b5,c5] = pca(qRadon(idx5,:));

%plot(cumsum(c)./sum(c),'o-')

f = figure('visible', 'off');


for k = 1:6
    %   figure
    imagesc(iradon(reshape(a(:,k),[size(R,1) size(R,2)]),0:179))
    colorbar
    set(gcf,'PaperOrientation','landscape');
    axis square
    title(sprintf('Eigenwing # %d',k))
    print(sprintf('plots/eig_%s_%s/diet_eigenface_%02d.pdf',...
        sex{jS},temp{kT},k),...
        '-dpdf','-bestfit')
    hold off
end

%% all

for j1 = 1:resolutioN
    for j2 = 1:resolutioN
        pImages(j1,j2,:) = double(histogram(squeeze(images(j1,j2,:)),15));
    end
end

normP = sum(pImages,3);

S = pImages.*log(pImages);
S(isnan(S)) = 0;

S = sum(S,3);
%we used nonNormalized prob distribution, so we have to take care of it.
S = -S./normP + log(normP);

imagesc(S)
colorbar
axis square
title('Entropy. All')


set(gcf,'PaperOrientation','landscape');
print('plots/entropy_bg_gr_all.pdf','-dpdf','-bestfit')
close all




%% 1


for jS = 1:2
    
    figure
    
    for j1 = 1:resolutioN
        for j2 = 1:resolutioN
            pImages(j1,j2,:) = ...
                double(histogram(squeeze(images(j1,j2,idxSex{jS})),15));
        end
    end
    
    normP = sum(pImages,3);
    
    S = -pImages.*log(pImages);
    S(isnan(S)) = 0;
    
    S = sum(S,3);
    % we used nonNormalized prob distribution,
    % so we have to take care of it.
    S = S./normP + log(normP);
    
    imagesc(S)
    colorbar
    axis square
    
    title(sprintf('Entropy. All %s',sex{jS}))
    
    
    set(gcf,'PaperOrientation','landscape');
    print(sprintf('plots/entropy_bg_gr_all_%s.pdf',sex{jS}),...
        '-dpdf','-bestfit')
    close all
    
end


for jC = 1:2
    
    figure
    
    for j1 = 1:resolutioN
        for j2 = 1:resolutioN
            pImages(j1,j2,:) = ...
                double(histogram(squeeze(images(j1,j2,idxChirality{jC})),15));
        end
    end
    
    normP = sum(pImages,3);
    
    S = -pImages.*log(pImages);
    S(isnan(S)) = 0;
    
    S = sum(S,3);
    % we used nonNormalized prob distribution,
    % so we have to take care of it.
    S = S./normP + log(normP);
    
    imagesc(S)
    colorbar
    axis square
    
    title(sprintf('Entropy. All %s',chirality{jC}))
    
    
    set(gcf,'PaperOrientation','landscape');
    print(sprintf('plots/entropy_bg_gr_all_%s.pdf',chirality{jC}),...
        '-dpdf','-bestfit')
    close all
    
end



for k = 1:length(mutants)
    
    figure
    
    for j1 = 1:resolutioN
        for j2 = 1:resolutioN
            pImages(j1,j2,:) = ...
                double(histogram(squeeze(images(j1,j2,idxMutants{k})),15));
        end
    end
    
    normP = sum(pImages,3);
    
    S = -pImages.*log(pImages);
    S(isnan(S)) = 0;
    
    S = sum(S,3);
    % we used nonNormalized prob distribution,
    % so we have to take care of it.
    S = S./normP + log(normP);
    
    imagesc(S)
    colorbar
    axis square
    
    title(sprintf('Entropy. All %s',mutants{k}))
    
    
    set(gcf,'PaperOrientation','landscape');
    print(sprintf('plots/entropy_bg_gr_all_%s.pdf',mutants{k}),...
        '-dpdf','-bestfit')
    close all
    
end


%% Kullback-Leibler divergence

%% 1
figure

for j1 = 1:resolutioN
    for j2 = 1:resolutioN
        pImages(j1,j2,:) = ...
            double(histogram(squeeze(images(j1,j2,idxSex{1})),15));
        qImages(j1,j2,:) = ...
            double(histogram(squeeze(images(j1,j2,idxSex{2})),15));
    end
end

normP = sum(pImages,3);
normQ = sum(qImages,3);

Dkl1 = -pImages.*log(qImages./pImages);
Dkl2 = -qImages.*log(pImages./qImages);

Dkl1(isnan(Dkl1)) = 0;
Dkl1(isinf(Dkl1)) = 0;

Dkl2(isnan(Dkl2)) = 0;
Dkl2(isinf(Dkl2)) = 0;

Dkl = Dkl1 + Dkl2;
%Dkl = Dkl2;

Dkl = sum(Dkl,3);
% we used nonNormalized prob distribution,
% so we have to take care of it.
Dkl = Dkl./normP - log(normP./normQ);

imagesc(Dkl)
colorbar
axis square

title(sprintf('D_{KL}. All %s vs %s',sex{1},sex{2}))

colorbar
set(gcf,'PaperOrientation','landscape');
print(sprintf('plots/dkl_bg_gr_all_%s_vs_%s.pdf',sex{1},sex{2}),...
    '-dpdf','-bestfit')
close all


figure

for j1 = 1:resolutioN
    for j2 = 1:resolutioN
        pImages(j1,j2,:) = ...
            double(histogram(squeeze(images(j1,j2,idxChirality{1})),15));
        qImages(j1,j2,:) = ...
            double(histogram(squeeze(images(j1,j2,idxChirality{2})),15));
    end
end

normP = sum(pImages,3);
normQ = sum(qImages,3);

Dkl1 = -pImages.*log(qImages./pImages);
Dkl2 = -qImages.*log(pImages./qImages);

Dkl1(isnan(Dkl1)) = 0;
Dkl1(isinf(Dkl1)) = 0;

Dkl2(isnan(Dkl2)) = 0;
Dkl2(isinf(Dkl2)) = 0;

Dkl = Dkl1 + Dkl2;
%Dkl = Dkl2;

Dkl = sum(Dkl,3);
% we used nonNormalized prob distribution,
% so we have to take care of it.
Dkl = Dkl./normP - log(normP./normQ);

imagesc(Dkl)
axis square

title(sprintf('D_{KL}. All %s vs %s',chirality{1},chirality{2}))

colorbar
set(gcf,'PaperOrientation','landscape');
print(sprintf('plots/dkl_bg_gr_all_%s_vs_%s.pdf',...
    chirality{1},chirality{2}), '-dpdf','-bestfit')
close all




figure

for j1 = 1:resolutioN
    for j2 = 1:resolutioN
        pImages(j1,j2,:) = ...
            double(histogram(squeeze(images(j1,j2,idxMutants{1})),15));
        qImages(j1,j2,:) = ...
            double(histogram(squeeze(images(j1,j2,idxMutants{2})),15));
    end
end

normP = sum(pImages,3);
normQ = sum(qImages,3);

Dkl1 = -pImages.*log(qImages./pImages);
Dkl2 = -qImages.*log(pImages./qImages);

Dkl1(isnan(Dkl1)) = 0;
Dkl1(isinf(Dkl1)) = 0;

Dkl2(isnan(Dkl2)) = 0;
Dkl2(isinf(Dkl2)) = 0;

Dkl = Dkl1 + Dkl2;
%Dkl = Dkl2;

Dkl = sum(Dkl,3);
% we used nonNormalized prob distribution,
% so we have to take care of it.
Dkl = Dkl./normP - log(normP./normQ);

imagesc(Dkl)
axis square

title(sprintf('D_{KL}. %s vs %s',mutants{1},mutants{2}))

colorbar
set(gcf,'PaperOrientation','landscape');
print(sprintf('plots/dkl_bg_gr_%s_vs_%s.pdf',mutants{1},mutants{2}),...
    '-dpdf','-bestfit')
close all


%% 2

for k = 1:length(mutants)
    for jS = 1:2
        %for jC = 1:2
        
        idxK1 = intersect(intersect(idxMutants{k},...
            idxChirality{1}), idxSex{jS});
        idxK2 = intersect(intersect(idxMutants{k},...
            idxChirality{2}), idxSex{jS});
        %idxK = intersect(idxMutants{k},idxSex{jS});
        
        figure
        
        for j1 = 1:resolutioN
            for j2 = 1:resolutioN
                pImages(j1,j2,:) = ...
                    double(histogram(squeeze(images(j1,j2,idxK1)),15));
                qImages(j1,j2,:) = ...
                    double(histogram(squeeze(images(j1,j2,idxK2)),15));
            end
        end
        
        normP = sum(pImages,3);
        normQ = sum(qImages,3);
        
        Dkl1 = -pImages.*log(qImages./pImages);
        Dkl2 = -qImages.*log(pImages./qImages);
        
        Dkl1(isnan(Dkl1)) = 0;
        Dkl1(isinf(Dkl1)) = 0;
        
        Dkl2(isnan(Dkl2)) = 0;
        Dkl2(isinf(Dkl2)) = 0;
        
        Dkl = Dkl1 + Dkl2;
        %Dkl = Dkl2;
        
        Dkl = sum(Dkl,3);
        % we used nonNormalized prob distribution,
        % so we have to take care of it.
        Dkl = Dkl./normP - log(normP./normQ);
        
        imagesc(Dkl)
        axis square
        
        title(sprintf('D_{KL}. %s %s || L vs R', mutants{k}, sex{jS}))
        
        colorbar
        set(gcf,'PaperOrientation','landscape');
        %end
        print(sprintf('plots/dkl_bg_gr_%s_%s_L_vs_R.pdf',...
            mutants{k}, sex{jS}), '-dpdf','-bestfit')
        close all
    end
end


for k = 1:length(mutants)
    for jC = 1:2
        %for jC = 1:2
        
        
        idxK1 = intersect(intersect(idxMutants{k},...
            idxChirality{jC}), idxSex{1});
        idxK2 = intersect(intersect(idxMutants{k},...
            idxChirality{jC}), idxSex{2});
        %idxK = intersect(idxMutants{k},idxSex{jS});
        
        figure
        
        for j1 = 1:resolutioN
            for j2 = 1:resolutioN
                pImages(j1,j2,:) = ...
                    double(histogram(squeeze(images(j1,j2,idxK1)),15));
                qImages(j1,j2,:) = ...
                    double(histogram(squeeze(images(j1,j2,idxK2)),15));
            end
        end
        
        normP = sum(pImages,3);
        normQ = sum(qImages,3);
        
        Dkl1 = -pImages.*log(qImages./pImages);
        Dkl2 = -qImages.*log(pImages./qImages);
        
        Dkl1(isnan(Dkl1)) = 0;
        Dkl1(isinf(Dkl1)) = 0;
        
        Dkl2(isnan(Dkl2)) = 0;
        Dkl2(isinf(Dkl2)) = 0;
        
        Dkl = Dkl1 + Dkl2;
        %Dkl = Dkl2;
        
        Dkl = sum(Dkl,3);
        % we used nonNormalized prob distribution,
        % so we have to take care of it.
        Dkl = Dkl./normP - log(normP./normQ);
        
        imagesc(Dkl)
        axis square
        
        title(sprintf('D_{KL}. %s %s || M vs F', mutants{k},...
            chirality{jC}))
        
        colorbar
        set(gcf,'PaperOrientation','landscape');
        %end
        print(sprintf('plots/dkl_bg_gr_%s_%s_M_vs_F.pdf',...
            mutants{k}, chirality{jC}), '-dpdf','-bestfit')
        close all
    end
end


%% Mutual information

figure

for j = 1:resolutioN^2
    
    [j1,j2] = ind2sub([resolutioN,resolutioN],j);
    pImg(j,:) = double(histogram(squeeze(images(j1,j2,:)),15));
end


j = 441;
for j1 = 90:110
    for j2 = 90:110
        
        j0 = sub2ind([resolutioN,resolutioN],j1,j2);
        fprintf('j = %d. j_0 = %d\n',j,j0)
        for k = resolutioN^2:-1:1
            
            %[j1,j2] = ind2sub([resolutioN,resolutioN],(j+j0));
            [k1,k2] = ind2sub([resolutioN,resolutioN],k);
            
            pqImg(j,k,:) = ...
                double(histogram(double(squeeze(images(j1,j2,:))).*...
                double(squeeze(images(k1,k2,:))),15));
        end
        
        for k =  resolutioN^2:-1:1
            
            miTest(j,k) = sum(squeeze(pqImg(j,k,:)).*...
                squeeze(log(squeeze(pqImg(j,k,:))./...
                (squeeze(pImg(j0,:).*pImg(k,:))'))));
            
        end
        j = j-1;
    end
end


normP = sum(pImages,3);

S = -pImages.*log(pImages);
S(isnan(S)) = 0;

S = sum(S,3);
% we used nonNormalized prob distribution,
% so we have to take care of it.
S = S./normP + log(normP);

imagesc(S)
colorbar
axis square

title(sprintf('Entropy. All %s',sex{jS}))


set(gcf,'PaperOrientation','landscape');
print(sprintf('plots/MI_bg_gr_all_%s.pdf',sex{jS}),...
    '-dpdf','-bestfit')
close all


%% 1
for jS = 1:2
    
    figure
    
    for j1 = 1:resolutioN
        for j2 = 1:resolutioN
            pImages(j1,j2,:) = ...
                double(histogram(squeeze(images(j1,j2,idxSex{jS})),15));
        end
    end
    
    normP = sum(pImages,3);
    
    S = -pImages.*log(pImages);
    S(isnan(S)) = 0;
    
    S = sum(S,3);
    % we used nonNormalized prob distribution,
    % so we have to take care of it.
    S = S./normP + log(normP);
    
    imagesc(S)
    colorbar
    axis square
    
    title(sprintf('Entropy. All %s',sex{jS}))
    
    
    set(gcf,'PaperOrientation','landscape');
    print(sprintf('plots/MI_bg_gr_all_%s.pdf',sex{jS}),...
        '-dpdf','-bestfit')
    close all
    
end

%%
%% Kolmogorov-Smirnov test


figure

ksSex = zeros(resolutioN,resolutioN);

for j1 = 1:resolutioN
    
    for j2 = 1:resolutioN
        
        hX = double(squeeze(images(j1,j2,idxSex{1})));
        gX = double(squeeze(images(j1,j2,idxSex{2})));
        
        [~,ksSex(j1,j2)] = kstest2(hX',gX');
    end
end


imagesc(ksSex)
colorbar
axis square

title(sprintf('Kolmogorov-Smirnov. %s vs %s. P-value',sex{1},sex{2}))
set(gcf,'PaperOrientation','landscape');

print(sprintf('plots/ks_all_sex.pdf'),...
    '-dpdf','-bestfit')


figure
ksMutants = zeros(resolutioN,resolutioN);


for jM1 = 1:length(mutants)
    for jM2 = jM1+1:length(mutants)
        
        for j1 = 1:resolutioN
            for j2 = 1:resolutioN
                
                hX = double(squeeze(images(j1,j2,idxMutants{jM1})));
                gX = double(squeeze(images(j1,j2,idxMutants{jM2})));
                
                [~,ksMutants(j1,j2)] = kstest2(hX',gX');
            end
        end
        
        
        imagesc(ksMutants)
        colorbar
        axis square
        
        title(sprintf('Kolmogorov-Smirnov. %s vs %s. P-value',...
            mutants{jM1},mutants{jM2}))
        set(gcf,'PaperOrientation','landscape');
        
        print(sprintf('plots/ks_all_%s_%s.pdf',mutants{jM1},...
            mutants{jM2}), '-dpdf','-bestfit')
    end
end
