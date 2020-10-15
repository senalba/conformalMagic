%% Init. We load images from an external drive.
workflow_0_dir_location


fileList = dir(strcat(dirInitiallyAligned,'*.tif'));
lFileList = length(fileList);


[populations, idxPop, lPopulation] = get_populations(fileList);

[idxL, idxR] = get_left_right(fileList);

nFields = length(strsplit(fileList(1).name,'_'));

labels = cell(lFileList,1);
sex = {'F','M'};
chirality = {'L','R'};

img = imread(strcat(dirAngularAligned,fileList(1).name));
indEmpty = find(img == 255);
R = radon(img);

imgRef = imread(strcat(dirInitiallyAligned,'G_outbred_24C_24_L.tif'));
imgRef(indEmpty) = 0;
histRef = imhist(imgRef);

images = uint8(zeros(size(img,1), size(img,2),lFileList));
pImages = zeros(size(img,1), size(img,2),15);
qImages = zeros(size(img,1), size(img,2),15);
%pImg = zeros(size(img,1)*size(img,2),15);
%pqImg = zeros(size(img,1)*size(img,2),size(img,1)*size(img,2),15);

%load('radon_distance_all_indices.mat')

for k = 1:lFileList
    
    name_str{k} = fileList(k).name;
    
    img = imread(strcat(dirInitiallyAligned,name_str{k}));
    %img = imread(strcat(dirInitiallyAligned,name_str{k}));
    img(indEmpty) = 0;
    %img = histeq(img, histRef);
    img = adapthisteq(img);
    %img  = imgradient(img);
    images(:,:,k)  = uint8(img);
end

%% RADON. we compute a radon image for each image


radonImages = zeros(size(R,1), size(R,2),lFileList);

for k = 1:lFileList
    
    radonImages(:,:,k) = radon(images(:,:,k));
end

R = radon(images(:,:,1));

qRadon = reshape(radonImages,[size(R,1)*size(R,2) lFileList])';
%qRadon1 = qRadon;
%qRadon = bsxfun(@rdivide,qRadon,sum(abs(qRadon),2));
%qRadon = bsxfun(@rdivide,qRadon,sum(qRadon,2));

[Z, zLoss] = tsne(qRadon,'Algorithm','exact','Distance','euclidean');




for indK = 1:length(populations)
    L(idxPop{indK}) = indK;
end

figure
%scatter(Z(:,1),Z(:,2),[],L,'filled')

for iP = 1:length(populations)
    hold on
    scatter(Z(idxPop{iP},1),Z(idxPop{iP},2),[],L(idxPop{iP}),'filled',...
        'DisplayName', populations{iP})
end
legend show
%legend(populatoins{1},populatoins{2},populatoins{3})

%% playing with centroids

%We choose female flies only.
%idxS = idxSex{1};

idx = cell(1,length(populations));

% wild type has idx=3
iW = 11;
%idxW = intersect(idxMutants{iW},idxSex{1});
idxW = idxPop{iW};

%we define set of wings for each mutants.
for iP = 1:length(populations)
    
    %idx{iM} = intersect(idxMutants{iM},idxS);
    idx{iP} = idxPop{iP};
    %lPopulation(iP) = length(idx{iP});
end

%we define vectors that points to the mutants of centroids
for iP = length(populations):-1:1
    
    [meanRadonImage{iP}, differenceRadonImage{iP},...
        directionMutant{iP}, unitDirectionMutant{iP}] = ...
        centroidsMutants(qRadon, idx{iP}, idxW);
end

%imagesc(iradon(reshape(meanRadonImage{iM},[size(R,1) size(R,2)]),0:179))

%for iM = [5,4,2,1]
%    figure
%    imagesc(iradon(reshape(meanRadonImage{iM},[size(R,1) size(R,2)]),0:179))
%    title(sprintf('mean %s wing',mutants{iM}))
%    axis square
%end
% we determine if directions are precise


nPerm = 20; %we perform 20 random labels assignments.

for iTr = nPerm:-1:1
    for iP = length(populations):-1:1
        
        idxP = randperm(lPopulation(iP));
        idxI{iTr,iP} = idx{iP}(idxP(1:floor(lPopulation(iP)/2)));
    end
end


%For each label assignments  we define directional vectors
for iTr = nPerm:-1:1
    for iP = length(populations)-1:-1:1
        [meanRadonImageI{iTr,iP}, differenceRadonImageI{iTr,iP},...
            directionMutantI{iTr,iP}, unitDirectionMutantI{iTr,iP}] = ...
            centroidsMutants(qRadon, idxI{iTr,iP}, idxI{iTr,3});
    end
end


%iM = 1;
%iTr = 1;
%dot(unitDirectionMutantI{iTr,iM},unitDirectionMutantI{10,iM})

%% Plots mean wings and directional vectors


for iP = 1:length(populations)

    figure
    %imagesc(iradon(reshape(meanRadonImage{iP},[size(R,1) size(R,2)]),0:179)) 
    %colorbar
    imshow(iradon(reshape(meanRadonImage{iP},[size(R,1) size(R,2)]),0:179),[])
    axis square
    %caxis([0 255])
    set(gca, 'XDir', 'reverse')
    set(gca, 'YDir', 'normal')
    title(sprintf('Mean wing for %s',populations{iP}))
    
    %savefig(sprintf('plots/mean_wing_img_%s.fig',populations{iP}))
    %print(sprintf('plots/mean_wing_img_%s.pdf',populations{iP}),'-dpdf','-bestfit')
end

for iP = length(populations)-1:-1:1

    figure
    imagesc(iradon(reshape(meanRadonImage{iP} - meanRadonImage{iW},...
        [size(R,1) size(R,2)]),0:179)) 
    %h = pcolor(iradon(reshape(meanRadonImage{iM} - meanRadonImage{3},...
    %    [size(R,1) size(R,2)]),0:179));
    %set(h, 'EdgeColor', 'none')
    set(gca, 'XDir', 'reverse')
    set(gca, 'YDir', 'normal')
    
    colorbar
    axis square
    caxis([-110 110])
    title(sprintf('Direction from Outbred to %s',populations{iP}))
    
    %savefig(sprintf('plots/inbred_direction_%s.fig',populations{iP}))
    %print(sprintf('plots/inbred_direction_%s.pdf',populations{iP}),'-dpdf','-bestfit')
end

close all

%% PCA prjection


[eigVectAll, bAll, eigValAll] = pca(qRadon);

[eigVectW,bW,eigValW] = pca(qRadon(idxW,:));


%eigVectW = eigVectW';

for iP = length(populations):-1:1
    
    [eigVect{iP},b{iP},eigVal{iP}] = pca(qRadon(idx{iP},:));
    %eigVect{iM} = eigVect{iM}';
end


%% Plot cov PCA

%load('data/eig_mutants_100_perm.mat','lM')

for iP = length(populations):-1:1
    
    hold off
    fIg = figure(iP);
    u = fIg.Units;
    fIg.Units = 'pixels';
    fIg.Position = [100 100 500 500];
    %figure('Units', 'pixels', ...
    %        'Position', [100 100 500 375]);
    hold on;
    %plot(cumsum(eigValW)/sum(eigValW),'.-')
    
    dataM{iP} = cumsum(eigVal{iP})/sum(eigVal{iP});
    
    wellSampled = line([1:lM(iP)], dataM{iP}(1:lM(iP)));
    
    notWellSampled = line([lM(iP) + 1 :length(eigVal{iP})],...
        dataM{iP}(lM(iP) + 1 :length(eigVal{iP})));
    
    set(wellSampled                         , ...
        'Color'             , [.7 .7 .7]    , ...
        'Marker'            , 'o'     , ...
        'MarkerSize'        , 6             , ...
        'LineStyle'         , '--'           , ...
        'LineWidth'         , 1.5             , ...
        'MarkerEdgeColor'   , [.2 .2 .2]    , ...
        'MarkerFaceColor'   , [.7 .7 .7]    );
    
    set(notWellSampled                                  , ...
        'Color'             , [0 0 .5]      , ...
        'LineStyle'         , '-'        , ...
        'Marker'            , '.'     , ...
        'MarkerSize'        , 4             , ...
        'LineWidth'         , 1.5             );
    
    ylim([0 1])
    xlim([1 length(eigVal{iP})])
    
    hTitle  = title (sprintf('%s mutant', mutants{iP})	);
    hXLabel = xlabel('Number of Eigenvectors'           );
    hYLabel = ylabel('Fraction of covariance '          );

    hLegend = legend( ...
                    [wellSampled,notWellSampled]    , ...
                    'Well-sampled eigenvectors'     , ...
                    'Non well-sampled eigenvectors'	);
    
          
    set(	gca                             , ...
            'FontName'   , 'Helvetica'      , ...
            'XScale','log'                  );
    set(    [hTitle, hXLabel, hYLabel]      , ...
            'FontName'   , 'Helvetica'     );
    set(    [hXLabel, hYLabel]              , ...
            'FontSize'   , 10               );
    set(    [hLegend, gca]                  , ...
            'FontSize'   , 18                );
    set(    hTitle                          , ...
            'FontSize'   , 24               ,  ...
            'FontWeight' , 'bold'           );
    ax = findobj(fIg,'type','axes', ...
            '-not','Tag','legend','-not','Tag','Colorbar');
	set(ax,'FontSize',10)     
    set(    [hXLabel, hYLabel]              , ...
            'FontSize'   , 18               );
    hLegend.Location = 'Best';
    
	set(    gca, ...
            'Box'         , 'off'       , ...
            'TickDir'     , 'out'       , ...
            'TickLength'  , [.02 .02]   , ...
            'XMinorTick'  , 'on'        , ...
            'YMinorTick'  , 'on'        , ...
            'YGrid'       , 'on'        , ...
            'XGrid'       , 'on'        , ...
            'XColor'      , [.3 .3 .3]  , ...
            'YColor'      , [.3 .3 .3]  , ...
            'YTick'       , 0:0.1:1    , ...
            'XTick'       , 0:50:length(eigVal{iP})    , ...
            'LineWidth'   , 1           );
   set(    hTitle                          , ...
            'FontSize'   , 24               ,  ...
            'FontWeight' , 'bold'           );        
    %print(fIg,sprintf('plots/cumulative_sum_%s.pdf',mutants{iM}),'-dpdf')
    %savefig(fIg,sprintf('plots/cumulative_sum_%s.fig',mutants{iM}))
end

%% Scalar Product for PC of the Wild type and mutants dirations

nPCAs = 20;
    
for iP = [5, 4, 2, 1]
    for iPCA = nPCAs:-1:1
        prodMutPcaW(iPCA,iP) = ...
            dot(unitDirectionMutant{iP},eigVectW(:,iPCA));
    end
end

figure
imagesc(abs(prodMutPcaW(:,[1,2,4,5])))
colorbar
colormap(esa)
caxis([0 1])
%set(gca, 'XTick', [1,2,4,5])
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})
xlabel('Mutants')
ylabel('Number of wild type PC')
set(gca, 'YTick', 1:nPCAs)

title('Scalar Product of the PC for Wild Type and directions to mutants')

%% scalar product random selection



for iTr = nPerm:-1:1
    for iP = [5, 4, 2, 1]
        for iPCA = nPCAs:-1:1
            prodMutPcaWI(iPCA,iP,iTr) = ...
                dot(unitDirectionMutantI{iTr,iP},eigVectW(:,iPCA));
        end
    end
end

figure
imagesc(abs(mean(prodMutPcaWI(:,[1,2,4,5]),3)))
colorbar
colormap(esa)
caxis([0 1])
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})
xlabel('Mutants')
ylabel('wild type PC')
set(gca, 'YTick', 1:nPCAs)
title('Abs Scalar Product of the random directional vector(random half) with PC')

% small divverence is  a good sign
figure
imagesc(abs(mean(prodMutPcaWI(:,[1,2,4,5],:),3)-...
    prodMutPcaW(:,[1,2,4,5])).^0.5)
colorbar
colormap(esa)
caxis([0 1])
set(gca, 'XTick', 1:5)
set(gca, 'YTick', 1:nPCAs)
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})
xlabel('Mutants')
ylabel('PC')
set(gca, 'YTick', 1:nPCAs)
title('Sqrt of the abs Difference between scalar products')
% small divverence is  a good sign

%% we would like to compute how consistent directions of the mutants are


for iTr = nPerm:-1:1
    for iP = [5, 4, 2, 1]
        
        mutProductI(iP,iTr) = ...
            dot(unitDirectionMutantI{iTr,iP},unitDirectionMutant{iP});
    end
end

figure
imagesc(mutProductI([1,2,4,5],:))
caxis([0 1])
colorbar
set(gca, 'YTick', 1:4)
yticks([1 2 3 4])
yticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})
ylabel('Mutants')
xlabel('Randomly selected directions')
title('Scalar product of a mutant directions with one for a random halves')

figure
bar([1:4],mean(mutProductI([1,2,4,5],:),2));
hold on
er = errorbar([1:4],mean(mutProductI([1,2,4,5],:),2),...
    std(mutProductI([1,2,4,5],:),0,2)/sqrt(nPerm));
er.Color = [0 0 0];
er.LineStyle = 'none';
ylim([0 1])
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})
title('Average scalar product of a mutant directions with one for a random halves')
hold off

%% we compute variability along the mutant directions



for iP = [5,4,2,1]
    
    %for l = numberMutants(3):-1:1
    for l = lPopulation(iP):-1:1
        projectionMut{iP}(l) = ...
            dot(differenceRadonImage{iP}(l,:),unitDirectionMutant{iP});
            %dot(qRadon(idx{iM}(l),:),unitDirectionMutant{iM});
            %dot(qRadon(idxW(l),:),unitDirectionMutant{iM});
        %dot(unitDirectionMutantI{iTr,iM},unitDirectionMutant{iM})
    end
    
    for l = lPopulation(3):-1:1
        projectionWT{iP}(l) = ...
            dot(differenceRadonImage{3}(l,:),unitDirectionMutant{iP});
            %dot(qRadon(idxW(l),:),unitDirectionMutant{iM});
    end
end


for iP = [5,4,2,1]
    
    fprintf('var for %s direction =  %d\n', ...
        mutants{iP},var(projectionMut{iP}))
end


for iP  = [5,4,2,1]
    fprintf('var for %s mutants  = %d\n',mutants{iP},...
        sum(var(differenceRadonImage{iP})))
end


for iP  = [5,4,2,1]
    
    rationM(iP) = var(projectionMut{iP})/sum(var(differenceRadonImage{iP}));
    %rationM(iM) = var(projectionMut{iM})/sum(var(differenceRadonImage{3}));
    fprintf('one dimensional var %s mutants/ (full var for %s) = %1.3f\n',...
        mutants{iP}, mutants{iP}, rationM(iP))
end

fprintf('\n')

for iP  = [5,4,2,1]
    
    ratioWT(iP) = var(projectionWT{iP})/sum(var(differenceRadonImage{3}));
    fprintf('one dimensional var WildType along %s/ (full var for WT) = %1.3f\n',...
        mutants{iP},  ratioWT(iP))
end
%% bootstrap for 1d Var(Mutant)/full Var(mutant) I


%unitDirectionMutantI{iTr,iM}



for iP = [5, 4, 2, 1]
    
    nM = length(idxI{1,iP});
    nW = length(idxI{1,3});
    
    for iTr = nPerm:-1:1
        for l = nM:-1:1
             
            projectionMutI(iTr,iP,l) = ...
                dot(differenceRadonImageI{iTr,iP}(l,:),...
                    unitDirectionMutantI{iTr,iP});
        end
        
        for k = nW:-1:1
             
            projectionMutWI(iTr,iP,k) = ...
                dot(differenceRadonImage{3}(k,:),...
                    unitDirectionMutantI{iTr,iP});
        end
    end   
end



%% bootstrap for 1d Var(Mutant)/full Var(mutant) II



for iP = [5, 4, 2, 1]

	rationM(iP) =   ...
        var(projectionMut{iP})/sum(var(differenceRadonImage{iP}));
    
    
    for iTr = nPerm:-1:1
        ratioMI(iP,iTr) = ...
            var(projectionMutI(iTr,iP,:),0,3)/ ...
                sum(var(differenceRadonImageI{iTr,iP}));
            
        ratioWI(iP,iTr) = ...
            var(projectionMutWI(iTr,iP,:),0,3)/ ...
                sum(var(differenceRadonImage{3}));    
    end
end


%%plot 
errorbar([1,2,3,4], mean(ratioMI([1,2,4,5],:),2),...
    std(ratioMI([1,2,4,5],:),0,2)./sqrt(nPerm),...
    'LineWidth', 2, 'DisplayName','Data')
hold on

plot([1,2,3,4], rationM([1,2,4,5]),'LineWidth',2,...
    'DisplayName','Data')
ylim([0 1])
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})
title('Bootstrap for mutants var ratios')


%% 1d var WT to 1st eaigval 




for iP = [5, 4, 2, 1]

	var2Eig(iP) = var(projectionWT{iP})/eigValW(1);
    
    for iTr = nPerm:-1:1
            
        var2EigI(iP,iTr) = ...
            var(projectionMutWI(iTr,iP,:),0,3)/eigValW(1);    
    end
end

%% bootstrap for variance

nBitr = 20;

for iTr = nBitr:-1:1
    
    idxS = idxSex{1};
    
    idxIR{iTr,3} = intersect(idxW, idxS);
    idxS = setdiff(idxS,idxIR{iTr,3});
    
    for iP = [5,4,2,1]
        
        idxIR{iTr,iP} = datasample(idxS,lPopulation(iP),'Replace',false);
        idxS = setdiff(idxS,idxIR{iTr,iP});
    end
    
    for iP = [5,4,2,1]
        [meanRadonImageR{iTr,iP}, differenceRadonImageR{iTr,iP},...
            directionMutantR{iTr,iP}, unitDirectionMutantR{iTr,iP}] = ...
            centroidsMutants(qRadon, idxIR{iTr,iP}, idxIR{iTr,3});
    end
    
    for iP = [5,4,2,1]
        %for l = numberMutants(3):-1:1
        for l = lPopulation(iP):-1:1
            projectionMutR{iTr,iP}(l) = ...
              dot(qRadon(idxIR{iTr,iP}(l),:),unitDirectionMutantR{iTr,iP});
                %dot(qRadon(idxW(l),:),unitDirectionMutantR{iTr,iM});
                
                
            %dot(unitDirectionMutantI{iTr,iM},unitDirectionMutant{iM})
        end
        
        for l = lPopulation(3):-1:1
            
            projectionWTR{iTr,iP}(l) = ...
                dot(qRadon(idxW(l),:),unitDirectionMutantR{iTr,iP});
                %dot(qRadon(idxW(l),:),unitDirectionMutantR{iTr,iM});
        end
        
        
        %fprintf('one dimensional var %s mutants/ (full var for %s) = %1.3f\n',...
        %    mutants{iM}, mutants{iM},...
        %    var(projectionMutR{iTr,iM})/sum(var(differenceRadonImageR{iTr,iM})))
        ratioMR(iTr,iP) =...
            var(projectionMutR{iTr,iP})/sum(var(differenceRadonImageR{iTr,iP}));
        ratioWR(iTr,iP) =...
            var(projectionWTR{iTr,iP})/sum(var(differenceRadonImage{3}));
       %var(projectionMutR{iTr,iM})/sum(var(differenceRadonImageR{iTr,iM}));
        var2EigR(iTr,iP) = var(projectionWTR{iTr,iP})/eigValW(1);
    end
    
    
end

%% plot Ratio of 1d Variance to Full Variance for each mutant
figure
btstrR = errorbar([1,2,3,4],mean(ratioMR([1,2,4,5],:),2),...
    std(ratioMR([1,2,4,5],:),0,2)./sqrt(nBitr),...
    'DisplayName','randomization');
hold on
%plot([1,2,3,4], rationM([1,2,4,5]),'LineWidth',2,...
%    'DisplayName','Data')

btstrI = errorbar([1,2,3,4], mean(ratioMI([1,2,4,5],:),2),...
    std(ratioMI([1,2,4,5],:),0,2)./sqrt(nPerm),...
    'LineWidth', 2, 'DisplayName','Data');
hold on
ylim([0 1])
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})


set(btstrI                        , ...
    'Color'             , [.7 .7 .7]    , ...
    'Marker'            , '.'     , ...
    'MarkerSize'        , 1             , ...
    'LineStyle'         , '-'           , ...
    'LineWidth'         , 1.5             , ...
    'MarkerEdgeColor'   , [.2 .2 .2]    , ...
    'MarkerFaceColor'   , [.7 .7 .7]    );

set(btstrR                                  , ...
    'Color'             , [0 0 .5]      , ...
    'LineStyle'         , '-.'        , ...
    'Marker'            , 'o'     , ...
    'MarkerSize'        , 2             , ...
    'LineWidth'         , 1             );


hTitle  = title('Ratio of 1d Variance to Full Variance for each mutant' );
hXLabel = xlabel('Mutants'                                              );
hYLabel = ylabel('Fraction of Variance '                                );

hLegend = legend(   [btstrI, btstrR]                , ...
                    'Ratio of 1d var to full var'	, ...
                    'Random label assignment'       );


set(	gca                             , ...
    'FontName'   , 'Helvetica'      , ...
    'XScale','log'                  );
set(    [hTitle, hXLabel, hYLabel]      , ...
    'FontName'   , 'Helvetica'     );
set(    [hXLabel, hYLabel]              , ...
    'FontSize'   , 10               );
set(    [hLegend, gca]                  , ...
    'FontSize'   , 18                );
set(    hTitle                          , ...
    'FontSize'   , 18               ,  ...
    'FontWeight' , 'bold'           );
ax = findobj(fIg,'type','axes', ...
    '-not','Tag','legend','-not','Tag','Colorbar');
set(ax,'FontSize',10)
set(    [hXLabel, hYLabel]              , ...
    'FontSize'   , 18               );
hLegend.Location = 'Best';

set(    gca, ...
    'Box'         , 'off'       , ...
    'TickDir'     , 'out'       , ...
    'TickLength'  , [.02 .02]   , ...
    'XMinorTick'  , 'on'        , ...
    'YMinorTick'  , 'on'        , ...
    'YGrid'       , 'on'        , ...
    'XGrid'       , 'on'        , ...
    'XColor'      , [.3 .3 .3]  , ...
    'YColor'      , [.3 .3 .3]  , ...
    'YTick'       , 0:0.2:1    , ...
    'LineWidth'   , 1           );
set(    hTitle                          , ...
        'FontSize'   , 18               ,  ...
        'FontWeight' , 'bold'           );
%set(gca, 'YScale','log')

%savefig(sprintf('plots/1d_var_2_full_var_each_mutant.fig'))
%print(sprintf('plots/1d_var_2_full_var_each_mutant.pdf'),'-dpdf')

%% Plot Ratio of 1d Variance for WT to Full Variance for WT


figure
btstrWR = errorbar([1,2,3,4],mean(ratioWR([1,2,4,5],:),2),...
    std(ratioWR([1,2,4,5],:),0,2)./sqrt(nBitr),...
    'DisplayName','randomization');
hold on
%plot([1,2,3,4], ratioWT([1,2,4,5]),'LineWidth',2,...
%    'DisplayName','Data')

btstrWI = errorbar([1,2,3,4], mean(ratioWI([1,2,4,5],:),2),...
    std(ratioWI([1,2,4,5],:),0,2)./sqrt(nPerm),...
    'LineWidth', 2, 'DisplayName','Data');

ylim([0 1])
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})




set(btstrWI                        , ...
    'Color'             , [.7 .7 .7]    , ...
    'Marker'            , '.'     , ...
    'MarkerSize'        , 1             , ...
    'LineStyle'         , '-'           , ...
    'LineWidth'         , 1.5             , ...
    'MarkerEdgeColor'   , [.2 .2 .2]    , ...
    'MarkerFaceColor'   , [.7 .7 .7]    );

set(btstrWR                                  , ...
    'Color'             , [0 0 .5]      , ...
    'LineStyle'         , '-.'        , ...
    'Marker'            , 'o'     , ...
    'MarkerSize'        , 2             , ...
    'LineWidth'         , 1             );


hTitle  = title('Ratio of 1d Variance for WT to Full Variance for WT');
hXLabel = xlabel('Mutants'                                              );
hYLabel = ylabel('Fraction of Variance'                                );

hLegend = legend(   [btstrWI, btstrWR]                , ...
                    'Ratio of 1d var to full var'	, ...
                    'Random label assignment'       );


set(	gca                             , ...
    'FontName'   , 'Helvetica'      , ...
    'XScale','log'                  );
set(    [hTitle, hXLabel, hYLabel]      , ...
    'FontName'   , 'Helvetica'     );
set(    [hXLabel, hYLabel]              , ...
    'FontSize'   , 10               );
set(    [hLegend, gca]                  , ...
    'FontSize'   , 18                );
set(    hTitle                          , ...
    'FontSize'   , 18               ,  ...
    'FontWeight' , 'bold'           );
ax = findobj(fIg,'type','axes', ...
    '-not','Tag','legend','-not','Tag','Colorbar');
set(ax,'FontSize',10)
set(    [hXLabel, hYLabel]              , ...
    'FontSize'   , 18               );
hLegend.Location = 'Best';

set(    gca, ...
    'Box'         , 'off'       , ...
    'TickDir'     , 'out'       , ...
    'TickLength'  , [.02 .02]   , ...
    'XMinorTick'  , 'on'        , ...
    'YMinorTick'  , 'on'        , ...
    'YGrid'       , 'on'        , ...
    'XGrid'       , 'on'        , ...
    'XColor'      , [.3 .3 .3]  , ...
    'YColor'      , [.3 .3 .3]  , ...
    'YTick'       , 0:0.2:1    , ...
    'LineWidth'   , 1           );
set(    hTitle                          , ...
        'FontSize'   , 18               ,  ...
        'FontWeight' , 'bold'           );
set(gca, 'YScale','log')

%savefig(sprintf('plots/1d_var_2_full_var_WT.fig'))
%print(sprintf('plots/1d_var_2_full_var_WT.pdf'),'-dpdf')

%%
figure

errorbar([1,2,3,4],mean(var2EigR([1,2,4,5],:),2),...
    std(var2EigR([1,2,4,5],:),0,2)./sqrt(nBitr),...
    'DisplayName','randomization')
hold on
plot([1,2,3,4], var2Eig([1,2,4,5]),'LineWidth',2,...
    'DisplayName','Data')
ylim([0 1])
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})
title('Bootstrap for WT Var to WT PCA Eigenvalue ratio')


%% Plot Ratio of 1d Variance for WT to first EigValue


fIg = figure

btstrW1eR = errorbar([1,2,3,4],mean(var2EigR([1,2,4,5],:),2),...
    std(var2EigR([1,2,4,5],:),0,2)./sqrt(nBitr),...
    'DisplayName','randomization');
hold on

btstrW1eI = errorbar( [1:4]                                       ,...
                    mean(var2EigI([1,2,4,5],:),2)               ,...
                    std(var2EigI([1,2,4,5],:),0,2)./sqrt(nPerm) ,...
                    'LineWidth',        2                       ,... 
                    'DisplayName',      'Data'                  );

ylim([0 1])
xticks([1 2 3 4])
xticklabels({mutants{1},mutants{2},mutants{4},mutants{5}})




set(btstrW1eI                           , ...
    'Color'             , [.7 .7 .7]	, ...
    'Marker'            , '.'           , ...
    'MarkerSize'        , 1             , ...
    'LineStyle'         , '-'           , ...
    'LineWidth'         , 1.5          	, ...
    'MarkerEdgeColor'   , [.2 .2 .2]    , ...
    'MarkerFaceColor'   , [.7 .7 .7]    );

set(btstrW1eR                             , ...
    'Color'             , [0 0 .5]      , ...
    'LineStyle'         , '-.'          , ...
    'Marker'            , 'o'           , ...
    'MarkerSize'        , 2             , ...
    'LineWidth'         , 1             );


hTitle  = title('Ratio of 1d Variance for WT to the 1st eigen value'    );
hXLabel = xlabel('Mutants'                                              );
hYLabel = ylabel('Fraction of Variance'                                 );

hLegend = legend(   [btstrW1eI, btstrW1eR]              , ...
                    'Ratio of 1d var to 1st eigen value', ...
                    'Random label assignment'           );


set(	gca                             , ...
    'FontName'   , 'Helvetica'          , ...
    'XScale','log'                      );

set(    [hTitle, hXLabel, hYLabel]      , ...
        'FontName'   , 'Helvetica'      );

set(    [hXLabel, hYLabel]              , ...
        'FontSize'   , 10               );

set(    [hLegend, gca]                  , ...
        'FontSize'   , 18              	);
    
set(    hTitle                          , ...
        'FontSize'   , 18               , ...
        'FontWeight' , 'bold'           );
    
ax = findobj(fIg,   'type'  , 'axes'    , '-not'    ,...
                    'Tag'   , 'legend'  , '-not'    , ...
                    'Tag'   , 'Colorbar'            );

set(ax,'FontSize',10)

set(    [hXLabel, hYLabel]              , ...
        'FontSize'   , 18               );
hLegend.Location = 'Best';

set(    gca, ...
        'Box'           , 'off'         , ...
        'TickDir'       , 'out'         , ...
        'TickLength'    , [.02 .02]     , ...
        'XMinorTick'    , 'on'          , ...
        'YMinorTick'    , 'on'          , ...
        'YGrid'         , 'on'          , ...
        'XGrid'         , 'on'          , ...
        'XColor'        , [.3 .3 .3]    , ...
        'YColor'        , [.3 .3 .3]    , ...
        'YTick'         , 0:0.2:1       , ...
        'LineWidth'     , 1             );
set(    hTitle                          , ...
        'FontSize'      , 18            ,  ...
        'FontWeight'	, 'bold'        );
set(gca,'YScale'        , 'log'         );

%savefig(sprintf('plots/1d_var_2_1_eigenvalue.fig'))
%print(sprintf('plots/1d_var_2_1_eigenvalue.pdf'),'-dpdf')


%% OLD

for iR = nPerm:-1:1
    
    idxS = idxSex{1};
    
    rIdx3 = intersect(idxW, idxS);
    idxS = setdiff(idxS,rIdx3);
    
    for iP = [1,2,4,5]
        
        rIdx{iR,iP} = datasample(idxS,lPopulation{iP},'Replace',false);
        idxS = setdiff(idxS,rIdx{iR,iP});
    end
    
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


%%

%imagesc(dotProdMat(a3(:,1:20),eImg))
for iP = [5, 4, 2, 1]
    
    eImg(iP,:) = unitDirectionMutant{iP};
end

imagesc(abs(dotProdMat(eigVectW(:,1:10),eImg([5,4,2,1],:)')))
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
    
    rIdx3 = idxW;
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



%% old
for iP = 1:5
    
    idx{iP} = intersect(idxMutants{iP},idxS);
    lPopulation(iP) = length(idx{iP});
    
    meanRadonImage{iP} = mean(qRadon(idx{iP},:),1);
    differenceRadonImage{iP} = qRadon(idx{iP},:) - meanRadonImage{iP};
end

%imagesc(iradon(reshape(meanRadonImage{iM},[size(R,1) size(R,2)]),0:179))
%imagesc(iradon(reshape(differenceRadonImage{iM}(1,:),[size(R,1) size(R,2)]),0:179))

%rmImg1 = mean(qRadon(idx1,:),1);
%qRadon1 = qRadon(idx1,:) - rmImg1;
%dImg1 = rmImg1 - rmImg3;
%eImg1 = dImg1./sqrt(dot(dImg1,dImg1));




%eImg = [eImg1; eImg2; eImg4; eImg5]';

%Q = orth([rmImg1;rmImg2;rmImg3;rmImg4;rmImg5]');
%Qr= orth([rmImg1-rmImg3;rmImg2-rmImg3;rmImg4-rmImg3;rmImg5-rmImg3]');

%imagesc(iradon(reshape(rmImg1,[size(R,1) size(R,2)]),0:179))
%imagesc(reshape(Q(:,2)',[size(R,1) size(R,2)]))

%sum((Q'*Qr).^2,2)
%sum(sum((Q'*Qr).^2,2).^2)/4


%imagesc(dotProdMat(eImg,eImg))
%imagesc(abs(dotProdMat(unitDirectionMutant,unitDirectionMutant)))
%colorbar
%title('Pairwise scalar product of mutant directions')
