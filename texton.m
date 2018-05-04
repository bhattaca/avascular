%read all the images 
close all; 
clear all;
baseINpath          = 'D:\\Arindam\\data\\3x3\\albili\\';
baseOUTpath         = 'D:\Arindam\data\3x3\garb\\albili\\';
outSummary = 'D:\Arindam\data\3x3\garb\albili\textons2\\';
%makeBlock (baseINpath, baseOUTpath);
%computeFeatures; 
vis (baseINpath, baseOUTpath, outSummary)

function vis( baseINpath, baseOUTpath, outSummary) 

%cluster 
data = csvread(fullfile([outSummary,'dataCompare.csv']));
indexMap = readtable(fullfile([outSummary,'index.csv']));
w = 1./var(data);
[wcoeff,score,latent,tsquared,explained] = pca(data,...
'VariableWeights',w);

boxplot(data,'Orientation','horizontal')

c3 = wcoeff(:,1:3)

coefforth = diag(sqrt(w))*wcoeff;
I = coefforth'*coefforth;
I(1:3,1:3)
%figure; plot ( score(:,1),score(:,2),'o' );

%opts = statset('Display','final');
%[idx,C] = kmeans(score,3,'Distance','cityblock',...
%    'Replicates',5,'Options',opts);

[idx,C,sumd,D] = kmeans(score (:, [1 2]),2,'MaxIter',10000,...
    'Display','final','Replicates',10);
X=score (:, [1 2]);
LL = table2array(indexMap(:,2));
figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off

clear selectedIndices;
delete ('D:\Arindam\data\3x3\garb\albili\SelectedTextures\*');
delete ('D:\Arindam\data\3x3\garb\albili\SelectedImages\*');
[l1, ~ ] = size ( indexMap ); 
distanceToTarget1 = 2.0;
distanceToTarget2 = 60;
loc = C(2,:);
loc = [-4.6, -0.9] 
p=1;
for ind=1:l1
    %if D(ind, 1) > distanceToTarget1 && D(ind, 2) > distanceToTarget2
    currentlocation = score (ind, [1 2]);
    l2dist  = sqrt ( (currentlocation(1)-loc(1))*(currentlocation(1)-loc(1)) ...
        + (currentlocation(2)-loc(2))*(currentlocation(2)-loc(2)));
    %l2dist
    if ( l2dist < distanceToTarget1)
        garb = [ind, " outof ", l1];
        disp ( garb );
        ss = indexMap(ind, 2); 
        cc=table2array(ss);
        prefix = cc{1,1};
        stringsplits = strsplit(prefix, {'-', '.'});
        %shift the texture 
        textureIn ='D:\Arindam\data\3x3\garb\albili\textons2\';
        texturePath = fullfile ( [textureIn, prefix]);
        textureOut  = 'D:\Arindam\data\3x3\garb\albili\SelectedTextures\';
        copyfile ( texturePath, textureOut);
        fldrName = stringsplits(2); 
        
        locX = stringsplits(5); locX = locX{1,1};locX = str2num(locX);
        locY = stringsplits(4); locY = locY{1,1};locY = str2num(locY);
        filename  = 'D:\\Arindam\\data\\3x3\\albili\\density_common_resultsDisease.csv';
        M         = readtable(filename);
        [m,n]  = size ( M );
        prefixArray = table2array(M(str2num(fldrName{1}),n-1));
        prefix      = prefixArray{1,1};
        imagepath   = fullfile ( [baseINpath, fldrName{1},'\\',prefix,'_SRLAngioEnface.png']);
        outpath     = fullfile ( ['D:\Arindam\data\3x3\garb\albili\SelectedImages','\\',prefix,'_SRLAngioEnface.png']);
        if exist(outpath, 'file') == 0
        copyfile ( imagepath, outpath );
        end
        lk     = imread ( imagepath );
        lk     = uint8(contrastEnhance ( lk ) ); 
        
        %for mm=1:128
        %    for nn=1:128
        %        lk(locY+mm,locX+nn)=255;
        %    end
        %end
        for offset = 1:128
            lk(locY,locX+offset)=255;
        end
        
        for offset = 1:128
            lk(locY+128,locX+offset)=255;
        end
        
        for offset = 1:128
            lk(locY+offset,locX)=255;
        end
        
        for offset = 1:128
            lk(locY+offset,locX+128)=255;
        end
        imwrite ( lk, outpath);
        
        
        
        selectedIndices{p}=ind; 
        p = p+1; 
    end
end

% show the selected points 
X=score (:, [1 2]);
figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
 plot(X(cell2mat(selectedIndices),1),X(cell2mat(selectedIndices),2),'g.','MarkerSize',12)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
%plot(X([3329:1:3456],1),X([3329:1:3456],2),'kx',...
%     'MarkerSize',15,'LineWidth',3) 
hold off
end

function computeFeatures
close all; 
clear all;
x = dir ('D:\Arindam\data\3x3\garb\albili\textons2\\*.png')

nfiles = size(x)

data = zeros ( 1,18);
for ii=1:nfiles
    f = fullfile (x(ii).folder,x(ii).name);
    currentimage = imread (f);
    %entropy 
    d1 = entropy ( currentimage/255);
    d2 = mean2   ( currentimage/255);
    d3 = std2    ( currentimage/255);
    currentimage = double(currentimage);
    %figure; imshow ( currentimage,[] ); 
    padcurrentimage = padarray(currentimage,[50 50],'both'); 
    
    num_iter = 5;
    delta_t = 1/7;
    kappa = 30;
    option = 2;
    ad = anisodiff2D(padcurrentimage,num_iter,delta_t,kappa,option);
    ad = uint8(mat2gray(ad)*255); %figure; imshow ( ad );
    
    
    [outIm, outScale, outDirection] = hess (1,4, double(ad));
    invI = (mat2gray(exp((outIm).^2)));
    invI255 = invI.*255;% figure; imshow ( invI255,[] );
    
    image = imcrop(invI255,[50 50 128 128]);
    outScale = imcrop(outScale,[50 50 128 128]);
    
    % scale information 
    d4 = histcounts ( image,10,'Normalization','probability');
    d5 = histcounts ( outScale,4,'Normalization','probability');
    d6 = Tamura_analysis(image);
    data(ii,:) = [d1, d2, d3, d4, d5, d6];
    M1{ii} = ii;
    M2{ii} = x(ii).name;
    csvwrite(fullfile([outSummary,'dataCompare.csv']),data);
end

fid = fopen( fullfile([outSummary,'index.csv']), 'w' );
M1temp = M1; 
M2temp = M2; 
for jj = 1 : length(M1 )
    l1 = M1temp {jj};
    l2 = M2temp {jj};
    fprintf( fid, '%d,%s\n', l1, l2 );
end
fclose( 'all' );


end

function [outIm1, whatScale1,Direction1]=hess ( scaleStart, scaleEnd, currentimage)

        Options.FrangiScaleRange = [scaleStart scaleEnd];
        Options.FrangiScaleRatio = 1;
        Options.BlackWhite = false;
        [outIm1, whatScale1,Direction1] = FrangiFilter2D(currentimage,Options);
end
function makeBlock(baseINpath, baseOUTpath)
    % find all the normal images. 
    filename  = 'D:\\Arindam\\data\\3x3\\albili\\density_common_resultsDisease.csv';
    M         = readtable(filename);
    [m,n]  = size ( M ); 
    data = zeros (1, 10);
    localStep = 1; % only for size 3x3; 
    for k=1:m
        scanSize = M(k,2);
        if scanSize{1,1} == 3
            prefixArray = table2array(M(k,n-1));
            prefix      = prefixArray{1,1};
            imagepath   = fullfile ( [baseINpath, num2str(k),'\\',[prefix,'_SRLAngioEnface.png']]);
            %structural image
            imagepathStruc    = fullfile ( [baseINpath, num2str(k),'\\',[prefix,'_SRLStructuralEnface.png']]);
            currentimage      = imread ( imagepath);
            %contrast enhance
            image     = contrastEnhance ( currentimage ) ; 
            %block apart
            blockImage ( uint8(image),128, baseOUTpath, k, prefix);
        end
    end
end

function [Ihq]= contrastEnhance ( currentimage)
        I          = currentimage;
        background = imopen(I,strel('disk',55));
        I          = I - background; 
        Ihq        = mat2gray(histeq(I)).*255;
end 

function blockImage ( currentimage, size, baseOUTpath,k, fname )
        fun = @(block_struct) imwrite(uint8(block_struct.data),...
            fullfile([baseOUTpath,'\\textons2\\row-', num2str(k),'-lc-',num2str(block_struct.location(1))...
            '-',num2str(block_struct.location(2)),'.png']));
        I2 = blockproc(currentimage,[size size],fun);
end 