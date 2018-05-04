%avascular detection

%read the images
close all; 
clear all;
%x = dir ('D:\Arindam\data\3x3\outputs\doheny-avascular\*.png') % this was
%the already extracted regions. 
%x = dir ('D:\Arindam\data\3x3\doheny\*.png'); %this is on the doheny data
x  = dir ( 'D:\\Arindam\\data\\3x3\\ZeissAvg\\normal\\*.png')
[nfiles,~] = size(x)
fileID = fopen('D:\\Arindam\\data\\3x3\\ZeissAvg\\expCompare.txt','w');


%%%%%%

numFeaturesPerImage = 100;
totalnumrows        = nfiles*numFeaturesPerImage;
garb = sprintf('totalnumrows %d',totalnumrows);
disp(garb)
data = zeros ( totalnumrows, 8);
%%%%%


r=1; 

globalStep = 1;
for ii=1:nfiles
    fprintf("image %d\n", ii);
    %for each image
    f            = fullfile (x(ii).folder,x(ii).name);
    fprintf(fileID,'%s maps to %d normal \n',f,ii);
    currentimage = imread (f);
    %preprocessing
    I = im2double(currentimage); 
    background = imopen(I,strel('disk',55));
    I2 = I - background; 
    Ihq = (histeq(I2));
    %
    % Segmentation
    scaleStart =1; 
    scaleEnd   =10;
    Options.FrangiScaleRange = [scaleStart scaleEnd];
    Options.FrangiScaleRatio = 1;
    Options.BlackWhite = false;
    sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
    [outIm1, whatScale1,Direction1] = Hessian_Vesselness(Ihq,Options,sigmas);
    invI = (mat2gray(exp((1-outIm1).^2)));
    %figure; imshow ( invI); title ( 'invI');
    %figure; imshow ( I,[]); title ( 'I' );
    %figure; imshow ( 1.0-invI);
    %figure; imshow (I.*( 1.0-invI),[]);
    weightedImage = I.*( 1.0-invI);
    %
    
    %coveredby large vessel
    scaleStart =4; 
    scaleEnd   =10;
    Options.FrangiScaleRange = [scaleStart scaleEnd];
    sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
    [outLarge, garb1,garb1] = Hessian_Vesselness(I,Options,sigmas);
    outLargeNL = (mat2gray(exp((outLarge).^2))); % nonlinear 
    %figure; imshow ( outLargeNL,[]);
    outLargeBin        = (imbinarize( outLargeNL,0.1)); %figure; imshow ( outLargeBin,[]);
    outLargeBin        = bwareaopen(outLargeBin, 100, 4); %figure; imshow ( outLargeBin,[]);
    outLargeBin        = imdilate ( outLargeBin, strel('disk',8)); %figure; imshow ( outLargeBin,[]);
    
    imwrite ( uint8(outLargeBin.*255),['D:\\Arindam\\data\\3x3\\ZeissAvg\\normalLabel\\',num2str(ii),'img-largevessel-normal.png']);
    
    %binarize
    imbin        = (imbinarize( invI,0.7));% figure; imshow ( [weightedImage, imbin],[]);
    %imbin        = imfill(imbin,'holes'); % figure; imshow ( imbin);
    %regionprops
    avascular_stats_per_img          = regionprops(imbin,'centroid','Area','PixelIdxList','Eccentricity','EquivDiameter','Orientation' );
    [m,n]         = size(avascular_stats_per_img);
    %reorder by area
    [~,idx]=sort([avascular_stats_per_img.Area],'descend');
    sortd_avascular_stats_per_img=avascular_stats_per_img(idx);
    
    numFeaturesPerImageLocal = numFeaturesPerImage;
    
    [l,p]= size(sortd_avascular_stats_per_img);
    if l <numFeaturesPerImage
        %imshow ( currentimage );
        numFeaturesPerImageLocal = min ( numFeaturesPerImage , l)
    end
    
    %top numFeaturesPerImage
    %tsortd_avascular_stats_per_img = sortd_avascular_stats_per_img((1:1:numFeaturesPerImageLocal),:);
    
    tsortd_avascular_stats_per_img = sortd_avascular_stats_per_img((2:1:numFeaturesPerImageLocal+1),:);
    %fill up data
    for kk =1:numFeaturesPerImageLocal
        data(globalStep,1)= ii;
        data(globalStep,2)= kk;
        data(globalStep,3)= tsortd_avascular_stats_per_img(kk).Area;
        data(globalStep,4)= abs(tsortd_avascular_stats_per_img(kk).Centroid(1)-512);
        data(globalStep,5)= abs(tsortd_avascular_stats_per_img(kk).Centroid(2)-512);
        data(globalStep,6)= tsortd_avascular_stats_per_img(kk).Eccentricity*100;
        data(globalStep,7)= tsortd_avascular_stats_per_img(kk).EquivDiameter;
        data(globalStep,8)= tsortd_avascular_stats_per_img(kk).Orientation;
        %percentage covered by large vessel.
        numLargeCovered =0;
        for ll=1:length(tsortd_avascular_stats_per_img(kk).PixelIdxList)
            if outLargeBin(tsortd_avascular_stats_per_img(kk).PixelIdxList(ll)) == 1
                numLargeCovered=numLargeCovered+1;
            end
        end
        
        data(globalStep,9)= numLargeCovered;
        data(globalStep,10)= 1;
        
        
        %
         redimage = zeros(1024,1024);
         redimage(tsortd_avascular_stats_per_img(kk).PixelIdxList)=128;
         outimage=zeros(1024,1024,3); 
         outimage(:,:,1)=uint8(currentimage).*0.5 + uint8(redimage); 
         outimage(:,:,2)=uint8(currentimage);
         outimage(:,:,3)=uint8(currentimage);

        [filepath,name,ext] = fileparts(x(ii).name);
         imwrite (uint8(imresize (cat(3,redimage,currentimage,currentimage),[128,128])) ...
        ,['D:\\Arindam\\data\\3x3\\ZeissAvg\\normalLabel\\',num2str(ii),'img',num2str(kk),'normal.png'])
        %
        
        
        
        globalStep= globalStep+1;
    end
   


 %   for  index=1:numFeaturesPerImageLocal
 
  %  end
end

%  disease data 
x  = dir ( 'D:\\Arindam\\data\\3x3\\ZeissAvg\\disease\\*.png')
[nfiles,~] = size(x)
offset = globalStep;
globalStep=1;
for ii=1:nfiles
    fprintf("image %d\n", ii);
    %for each image
    f            = fullfile (x(ii).folder,x(ii).name);
    fprintf(fileID,'%s maps to %d disease\n',f,ii);
    currentimage = imread (f);
    %preprocessing
    I = im2double(currentimage); 
    background = imopen(I,strel('disk',55));
    I2 = I - background; 
    Ihq = (histeq(I2));
    %
    % Segmentation
    scaleStart =1; 
    scaleEnd   =10;
    Options.FrangiScaleRange = [scaleStart scaleEnd];
    Options.FrangiScaleRatio = 1;
    Options.BlackWhite = false;
    sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
    [outIm1, whatScale1,Direction1] = Hessian_Vesselness(Ihq,Options,sigmas);
    invI = (mat2gray(exp((1-outIm1).^2)));
    %figure; imshow ( invI);
    %figure; imshow ( 1.0-invI);
    %figure; imshow (I.*( 1.0-invI),[]);
    weightedImage = I.*( 1.0-invI);
    %
        %coveredby large vessel
    scaleStart =4; 
    scaleEnd   =10;
    Options.FrangiScaleRange = [scaleStart scaleEnd];
    sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
    [outLarge, garb1,garb1] = Hessian_Vesselness(I,Options,sigmas);
    outLargeNL = (mat2gray(exp((outLarge).^2))); % nonlinear 
    %figure; imshow ( outLargeNL,[]);
    outLargeBin        = (imbinarize( outLargeNL,0.1)); %figure; imshow ( outLargeBin,[]);
    outLargeBin        = bwareaopen(outLargeBin, 100, 4); %figure; imshow ( outLargeBin,[]);
    outLargeBin        = imdilate ( outLargeBin, strel('disk',8)); %figure; imshow ( outLargeBin,[]);
    imwrite ( uint8(outLargeBin.*255),['D:\\Arindam\\data\\3x3\\ZeissAvg\\diseaseLabel\\',num2str(ii),'img-largevessel-disease.png']);
    
    %binarize
    %imbin        = imcomplement(imbinarize( weightedImage,0.1));% figure; imshow ( [weightedImage, imbin],[]);
    %imbin        = (imbinarize( invI,0.7));% figure; imshow ( [weightedImage, imbin],[]);
    imbin        = (imbinarize( invI,0.7));
    %regionprops
    avascular_stats_per_img          = regionprops(imbin,'centroid','Area','PixelIdxList','Eccentricity','EquivDiameter','Orientation' );
    [m,n]         = size(avascular_stats_per_img);
    %reorder by area
    [~,idx]=sort([avascular_stats_per_img.Area],'descend');
    sortd_avascular_stats_per_img=avascular_stats_per_img(idx);
    
    numFeaturesPerImageLocal = numFeaturesPerImage;
    
    [l,p]= size(sortd_avascular_stats_per_img);
    if l <numFeaturesPerImage
        %imshow ( currentimage );
        numFeaturesPerImageLocal = min ( numFeaturesPerImage , l)
    end
    
    %top numFeaturesPerImage
    %tsortd_avascular_stats_per_img = sortd_avascular_stats_per_img((1:1:numFeaturesPerImageLocal),:);
    
    tsortd_avascular_stats_per_img = sortd_avascular_stats_per_img((2:1:numFeaturesPerImageLocal+1),:);
    %fill up data
    for kk =1:numFeaturesPerImageLocal
        data(offset + globalStep,1)= ii;
        data(offset + globalStep,2)= kk;
        data(offset + globalStep,3)= tsortd_avascular_stats_per_img(kk).Area;
        data(offset + globalStep,4)= abs(tsortd_avascular_stats_per_img(kk).Centroid(1)-512);
        data(offset + globalStep,5)= abs(tsortd_avascular_stats_per_img(kk).Centroid(2)-512);
        data(offset + globalStep,6)= tsortd_avascular_stats_per_img(kk).Eccentricity*100;
        data(offset + globalStep,7)= tsortd_avascular_stats_per_img(kk).EquivDiameter;
        data(offset + globalStep,8)= tsortd_avascular_stats_per_img(kk).Orientation;
        
                %percentage covered by large vessel.
        numLargeCovered =0;
        for ll=1:length(tsortd_avascular_stats_per_img(kk).PixelIdxList)
            if outLargeBin(tsortd_avascular_stats_per_img(kk).PixelIdxList(ll)) == 1
                numLargeCovered=numLargeCovered+1;
            end
        end
        
        data(offset + globalStep,9)= numLargeCovered;
        data(offset + globalStep,10)= 0; %disease
        
        
        %
         redimage = zeros(1024,1024);
         redimage(tsortd_avascular_stats_per_img(kk).PixelIdxList)=128;
         outimage=zeros(1024,1024,3); 
         outimage(:,:,1)=uint8(currentimage).*0.5+  uint8(redimage); 
         outimage(:,:,2)=uint8(currentimage);
         outimage(:,:,3)=uint8(currentimage);

        [filepath,name,ext] = fileparts(x(ii).name);
         imwrite (uint8(imresize (cat(3,redimage,currentimage,currentimage),[128,128])) ...
        ,['D:\\Arindam\\data\\3x3\\ZeissAvg\\diseaseLabel\\',num2str(ii),'img',num2str(kk),'disease.png'])
        %
        
        globalStep= globalStep+1;
    end
end
fclose(fileID);


datasmall = data ( :,[3,4,5,6,7,8,9]);
%data standardization 
[m,numFeatures] = size ( datasmall );
for l = 1:numFeatures
    temp = datasmall ( :,l);
    %figure; hist ( temp ); title ( 'old');
    meanTemp = mean ( temp ); 
    stdTemp  = std  ( temp );
    temp  = temp-meanTemp; 
    temp  = temp/stdTemp;
    datasmall(:,l) = temp;
    garb = ['l ',num2str(l), ' mean ', num2str(meanTemp), ' std ', num2str(stdTemp)];
    disp ( garb );
    %figure; hist ( temp ); title ( 'new');
end


%[W,H] = nnmf(datasmall,2);
%biplot(H','scores',W,'varlabels',{'1','2','3','4','5','6'});

csvwrite('D:\\Arindam\\data\\3x3\\ZeissAvg\\dataCompare.csv',data);
data = csvread ( 'D:\\Arindam\\data\\3x3\\ZeissAvg\\dataCompare.csv' );


%Save the tsne results 
Y = tsne (datasmall); 
csvwrite('D:\\Arindam\\data\\3x3\\ZeissAvg\\tsne.csv',Y);
%sanity check save the datasmall
csvwrite('D:\\Arindam\\data\\3x3\\ZeissAvg\\datasmall.csv',datasmall);


%clustering 
figure; gscatter( Y(:,1), Y(:,2), data(:,10));title ( 'Combined, 0-diseases'); hold on; 
plot(Y( (4601:1:4700),1),Y( (4601:1:4700),2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
 
 gscatter( Y( (13502:1:13601),1), Y((13502:1:13601) ,2), data((13502:1:13601),1),'Black','+',12);
 
figure; gscatter(Y(1:1:9000,1),Y(1:1:9000,2),data(1:1:9000,10), 'cyan'); title ( 'Normal');
figure; gscatter(Y(9001:1:23801,1),Y(9001:1:23801,2),data(9001:1:23801,10), 'red'); title ( 'Disease');

figure; gscatter(Y(9001:1:23801,1),Y(9001:1:23801,2),data(9001:1:23801,10), 'red'); title ( 'Disease');
%kmeans
Ynormal = Y(1:1:9000,:);
[idx,C] = kmeans(Ynormal,30,'Replicates',5);

cmap = jet(30);
figure;
plot(Ynormal(idx==1,1),Ynormal(idx==1,2),'.','Color', cmap(1, :),'MarkerSize',12)
hold on
for ii =2:30
    plot(Ynormal(idx==ii,1),Ynormal(idx==ii,2),'.','Color', cmap(ii, :),'MarkerSize',12)
end
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
title 'Normal Cluster Assignments and Centroids'
hold off
%kmeans disease
Ydisease = Y(9001:1:23801,:);
[idx,C] = kmeans(Ydisease,30,'Replicates',5);

cmap = jet(30);
figure;
plot(Ydisease(idx==1,1),Ydisease(idx==1,2),'.','Color', cmap(1, :),'MarkerSize',12)
hold on
for ii =2:30
    plot(Ydisease(idx==ii,1),Ydisease(idx==ii,2),'.','Color', cmap(ii, :),'MarkerSize',12)
end
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
title 'Disease Cluster Assignments and Centroids'
hold off




%
% read the disease image 
%{
x = dir ('D:\Arindam\data\3x3\ZeissAvg\normalBMP\*.bmp');
[nfiles,~] = size(x)
r=1;
for ii=1:nfiles
    f            = fullfile (x(ii).folder,x(ii).name);
    currentimage = imread (f);
    imwrite ( rgb2gray ( currentimage), ['D:\Arindam\data\3x3\ZeissAvg\normal\', num2str( r ) , '.png'])
    r = r+1;
end
%}

Y = tsne(datasmall,'Algorithm','exact','Distance','mahalanobis');
subplot(2,2,1)
gscatter(Y(:,1),Y(:,2),data(:,8))
title('Mahalanobis')

Y = tsne(datasmall,'Algorithm','exact','Distance','cosine');
subplot(2,2,2)
gscatter(Y(:,1),Y(:,2),data(:,8))
title('Cosine')


Y = tsne(datasmall,'Algorithm','exact','Distance','chebychev');
subplot(2,2,3)
gscatter(Y(:,1),Y(:,2),data(:,8))
title('Chebychev')






rng('default') % for fair comparison
Y = tsne(datasmall,'Algorithm','exact','Distance','euclidean');
subplot(2,2,4)
gscatter(Y(:,1),Y(:,2),data(:,8))
title('Euclidean')

%PCA
w = 1./var(datasmall);
[wcoeff,score,latent,tsquared,explained] = pca(datasmall,...
'VariableWeights',w);

boxplot(datasmall,'Orientation','horizontal')

c3 = wcoeff(:,1:3)

coefforth = diag(sqrt(w))*wcoeff;
I = coefforth'*coefforth;
I(1:3,1:3)

figure()
plot(score(:,1),score(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

figure; gscatter( score(:,1),score(:,2), data(:,10));title ( 'Combined, 0-diseases'); hold on; 
plot(score( (4601:1:4700),1),score( (4601:1:4700),2),'kx',...
     'MarkerSize',10,'LineWidth',3) ; hold on; 
 plot(score( (14502:1:14601),1),score( (14502:1:14601),2),'gx',...
     'MarkerSize',10,'LineWidth',3) 

 GMModel = fitgmdist(score,3);
 gmPDF = @(x,y)pdf(GMModel,[x y]);
 
h = gscatter(score(:,1),score(:,2),data(:,10));
hold on
ezcontour(gmPDF,[-10 10],[-10 10]);
ezcontour(@(x1,x2)pdf(GMModel,[x1 x2]),get(gca,{'XLim','YLim'}))
 