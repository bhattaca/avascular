clear all;
close all;
baseOUTpath         = 'D:\Arindam\data\3x3\garb\\albili\\';

test();
%clusteringAlgo(baseOUTpath);


function clusteringAlgo(baseOUTpath)
data = csvread(fullfile([baseOUTpath,'dataCompare.csv']));
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
w = 1./var(datasmall);
[wcoeff,score,latent,tsquared,explained] = pca(datasmall,...
'VariableWeights',w);

boxplot(datasmall,'Orientation','horizontal')

c3 = wcoeff(:,1:3)

coefforth = diag(sqrt(w))*wcoeff;
I = coefforth'*coefforth;
I(1:3,1:3)

figure; gscatter( score(:,1),score(:,2), data(:,10),'kx');title ( 'Combined, 0-diseases'); 



%hold on;
%plot(score( (2451:1:2500),1), score( (2451:1:2500),2),'kx',...
%     'MarkerSize',15,'LineWidth',3) 


% GMModel = fitgmdist(score ( :,[1 2]),1);
% gmPDF = @(x,y)pdf(GMModel,[x y]);
 
%h = gscatter(score(:,1),score(:,2),data(:,10));
%hold on
%ezcontour(gmPDF,[-10 10],[-10 10]);
%ezcontour(@(x1,x2)pdf(GMModel,[x1 x2]),get(gca,{'XLim','YLim'}))
end

function test()
baseINpath          = 'D:\\Arindam\\data\\3x3\\albili\\';
baseOUTpath         = 'D:\Arindam\data\3x3\garb\\albili\\';
numFeaturesPerImage = 100;

% find all the normal images. 
filename  = 'D:\\Arindam\\data\\3x3\\albili\\density_common_resultsDisease.csv';

M         = readtable(filename);
% all 3x3 images
%temp3x3   = M ( M.ScanMmX==3,:);
%writetable(temp3x3,'density_common_results_3x3.csv','Delimiter',',');

[m,n]  = size ( M ); 

data = zeros (1, 10);

localStep = 1; % only for size 3x3; 
%iterate data
for k=1:m
    scanSize = M(k,2);
    if scanSize{1,1} == 3
        prefixArray = table2array(M(k,n-1));
        prefix      = prefixArray{1,1};
        imagepath   = fullfile ( [baseINpath, num2str(k),'\\',[prefix,'_SRLAngioEnface.png']]);
        %structural image
        imagepathStruc    = fullfile ( [baseINpath, num2str(k),'\\',[prefix,'_SRLStructuralEnface.png']]);
        currentimage      = imread ( imagepath);
        %currentimage     = imread ( 'D:\Arindam\data\3x3\albili\\22\\P430824AJ_Angiography 3x3 mm_11-7-2016_11-48-5_OS_sn2646_FlowCube_z.img_SRLAngioEnface.png');
        currentimageStruc     = imresize ( imread ( imagepathStruc),[1024,1024]);
        %currentimageMeanSub = currentimage - mean ( currentimage(:))*0.5;
        
        currentimage = im2double( currentimage);
        meanImg = mean2 ( currentimage)
        stdImg  = std2  ( currentimage)
        
        %currentimage2  = (currentimage - meanImg)/stdImg; 
        %currentimage3  = mat2gray ( currentimage2.*255);
        %figure; imshow ( currentimage3,[]);
        %figure; imshow ( currentimage,[]);
        
        garb = [k, " out of ", m];
        disp ( garb );
        %processimage ( i, localStep, k, baseOUTpath, numFeaturesPerImage,data);
        %
            %preprocessing
          
        I          = currentimage;
        background = imopen(I,strel('disk',55));
        I         = I - background; 
        Ihq        = (histeq(I)).*255;
        
        %anisotropic diffusion of Ihq
        num_iter = 35;
        delta_t = 1/7;
        kappa = 30;
        option = 2;
        %IhqAnisoSmooth = anisodiff2D(Ihq,num_iter,delta_t,kappa,option);
          
        %
        % Segmentation
        scaleStart =1; 
        scaleEnd   =10;
        Options.FrangiScaleRange = [scaleStart scaleEnd];
        Options.FrangiScaleRatio = 1;
        Options.BlackWhite = false;
        sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
        [outIm1, whatScale1,Direction1] = FrangiFilter2D(Ihq,Options);
        %figure; imshow (outIm1,[]);
        invI = (mat2gray(exp((outIm1).^2)));
        invI255 = invI.*255;
        
        num_iter = 35;
        delta_t = 1/7;
        kappa = 30;
        option = 2;
        ad = anisodiff2D(invI255,num_iter,delta_t,kappa,option);
        ad = mat2gray(ad)*255;
        %figure; imshow ( ad,[]);
        ad2 = uint8(ad);% figure; imshow (ad2,[]);
        %ad2log = logsample(ad2, 100, 500, 512, 512, 400, 360);
        %figure; imshow ( ad2log,[]); 
        %ad2logInv = logsampback(ad2log, 100, 500);
        %figure; imshow ( ad2logInv,[]);
        
        %break down each image to 600,2000
        imP = ImToPolar(ad2, 0.1, 1, 640, 1024);
        fun = @(block_struct) {
            imwrite(block_struct.data,[baseOUTpath, num2str(k),'-',num2str(rand()*1000),'.png']);};
        I2 = blockproc(uint8(Ihq),[128 128],fun);
        
              
        
        %imbin = imbinarize ( ad2);% figure; imshow ( imbin,[]);
        %skel1 =  bwmorph ( imbin, 'thin', 'inf');% figure; imshow ( skel1,[]);
        %skel1 = skel1.*255;
                
                
        %shock filter
        %outputImage = uint8( mat2gray(shock_filter(ad2,2,0.5,3)).*255);
        
         
        %imbin2 = imbinarize ( outputImage); figure; imshow ( imbin2,[]);
        

        
        %outputImage2 = shock_filter(skel1,2,15,15); figure; imshow  (histeq(mat2gray(outputImage2)),[]);
        
        %figure; imshow ( invI);
        
        % Smoothing parameters

    

        %imbin = imbinarize ( outIm1, 0.1);
        %figure; imshow (imbin ,[]);
        %imbin = bwareaopen  ( imbin,200 );
        %c = currentimage; 
        %currentimage = impyramid(currentimage, 'reduce');
        

        %D = mat2gray(bwdist(skel1,'quasi-euclidean')); figure; imshow ( (D),[]);
%        D = -D;
%D(~imbin) = Inf;
        %DBin = imbinarize ( uint8(imP), 0.1); 
        %figure; imshow ( DBin,[]);
        
        %bw = activecontour(outIm1, DBin, 500, 'Chan-Vese','SmoothFactor',2.5);
        %bw11 = activecontour(outIm1, DBin, 200, 'Chan-Vese','SmoothFactor',2.5);
        %bw22 = activecontour(outIm1, DBin, 600, 'Chan-Vese','SmoothFactor',2.5);
        %  imwrite (uint8(imresize (bw,[256,256])) ...
        %    ,[baseOUTpath,'bw',num2str(k),'.png']);
        %          imwrite (uint8(imresize (bw11,[256,256])) ...
        %    ,[baseOUTpath,'bw11',num2str(k),'.png']);
        %          imwrite (uint8(imresize (bw22,[256,256])) ...
        %    ,[baseOUTpath,'bw22',num2str(k),'.png'])
        %figure; imshow ([currentimage, uint8(bw).*255],[]);
        %{

        % coveredby large vessel
        scaleStart = 4; 
        scaleEnd   = 10;
        Options.FrangiScaleRange = [scaleStart scaleEnd];
        sigmas=Options.FrangiScaleRange(1):Options.FrangiScaleRatio:Options.FrangiScaleRange(2);
        [outIm1, whatScale1,Direction1] = FrangiFilter2D(imP,Options);
        outLargeNL = (mat2gray(exp((outLarge).^2))); % nonlinear 
        %figure; imshow ( outLargeNL,[]);
        outLargeBin        = (imbinarize( outLargeNL,0.1)); %figure; imshow ( outLargeBin,[]);
        outLargeBin        = bwareaopen(outLargeBin, 100, 4); %figure; imshow ( outLargeBin,[]);
        outLargeBin        = imdilate  ( outLargeBin, strel('disk',8)); %figure; imshow ( outLargeBin,[]);

        largevesselpath = [baseOUTpath, num2str(k),'img-largevessel-normal.png']; 
        imwrite ( uint8(outLargeBin.*255), largevesselpath);
        %}
        %binarize
        %imbin        = (imbinarize( invI,0.7));% figure; imshow ( [weightedImage, imbin],[]);
        %{ 
        imbin        = (imbinarize( Ihq,0.3));
         imbin        = imcomplement ( imbin);
         imbin        = imfill ( imbin, 'holes');
         imbin        = imopen(imbin, strel('disk',1)); %figure; imshow ( imbin,[]);

         %}
         %{
         background = imopen(currentimageStruc,strel('disk',55));
          strucimage         = currentimageStruc - background; 
         strucimage        = (histeq(strucimage));
         % figure; imshow ( imcomplement(strucimage));
         outFAZ = imbinarize ( imcomplement(strucimage),0.9);% figure; imshow ( outFAZ,[]); 
         outFAZmask = bwareaopen(outFAZ, 1000); %figure; imshow ( 1- outFAZmask,[]);
         
         %mix = imbin.*(1-outFAZ); figure; imshow ( mix,[]);
        %regionprops
        avascular_stats_per_img          = regionprops (bw,'centroid','Area','PixelIdxList','Eccentricity','EquivDiameter','Orientation' );

        %reorder by area
        [~,idx] = sort([avascular_stats_per_img.Area],'descend');
        sortd_avascular_stats_per_img=avascular_stats_per_img(idx);

        numFeaturesPerImageLocal = numFeaturesPerImage;

        [l,~]= size(sortd_avascular_stats_per_img);
        if l < numFeaturesPerImage
            %imshow ( currentimage );
            numFeaturesPerImageLocal = min ( numFeaturesPerImage , l);
        end

        tsortd_avascular_stats_per_img = sortd_avascular_stats_per_img((2 : 1 : numFeaturesPerImageLocal),:);
        if tsortd_avascular_stats_per_img(1,:).Area > 2000
            continue; 
        end

        %fill up data
        perimage = zeros(1024,1024);
        for kk = 1:numFeaturesPerImageLocal
            data(localStep,1)= k;
            data(localStep,2)= kk;
            data(localStep,3)= tsortd_avascular_stats_per_img(kk).Area;
            data(localStep,4)= abs(tsortd_avascular_stats_per_img(kk).Centroid(1)-512);
            data(localStep,5)= abs(tsortd_avascular_stats_per_img(kk).Centroid(2)-512);
            data(localStep,6)= tsortd_avascular_stats_per_img(kk).Eccentricity;
            data(localStep,7)= tsortd_avascular_stats_per_img(kk).EquivDiameter;
            data(localStep,8)= tsortd_avascular_stats_per_img(kk).Orientation;
            %percentage covered by large vessel.
            numLargeCovered = 0;
            for ll=1:length(tsortd_avascular_stats_per_img(kk).PixelIdxList)
                if outLargeBin(tsortd_avascular_stats_per_img(kk).PixelIdxList(ll)) == 1
                    numLargeCovered=numLargeCovered+1;
                end
            end

            data(localStep,9) = numLargeCovered;
            data(localStep,10)= 1;


            %
             redimage = zeros(1024,1024);
             redimage(tsortd_avascular_stats_per_img(kk).PixelIdxList)=128;
             perimage(tsortd_avascular_stats_per_img(kk).PixelIdxList)=255;
             outimage=zeros(1024,1024,3); 
             outimage(:,:,1)=uint8(currentimage).*0.5 + uint8(redimage); 
             outimage(:,:,2)=uint8(currentimage);
             outimage(:,:,3)=uint8(currentimage);

             imwrite (uint8(imresize (outimage,[256,256])) ...
            ,[baseOUTpath,num2str(k),'img',num2str(kk),'.png'])
            %
            
            %
        localStep = localStep+1;
        end 
        %per image make one image of all the chosen avascular regions 
        outimage2=zeros(1024,1024,3); 
        outimage2(:,:,1)=uint8(perimage); 
        outimage2(:,:,2)=uint8(currentimage);
        outimage2(:,:,3)=uint8(currentimage);
        imwrite ([uint8(cat(3,currentimage,currentimage,currentimage)), uint8(outimage2)] ...
            ,[baseOUTpath,'total-', num2str(k),'.png'])
        %}
    end
    
end

%csvwrite(fullfile([baseOUTpath,'dataCompare.csv']),data);

end
