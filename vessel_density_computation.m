%Algorithm for computing vessel density. 
% 5/3/2018
clear all; close all; 
%img = imread ( 'D:\Arindam\data\3x3\doheny\P5000-4343-DRR-101B_Angiography 3x3 mm_8-29-2016_10-22-12_OD_sn1380_FlowCube_z.png'); 
img = rgb2gray( imread ( 'C:\Users\u6abhatt\Downloads\garb.png'));
img = im2double ( img );
%figure; imshow ( img ) ; 

%pre-processing 
preprocessedImage = preProcess ( img ) ;% figure; imshow ( preprocessedImage,[]);
%Hessian Analysis 

scaleStart =4; 
scaleEnd   =4;
Options.FrangiScaleRange = [scaleStart scaleEnd];
Options.FrangiScaleRatio = 1;
Options.BlackWhite = false;
[outIm1, whatScale1,Direction1] = FrangiFilter2D(preprocessedImage,Options);
%figure; imshow (outIm1,[]);
outIm1 = mat2gray(outIm1)*255;
%invI = (mat2gray(exp((outIm1).^2)));
%binarize the hessian image
T = adaptthresh(outIm1, 0.4);
bin1 = imbinarize(outIm1,T);% figure; imshow ( bin1,[]);title ( 'binary');
%figure; imshow ( Direction1.*bin1,[]); title ( 'direction');
%masked direction
maskDir =  Direction1.*bin1;
figure; imshow ( maskDir,[]);title ( 'maskDir');

temp = zeros ( 1024, 1024);
temp = outIm1; 
nextrow = 440; 
nextcol = 321;
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
nextrow = 246;
nextcol = 601;
[temp,nextrow, nextcol]= metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
nextrow =769
nextcol =327
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
[temp,nextrow, nextcol] = metaTractsCombined (outIm1,temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );

%for k=1:15
%[outIm1] =  metaTractsNormal ( outIm1, 439,  475, maskDir(439,475) );
%[nextrow, nextcol] = metaTractsTangent ( outIm1, 439,  475,  maskDir, maskDir(439,475) );
%temp   =  metaTractsNormal ( temp, nextrow,  nextcol, maskDir(nextrow,nextcol) );
%[nextrow, nextcol, temp] = metaTractsTangent ( temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
%end


%temp = metaTractsCombined (temp, nextrow,  nextcol,  maskDir, maskDir(nextrow,nextcol) );
%[temprowNormal, tempcolNormal, temp]=metaTractsTangent( temp, nextrow, nextcol,  maskDir, maskDir(nextrow,nextcol));
%outIm1   =  metaTractsNormal ( outIm1, 440,  321, maskDir(440,321) );
%outIm1   =  metaTractsNormal ( outIm1, 580,  473, maskDir(580,473) );
figure; imshow( temp,[]);
% find subimage 

%gabor filter 
%{
wavelength = 4;
orientation = 90;
preprocessedImage = uint8(preprocessedImage);
[mag,phase] = imgaborfilt(preprocessedImage,wavelength,orientation);

figure
subplot(1,3,1);
imshow(preprocessedImage);
title('Original Image');
subplot(1,3,2);
imshow(mag,[])
title('Gabor magnitude');
subplot(1,3,3);
imshow(phase,[]);
title('Gabor phase');
%}

function [out]=metaTracts ( image, xind, yind, angle)

w = 1024; 
h = 1024;
out = image;

perp = 0; 
if angle > 90 
    perp = angle -90; 
else 
    perp = angle + 90; 
end
disp ( perp )
if perp <= 45
    xinc = tan ( degtorad(perp));
    yinc =1;
else 
    yinc = cot ( degtorad(perp));
    xinc =1;
end
for l =1:50
    
    x = xind +  xinc*l ; 
    y = yind +  yinc*l;
    %disp ( ['xinc', num2str(xinc),', x', num2str(x),',y ', num2str(y)]);
    x = round ( x ); 
    y = round ( y );
    %disp ( ['round x', num2str(x),',y ', num2str(y)]);
    out (x,y) = 255;
end
end

function [outimage, nextrow, nextcol]=metaTractsCombined   ( inimage, outimage, row, col , angleMatrix, angle) 
tangentUnitVectorX = cos(degtorad(angle));
tangentUnitVectorY = sin(degtorad(angle));

normalUnitVectorX = -1*tangentUnitVectorX; 
normalUnitVectorY = tangentUnitVectorY;

nextrow = row; 
nextcol = col;

currentMin = 180;
for l=-30:30
    temprow = row +  normalUnitVectorX*l; 
    tempcol = col -  normalUnitVectorY*l;   %% AB_DEBUG deal with matlab axis ( negating y ) 
    %disp ( ['xinc', num2str(xinc),', x', num2str(x),',y ', num2str(y)]);
    temprow = round ( temprow ); 
    tempcol = round ( tempcol );
    %disp ( ['round x', num2str(x),',y ', num2str(y)]);
    % do the tangent
    
    [temprowNormal, tempcolNormal, outimage]=metaTractsTangent( inimage, outimage,  temprow, tempcol, angleMatrix, angle);
    
    %see if this is the closest angle
    tempangle = angleMatrix ( temprowNormal, tempcolNormal);
    
    if tempangle <= currentMin && inimage ( temprow, tempcol) > 0 && temprowNormal ~=row && tempcolNormal ~=col
        currentMin = angle;
        nextrow = temprowNormal; 
        nextcol = tempcolNormal;
    end
    
    
    outimage (temprow,tempcol) = 255;
end

end
function [image]=metaTractsNormal( image, row, col, angle)

tangentUnitVectorX = cos(degtorad(angle));
tangentUnitVectorY = sin(degtorad(angle));

normalUnitVectorX = -1*tangentUnitVectorX; 
normalUnitVectorY = tangentUnitVectorY;
angle
for l=-20:20
    temprow = row +  normalUnitVectorX*l; 
    tempcol = col -  normalUnitVectorY*l;   %% AB_DEBUG deal with matlab axis ( negating y ) 
    %disp ( ['xinc', num2str(xinc),', x', num2str(x),',y ', num2str(y)]);
    temprow = round ( temprow ); 
    tempcol = round ( tempcol );
    %disp ( ['round x', num2str(x),',y ', num2str(y)]);
    image (temprow,tempcol) = 255;
end

end


function [nextrow, nextcol, outimage]=metaTractsTangent( inimage, outimage, row, col, angleMatrix, angle)

tangentUnitVectorX = cos(degtorad(angle));
tangentUnitVectorY = sin(degtorad(angle));

currentMin = 180; 
t = 0; 
%temp = zeros ( 1024, 1024);
for l=1:50
    temprow = row +  tangentUnitVectorY*l; 
    tempcol = col -  tangentUnitVectorX*l;   %% AB_DEBUG deal with matlab axis ( negating y ) 
    %disp ( ['xinc', num2str(xinc),', x', num2str(x),',y ', num2str(y)]);-a
    temprow = round ( temprow ); 
    tempcol = round ( tempcol );
    tempAngle = angleMatrix ( temprow, tempcol); 
    if abs(tempAngle - angle ) < currentMin && inimage ( temprow, tempcol) > 0 
        currentMin = abs(tempAngle - angle );
        t = l;
    end
    %disp ( ['round x', num2str(x),',y ', num2str(y)]);
    outimage (temprow,tempcol) = 255;
    %temp(temprow, tempcol) = 255; 
    %figure; imshow ( temp,[]);
end
nextrow = round (row +  tangentUnitVectorX*t); 
nextcol = round (col -  tangentUnitVectorY*t);
end



function [out]=preProcess(img)
%preprocess image
%Contrast-Limited Adaptive Histogram Equalization (CLAHE)
out = adapthisteq(img,'NumTiles',[16 16],'ClipLimit',0.005);
out = mat2gray ( out )*255;

end