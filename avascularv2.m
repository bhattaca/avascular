clear all;
close all; 
% find all the normal images. 
filename  = 'D:\\Arindam\\data\\3x3\\albili\\density_common_results.csv';
M    = readtable(filename);
% all 3x3 images
temp3x3 = M ( M.ScanMmX==3,:)
writetable(temp3x3,'density_common_results_3x3.csv','Delimiter',',');


[m,n]  = size ( temp3x3); 
%iterate data
for k=1:m
    temp
end




pathToImages = 'P360613JB_Angiography 3x3 mm_5-4-2016_10-0-39_OS_sn2625_FlowCube_z.img'

directory       = 'D:\\Arindam\\data\\3x3\\albili\\';      % Full path of the directory to be searched in
filesAndFolders = dir(directory);     % Returns all the files and folders in the directory
filesInDir      = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory                    