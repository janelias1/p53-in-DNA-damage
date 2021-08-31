function dat=FFtoMatlab_importfilesol(fileToRead1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

fileID = fopen(fileToRead1,'r');

dataArray = textscan(fileID, '%[^\n\r]', ...
         'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
     
fclose(fileID);
        
dat = dataArray{:};
dat=dat(6:end-1);
dat=str2double(dat);
dat=dat';