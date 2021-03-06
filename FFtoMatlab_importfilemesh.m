function [points seg tri] = FFtoMatlab_importfilemesh(fileToRead1)
% Import the file
rawData1 = importdata(fileToRead1);

[unused,name] = fileparts(fileToRead1);
newData1.(genvarname(name)) = rawData1;

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

np=rawData1(1,1);
k=0;
for i=2:np+1
    k=k+1;
points(1,k)=rawData1(i,1);
points(2,k)=rawData1(i,2);
end

nt=rawData1(1,2);
k=0;

for i=np+2:2:np+1+2*nt
k=k+1;
tri(1,k)=rawData1(i,1);
tri(2,k)=rawData1(i,2);
tri(3,k)=rawData1(i,3);
tri(4,k)=rawData1(i+1,1);
end

k=0;
for i=1:nt
    k=k+1;
    lecseg(k,1)=tri(1,i);
    lecseg(k,2)=tri(2,i);
    k=k+1;
    lecseg(k,1)=tri(2,i);
    lecseg(k,2)=tri(3,i);
    k=k+1;
    lecseg(k,1)=tri(3,i);
    lecseg(k,2)=tri(1,i);
    
end

nlecseg=k;

k=0;
for i=1:nlecseg
    sw=0;
    for j=1:i-1
        if((lecseg(i,1)==lecseg(j,1) && lecseg(i,2)==lecseg(j,2)) || (lecseg(i,1)==lecseg(j,2) && lecseg(i,2)==lecseg(j,1)))
            sw=1;
        end
    end
    if(sw<1)
        k=k+1;
        seg(1,k)=lecseg(i,1);
        seg(2,k)=lecseg(i,2);
        seg(3,k)=0;
        seg(4,k)=1;
        seg(5,k)=k;
        seg(6,k)=1;
        seg(7,k)=0;

        
    end 
end

nseg=k;

