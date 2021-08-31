% Use this function to plot 2D concentrations of p53, Mdm2 and Mdm2 mRNA

clear;

% Plot mesh
[points1, seg1, tri1]=FFtoMatlab_importfilemesh('Cell2Dnuc.msh');
[points2, seg2, tri2]=FFtoMatlab_importfilemesh('Cell2Dcyt.msh');

figure(1);
xlim1 = min(min(points1(1,:)),min(points2(1,:)));
xlim2 = max(max(points1(1,:)),max(points2(1,:)));
ylim1 = min(min(points1(2,:)),min(points2(2,:)));
ylim2 = max(max(points1(2,:)),max(points2(2,:)));
xlim([xlim1 xlim2]);
ylim([ylim1 ylim2]);
hold on;
pdemesh(points1,seg1,tri1);
pdemesh(points2,seg2,tri2);
hold off;
box on;

% path to the folder with .msh and .sol files, for example:
path = '/Users/janelias/FreeFem++/solfiles/';
av_files = dir(path);

for i=3:length(av_files)
    filename = av_files(i,1).name;
    if length(filename)==15 && strcmp(filename(1:4),'p53N')==1 ...
            && strcmp(filename(12:15),'.sol')==1
        
        p53nuc=FFtoMatlab_importfilesol(horzcat(path,filename));
        
        timestr=filename(7:11);
        timea=str2double(timestr)/100;
        
        filename2 = horzcat('p53C.',filename(6:end));
        p53cyt=FFtoMatlab_importfilesol(horzcat(path,filename2));
        
        figure(2); cla;
        hold on;
        pdeplot(points1,seg1,tri1,'zdata',p53nuc,'xydata',p53nuc,'mesh','off','colormap','jet','colorbar','off'); hold on;
        pdeplot(points2,seg2,tri2,'zdata',p53cyt,'xydata',p53cyt,'mesh','off','colormap','jet','colorbar','off');
        view(2);
        xlim([xlim1 xlim2]);
        ylim([ylim1 ylim2]);
        caxis([0 2]);
        title(horzcat('t = ',num2str(timea)),'FontSize',36 );
        hold off;
        box off; grid off;
        set(gca,'visible','off')
        drawnow;
        pause(0.2);
        
        if ismember(timea, [0,2,4,6,8,10,12,14,16,18])
            pause;
        end
                   
    end
end