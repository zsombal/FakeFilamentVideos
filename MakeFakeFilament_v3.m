clear all;
close all;
tic

speed_mean = 5
speed_std = 2
fmot = 1; % between [0, 1.0], at 0.5, 50% of filaments starts moving
fmot_change = 0; % between [-0.5, 0.5]. E.g. for 0.5 all filaments will move by the end
fmot_change_rate = 0; % between [1,100]. 1 changes fastest, 5 changes slower, and so on
random = 0; % random motion of the non-moving filaments
speckle = 0.01; % speckle noise on the hole image

height = 480;
width = 720;
outside = 400;

bend=0.5;
filamentcount=30;
speed = abs(normrnd(speed_mean, speed_std,[1, filamentcount]))
% speed=[1.5, 10.5]
thickness=1;
minimumlength=10;
maximumlength=100;

output = ['speed', num2str(speed_mean), '_fmot',num2str(fmot),...
    '_fmotChange',num2str(fmot_change), '_fmotChangeRate',num2str(fmot_change_rate),...
     '_RandomMot',num2str(random),'_Speckle',num2str(speckle),'.avi']

aviobj = VideoWriter(output);
open(aviobj)


for i=1:filamentcount
    xstart(i)=round(rand*width)+outside/2;
    ystart(i)=round(rand*height)+outside/2;
    xylength(i)=round(rand*(maximumlength-minimumlength))+minimumlength;
    filament(i,1,1)=xstart(i);
    filament(i,1,2)=ystart(i);
    theta1=rand*2*pi;
    filament(i,2,1)=xstart(i)+speed(i).*cos(theta1);
    filament(i,2,2)=ystart(i)+speed(i).*sin(theta1);
    
    for j=3:xylength(i)./speed(i)+1
        theta2=atan2((filament(i,j-1,2)-filament(i,j-2,2)),(filament(i,j-1,1)-filament(i,j-2,1)));
        newtheta=theta2+(0.5-rand)*bend;
        filament(i,j,1)=filament(i,j-1,1)+speed(i).*cos(newtheta);
        filament(i,j,2)=filament(i,j-1,2)+speed(i).*sin(newtheta);
    end
end
filament(:,:,1);
act_maximumlength = max(xylength);

% filament(1,end,1)-filament(1,1,1)
% filament(2,2,1)
% filament(2,1,1)
% filament(2,3,1)-filament(2,1,1)
%%
ismoving = rand(filamentcount,1) + fmot - 0.5;

figure
M.colormap = [];
Mwrite.colormap = []
for t=1:150
    M.cdata = zeros(height+outside, width+outside, 3);
    max(max(M.cdata))
    if rem(t, fmot_change_rate) == 0
        ismoving = (ismoving + fmot_change +  rand(filamentcount,1))/2;
    end
    for i=1:filamentcount
        if i <= filamentcount
        ender=length(filament(i,:,1))-find(abs(filament(i,end:-1:1,1))>0,1)+1;
        ender_list(i) = ender;
        theta2=atan2((filament(i,ender,2)-filament(i,ender-1,2)),(filament(i,ender,1)-filament(i,ender-1,1)));
        newtheta=theta2+(0.5-rand)*bend;

        if round(ismoving(i,1))
            filament(i,1:ender-1,1)=filament(i,2:ender,1);
            filament(i,1:ender-1,2)=filament(i,2:ender,2);
            filament(i,ender,1)=filament(i,ender,1)+speed(i).*cos(newtheta);
            filament(i,ender,2)=filament(i,ender,2)+speed(i).*sin(newtheta);
            if filament(i,ender,1) >= width+outside
                filament(i,ender,1) = 0;
                filament(i,ender,2) = 0;
                ender_list(i) = ender - 1;
            end
            if filament(i,ender,2) >= height+outside
                filament(i,ender,2) = 0;
                filament(i,ender,1) = 0;
                ender_list(i) = ender - 1;
            end
            
        else
            filament(i,1:xylength(i),1) = filament(i,1:xylength(i),1) + (rand()-0.5) * random;
            filament(i,1:xylength(i),2) = filament(i,1:xylength(i),2) + (rand()-0.5) * random;
        end
        
        if ender_list(i) == 1
            filament(i,:,:) = [];
            speed(i) = [];
            xylength(i) = [];
            ismoving(i) = [];
            filamentcount = filamentcount - 1;
            ender_list(i) = [];
        end
        end
    end
    
    set(gca,'NextPlot','replacechildren');

    for i=1:filamentcount
        plot(filament(i,1:ender_list(i),1),filament(i,1:ender_list(i),2),'w','LineWidth',1.8)
        hold on
    end
    hold off
    axis([0 width+outside 0 height+outside]);
    set(gca,'Units','Pixels');
    set(gcf,'Units','Pixels');
%     A=get(gca,'Position');
%     B=get(gcf,'Position');
%     set(gca,'Position',[70 70 639 479]);
%     set(gcf,'Position',[30 30 800 700]);
    set(gca,'Color',[0 0 0]);
    M=getframe;
    temp=uint8(double(M.cdata(:,:,1)));
    temp = imresize(temp,[480 720]);
    temp = 255-imnoise(255-temp, 'speckle',0.1);

    Mwrite.cdata(:,:,1)=temp;
    Mwrite.cdata(:,:,2)=temp;
    Mwrite.cdata(:,:,3)=temp;
    writeVideo(aviobj,Mwrite);
end
close
close(aviobj);

fprintf('DONE\n')
toc

