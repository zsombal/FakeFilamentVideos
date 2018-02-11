clear all;
close all;
tic

speed=1; %inverted
fmot = 1; % between [0, 1.0], at 0.5, 50% of filaments starts moving
fmot_change = 0; % between [-0.5, 0.5]. E.g. for 0.5 all filaments will move by the end
fmot_change_rate = 0; % between [1,100]. 1 changes fastest, 5 changes slower, and so on
random = 0; % random motion of the non-moving filaments
speckle = 0.0; % speckle noise on the hole image

height = 480;
width = 720;
outside = 400;

bend=0;
filamentcount=50;
thickness=1;
minimumlength=20;
maximumlength=20;

output = ['speed', num2str(speed), '_fmot',num2str(fmot),...
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
    theta1=0%pi/4%rand*2*pi;
    filament(i,2,1)=xstart(i)+cos(theta1);
    filament(i,2,2)=ystart(i)+sin(theta1);
    
    for j=3:xylength(i)
        theta2=atan2((filament(i,j-1,2)-filament(i,j-2,2)),(filament(i,j-1,1)-filament(i,j-2,1)));
        newtheta=theta2+(0.5-rand)*bend;
        filament(i,j,1)=filament(i,j-1,1)+cos(newtheta);
        filament(i,j,2)=filament(i,j-1,2)+sin(newtheta);
    end
end

act_maximumlength = max(xylength);

%%
% ismoving = zeros(filamentcount,1);
ismoving = rand(filamentcount,1) + fmot - 0.5;

M.colormap = [];

for t=1:150
    M.cdata = zeros(height+outside, width+outside, 3);

    if rem(t, fmot_change_rate) == 0
        ismoving = (ismoving + fmot_change +  rand(filamentcount,1))/2;
    end
    for i=1:filamentcount
        ender=length(filament(i,:,1))-find(abs(filament(i,end:-1:1,1))>0,1)+1;
        theta2=atan2((filament(i,ender,2)-filament(i,ender-1,2)),(filament(i,ender,1)-filament(i,ender-1,1)));
        newtheta=theta2+(0.5-rand)*bend;
        if round(ismoving(i,1))
            filament(i,1:ender-1,1)=filament(i,2:ender,1);
            filament(i,1:ender-1,2)=filament(i,2:ender,2);
            filament(i,ender,1)=filament(i,ender,1)+cos(newtheta);
            filament(i,ender,2)=filament(i,ender,2)+sin(newtheta);
        else
            filament(i,1:xylength(i),1) = filament(i,1:xylength(i),1) + (rand()-0.5) * random;
            filament(i,1:xylength(i),2) = filament(i,1:xylength(i),2) + (rand()-0.5) * random;
        end

    end
    
%     set(gca,'NextPlot','replacechildren');
    M_lin = reshape(M.cdata(:,:,1), 1, (height+outside)*(width+outside));
    for i=1:filamentcount
%         xylength(i)
%         lengthfil= length(filament(i,:,1))-find(filament(i,end:-1:1,1)>0,1)+1
        M_lin(sub2ind(size(M.cdata(:,:,1)), round(filament(i,1:xylength(i)-1,2)),round(filament(i,1:xylength(i)-1,1)))) = 255;
%         
%         plot(filament(i,1:lengthfil,1),filament(i,1:lengthfil,2),'w','LineWidth',thickness)
%         hold on
    end
    
    out = imgaussfilt(reshape(M_lin, height+outside, width+outside), 1);
    M.cdata(:,:,1) = out;
%     M_speckle = imnoise(zeros(height+outside, width+outside, 1)+1/2, 'speckle', 1);
%     M.cdata(:,:,1) = 255-imnoise(255-uint8(double(out)), 'speckle',0.000005)
%     imshow(M_speckle)
%     imshow(M.cdata(:,:,1));
    %%
    if t>1
        Mlast=M;
    end
    if t>1
        for i=1:speed
            
%             M_speckle = imnoise(zeros(height, width, 1), 'speckle', 1);
%             imshow(double(M_speckle))
%             temp=uint8(((double(M.cdata(:,:,1))*i+double(Mlast.cdata(:,:,1))*(speed-i)))/speed);
            temp=uint8(((double(M.cdata(outside/2:outside/2+height-1,outside/2:outside/2+width-1,1))*i+double(Mlast.cdata(outside/2:outside/2+height-1,outside/2:outside/2+width-1,1))*(speed-i)))/speed);
            temp = 255-imnoise(255-temp, 'speckle',speckle);
            %             Mwrite.cdata=M;
%             temp = imadd(temp, uint8(double(M_speckle)));
%             figure
%             imshow(temp);
            
            Mwrite.cdata(:,:,1)=temp;
            Mwrite.cdata(:,:,2)=temp;
            Mwrite.cdata(:,:,3)=temp;
            Mwrite.colormap = [];
%             imshow(Mwrite.cdata);
  
            writeVideo(aviobj,Mwrite);
        end
    end
end
close
close(aviobj);


fprintf('DONE\n')
toc

