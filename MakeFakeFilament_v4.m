clear all;
close all;

% --- This part creates traces with bent random bending, then fits them to 
% splines, which it then stores.

filamentcount=100;
trace_speed = 3
thickness=1;

height = 480;
width = 720;

tracelength = 300
time_trace = 1:tracelength+2
bend = 1;

h1 = figure
h2 = figure

for i=1:filamentcount
    xstart(i)=rand*width;
    ystart(i)=rand*height;
    traces(i,1,1)=xstart(i);
    traces(i,1,2)=ystart(i);
    theta1=rand*2*pi;
    traces(i,2,1)=xstart(i)+trace_speed.*cos(theta1);
    traces(i,2,2)=ystart(i)+trace_speed.*sin(theta1);

for t = 1:tracelength
    
    ender=2;
    theta2=atan2((traces(i,t+1,2)-traces(i,t+0,2)),(traces(i,t+1,1)-traces(i,t+0,1)));
    newtheta=theta2+(0.5-rand)*bend;
    
    traces(i,t+2,1)=traces(i,t+1,1)+trace_speed.*cos(newtheta);
    traces(i,t+2,2)=traces(i,t+1,2)+trace_speed.*sin(newtheta);
    
    axis([0 width 0 height]);
        
end

figure(h2)
subplot(2,1,1)
plot(time_trace,traces(i,:,1))
hold on
s1 = spline(time_trace,traces(i,:,1));
plot(time_trace,ppval(s1, time_trace),'r+')

subplot(2,1,2)
plot(time_trace,traces(i,:,2))
hold on
s2 = spline(time_trace,traces(i,:,2));
plot(time_trace,ppval(s2,time_trace),'r+')

figure(h1)
plot(ppval(s1,time_trace),ppval(s2,time_trace),'b','LineWidth',1.8)
hold on
plot(traces(i,:,1),traces(i,:,2),'b+')
splines{1,i} = s1;
splines{2,i} = s2;
[num2str(i) ' out of ' num2str(filamentcount) ' is done']
end

%%
% --- Here we actually let filaments walk on these traces. Here we vary the
% size and the speed of the filaments
tic


filamentcount= 100
max_frame = 50
min_fil_length = 2
max_fil_length = 20
filament_length = zeros(filamentcount,1) + 10
filament_length = rand(filamentcount,1)*(max_fil_length - min_fil_length) + min_fil_length

speed = zeros(filamentcount, max_frame);

% define the change in speed in terms of hill curves
hill = @(is_inc, vmax, halftime, n, time)...
    vmax*((is_inc==0) + ( (is_inc==1) - 0.5 ) * 2 .* time.^n ./( halftime.^n + time.^n));

vmax = 2;
n = 2;

for i = 1:filamentcount
    time{i} = 1:filament_length(i,1);
    length(time{i});
    halftime_list(i) = max_frame./2*(rand*0.5+0.75)+10
    speed(i,:) = hill(0, vmax/5, halftime_list(i), n, 0:1:max_frame-1);
%     speed(i,:) = 1;
    plot(speed(i,:)*5)
    hold on
end

mean_halftime = mean(halftime_list);


fmot = 1; % between [0, 1.0], at 0.5, 50% of filaments starts moving
fmot_change = 0; % between [-0.5, 0.5]. E.g. for 0.5 all filaments will move by the end
fmot_change_rate = 0; % between [1,100]. 1 changes fastest, 5 changes slower, and so on
speckle = 0.05; % speckle noise on the hole image
ismoving = rand(filamentcount,1) + fmot - 0.5;


output = ['vmax', num2str(vmax), '_vHalf',num2str(mean_halftime),...
    '_n', num2str(n), '_fmot',num2str(fmot),...
    '_fmotChange',num2str(fmot_change), '_fmotChangeRate',num2str(fmot_change_rate),...
     '_speckle',num2str(speckle)]

% output = ['file.avi']

aviobj = VideoWriter([output '.avi']);
open(aviobj)
M.colormap = [];
Mwrite.colormap = [];

figure


for tt = 1:max_frame
    
    if rem(tt, fmot_change_rate) == 0
        ismoving = (ismoving + fmot_change +  rand(filamentcount,1))/2;
    end
    set(gca,'NextPlot','replacechildren');
    
    for i = 1:filamentcount

        plot(ppval(splines{1,i},time{i}),ppval(splines{2,i},time{i}),'w','LineWidth',1.2)
        hold on
        speed(i,tt) = round(ismoving(i,1)) * speed(i,tt);
        % move
        time{i} = time{i} + speed(i,tt);
        
    end
    
    axis([0 width 0 height]);
    set(gca,'Units','Pixels');
    set(gcf,'Units','Pixels');
    set(gca,'Color',[0 0 0]);
    
%     pause(0.1)
    
    M=getframe;
    temp=uint8(double(M.cdata(:,:,1)));
    temp = imresize(temp,[480 720]);
    temp = 255-imnoise(255-temp, 'speckle',speckle);

    Mwrite.cdata(:,:,1)=temp;
    Mwrite.cdata(:,:,2)=temp;
    Mwrite.cdata(:,:,3)=temp;
    
    writeVideo(aviobj,Mwrite);
    
end


% close
close(aviobj);

fprintf('DONE\n')
toc

figure
for i = 1:filamentcount
    plot(speed(i,:))
    hold on
end

ylim([0 vmax/5*1.1])
print('-depsc',['vel_', output '.eps'])
savefig(['velFigs_', output '.fig'])



