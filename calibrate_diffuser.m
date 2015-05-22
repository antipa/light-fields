%% Display pattern and trigger camera

im = zeros(501,101);        %Make empty image
h1 = figure(1)  %Create figure window with handle
set(h1,'MenuBar','none')    %Clear menu bar from window to reduce stray light
c = 50; %Column about which the pattern will be centered
r = 450   %Row center
M = 50; %Pattern size in row direction (first dimension)
N = 50; %Size in column direction
dwell = 2;  %Wait time between images (in seconds)
for m = 1:M     
    for n = 1:N
        set(0,'CurrentFigure',h1)    %Set figure 1 to active without bringing into focus
        im(round(r+M/2-m),round(c+N/2-n)) = 1;  %Turn on pixel we are interested in
        imshow(im);%    %Display
        set(h1,'Color',[0 0 0]) %Set background to black
        axis image %Image axis makes it display pixel-for-pixel
        pause(dwell)    %Wait
        
        %---------Trigger  camera here----------
        
        im(round(r+M/2-m),round(c+N/2-n)) = 0;  %Turn off pixel
    end
end