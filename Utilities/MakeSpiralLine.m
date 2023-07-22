function Spiral=MakeSpiralLine(a,b,maxSize)

% maxSize=4*pi+0.5*pi;
% a=3;
% b=6;
theta=0.5*pi:0.001*pi:maxSize;
coe=a+b.*theta;
x=coe.*cos(-theta)+95;
y=coe.*sin(-theta)+105;

% close all;
% figure(1);imagesc(ones(191,191),[0,1]);colormap gray;axis image;hold on;
% plot(x,y,'r');

Spiral=[x;y;-theta+pi/2];

end