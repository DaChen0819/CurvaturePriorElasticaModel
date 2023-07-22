clc;close all;clear variables;
numOriens=72;

load('data/rawImg.mat');
load('data/B.mat');load('data/X.mat');load('data/Y.mat');load('data/Z.mat');
normalizedImg=im2double(rawImg(:,:,1));

% for better visualization, we only use 72 orientations.
normalizedOS=normalizedOS(:,:,1:2:end);

mSampledPt=[560;469];
mSampledPt2=[385;670];
figure(1);imagesc(normalizedImg,[0 1.05]);axis image;axis off;colormap(gray);hold on;
h1=plot(mSampledPt(1),mSampledPt(2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',8);
h2=plot(mSampledPt2(1),mSampledPt2(2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','Markersize',8);


[numRows,numCols,~]=size(normalizedOS);
cost=exp(-3.*normalizedOS);
figure(2);list=[6,18,36];list2=[0.4,0.4,0.5];

for i=1:3
    s=surf(X(:,:,i),Y(:,:,i),Z(:,:,i),B(:,:,i),'EdgeColor','None','HandleVisibility','on','FaceAlpha',list2(i));
    hold on;
end
view([-33 26]);shading interp;colormap(othercolor('Mlakecolors'));%Mlakecolors %Msunsetcolors %Msouthwestcolors
h=colorbar('AxisLocation','in','Ticks',[0.1  0.95],'TickLabels',{'Low','High'});
h.TickDirection='in';
zlim([6 36]);zticks([6 18 36]);zticklabels({'\pi/6','\pi/2','\pi'});
xlim tight;xticks([]);xticklabels('');
ylim tight;yticks('');yticklabels('');
xlabel('x','FontSize',18);ylab=ylabel('y','FontSize',18);z=zlabel('$\theta$','FontSize',18);
z.Rotation=0;
z.Interpreter='latex';
axis square;grid on;
ax = gca;ax.FontSize=18;
box on;
set(gca,'Linewidth',1)



