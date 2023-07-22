clear variables; clc; close all;

numGrids=151;
numOriens=72;
activeRegion=ones(numGrids,numGrids,numOriens);
lambda=2*pi;

% Defining the domain
input.order=1.0;
input.seedRadius=0;
input.eps=0.1;
input.dims = [numGrids;numGrids;numOriens];
input.arrayOrdering='YXZ_ColumnMajor';
input.verbosity=0;
input.origin=[0;0]; 
input.gridScale=2*pi/numOriens; 
input.cost = 1; 
input.xi = 4; 

seedPos=round([0.9*numGrids;0.1*numGrids]);
seedOrien=2*pi/5;
rescaledSeeds=RescaledCoords(seedPos,input.origin,[input.gridScale;input.gridScale],false);
input.seeds=[rescaledSeeds;seedOrien];

tipPos=round([0.4*numGrids;0.6*numGrids]);
tipOrien=pi/4;
rescaledTipPos=RescaledCoords(tipPos,input.origin,[input.gridScale;input.gridScale],false);
input.tips = [rescaledTipPos',tipOrien.*ones(size(rescaledTipPos,2),1)]'; % where the geodesics end. [x;y;theta]

input.exportValues=1; % distance table, of size [n,n,numberOfDirections]
input.model='ElasticaExt2'; %Alternatively : 'ReedsSheppExt2','ReedsSheppForwardExt2', 'ElasticaExt2', 'DubinsExt2'
input.geodesicSolver='ODE';
input.geodesicStep=0.5;
input.geodesicCausalityTolerance=4;
input.geodesicTargetTolerance=6;

try
input.kappa=0.*activeRegion; 
output=HamiltonCurvaturePriorMinimalPaths(input);
cGeodesics = mat2cell(output.geodesicPoints,3,output.geodesicLengths); 

input.kappa=(1/6).*activeRegion; 
output1=HamiltonCurvaturePriorMinimalPaths(input);
cGeodesics1 = mat2cell(output1.geodesicPoints,3,output1.geodesicLengths); 

input.kappa=(1/3)*activeRegion; 
output2=HamiltonCurvaturePriorMinimalPaths(input);
cGeodesics2 = mat2cell(output2.geodesicPoints,3,output2.geodesicLengths);

input.kappa=(1/2)*activeRegion; 
output3=HamiltonCurvaturePriorMinimalPaths(input);
cGeodesics3 = mat2cell(output3.geodesicPoints,3,output3.geodesicLengths);

catch
load('data/cGeodesics.mat');
load('data/cGeodesics1.mat');
load('data/cGeodesics2.mat');
load('data/cGeodesics3.mat');
end

green1=addcolor(198);
yellow1=addcolor(158);

l=18;
figure(1);imagesc(activeRegion(:,:,1),[0 1]); axis image;hold on;colormap gray;grid off;
rescaledGeodesic=RescaledCoords(cGeodesics{1}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
rescaledGeodesic1=RescaledCoords(cGeodesics1{1}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
rescaledGeodesic2=RescaledCoords(cGeodesics2{1}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
rescaledGeodesic3=RescaledCoords(cGeodesics3{1}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
h1=plot(rescaledGeodesic(1,:),rescaledGeodesic(2,:),  'Color','b','LineWidth',2);
h2=plot(rescaledGeodesic1(1,:),rescaledGeodesic1(2,:),'r','LineWidth',2);
h3=plot(rescaledGeodesic2(1,:),rescaledGeodesic2(2,:),'Color',yellow1,'LineWidth',2);
h4=plot(rescaledGeodesic3(1,:),rescaledGeodesic3(2,:),'Color',green1,'LineWidth',2);
h5=func_DrawArrow(seedPos,seedPos+l*[cos(seedOrien);sin(seedOrien)],25,'BaseAngle',30,'Width',3); 
set(h5,'Facecolor','r','EdgeColor','r');
h6=func_DrawArrow(tipPos,tipPos+l*[cos(tipOrien);sin(tipOrien)],25,'BaseAngle',30,'Width',3);
set(h6,'Facecolor','b','EdgeColor','b');
plot(seedPos(1),seedPos(2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',10);
plot(tipPos(1),tipPos(2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','Markersize',10);
h=legend([h1 h2 h3 h4],'$\omega\equiv0$','$\omega\equiv 1/6$','$\omega\equiv 1/3$','$\omega\equiv1/2$');
set(h,'Interpreter','latex','fontsize',22,'Location','northwest');
set(gca,'XTickLabel','');set(gca,'YTickLabel','');
seedPosShifted=seedPos+[-10;0];
tipPosShifted =tipPos +[3;-10];
text(seedPosShifted(1),seedPosShifted(2),'$s$','FontSize',22,'Interpreter','latex');
text(tipPosShifted(1),tipPosShifted(2),'$y$','FontSize',22,'Interpreter','latex');
alpha(h5,0.4)



smoothWidth=1;
figure(2);imagesc(repmat(activeRegion(:,:,1),[1 1 3])); axis image;hold on;grid off;
c=CurvatureEstimationFromPath(cGeodesics{1},smoothWidth);
c1=CurvatureEstimationFromPath(cGeodesics1{1},smoothWidth);
c2=CurvatureEstimationFromPath(cGeodesics2{1},smoothWidth);
c3=CurvatureEstimationFromPath(cGeodesics3{1},smoothWidth);
patch(cat(2,rescaledGeodesic(1,:),nan),cat(2,rescaledGeodesic(2,:),nan),[c 0]',...
        'edgecolor','interp','MarkerFaceColor','flat','LineWidth',2);
patch(cat(2,rescaledGeodesic1(1,:),nan),cat(2,rescaledGeodesic1(2,:),nan),[c1 0]',...
        'edgecolor','interp','MarkerFaceColor','flat','LineWidth',2);
patch(cat(2,rescaledGeodesic2(1,:),nan),cat(2,rescaledGeodesic2(2,:),nan),[c2 0]',...
        'edgecolor','interp','MarkerFaceColor','flat','LineWidth',2);
patch(cat(2,rescaledGeodesic3(1,:),nan),cat(2,rescaledGeodesic3(2,:),nan),[c3 0]',...
        'edgecolor','interp','MarkerFaceColor','flat','LineWidth',2);
plot(seedPos(1),seedPos(2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',7);
plot(tipPos(1),tipPos(2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','Markersize',7);
colorbar('south');
colormap jet;
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');





