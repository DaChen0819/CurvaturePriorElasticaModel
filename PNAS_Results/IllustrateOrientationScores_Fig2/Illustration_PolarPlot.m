clc;close all;clear variables;

numOriens=144;
mSampledPt=[560;469];
mSampledPt2=[385;670];

load('data/listOS1.mat');
load('data/listOS2.mat');
theta=linspace(0,2*pi,numOriens+1);
h=polarplot(theta(1:numOriens+1),listOS1(1:numOriens+1));
h.Color='red';
h.LineWidth=2;
pax=gca;
thetaticks([0 90 180 270])
thetaticklabels({'0','\pi/2','\pi','3\pi/2'})
pax.FontSize=15;
pax.RTick=[0,0.5,0.8,1.1];
pax.RTickLabel={'0';'0.5';'0.8';'1.1'};
pax.RColor=addcolor(258);
pax.GridAlpha=0.5;
hold on;
h2=polarplot(theta(1:numOriens+1),listOS2(1:numOriens+1));
h2.Color='blue';
h2.LineWidth=2;

