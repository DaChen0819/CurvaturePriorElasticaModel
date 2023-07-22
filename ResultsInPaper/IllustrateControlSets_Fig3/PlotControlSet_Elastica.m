close all;clear variables;clc;
beta=4;
s=0:0.002:1;
t=1:-0.002:0;

figure(1);hold on;
a1=sqrt(1/4-(s-0.5).^2)/beta;
b1=-sqrt(1/4-(t-0.5).^2)/beta;
h1=plot(cat(2,s,t,0),cat(2,a1,b1,0));
h1.LineWidth=3;
h1.Color=addcolor(195);

beta=2;
a2=sqrt(1/4-(s-0.5).^2)/beta;
b2=-sqrt(1/4-(t-0.5).^2)/beta;
h2=plot(cat(2,s,t,0),cat(2,a2,b2,0));
h2.LineWidth=3;
h2.Color=addcolor(165);

beta=1;
a3=sqrt(1/4-(s-0.5).^2)/beta;
b3=-sqrt(1/4-(t-0.5).^2)/beta;
h3=plot(cat(2,s,t,0),cat(2,a3,b3,0));
h3.LineWidth=3;
h3.Color=addcolor(260);

ax=gca;
ax.XAxisLocation='origin';
axis image;
ax.YLim=[-0.55 0.55];
ax.XLim=[0 1.2];
ax.LineWidth=1;
ax.Layer='top';
ax.FontSize=12;

yticks([-0.5,-0.2,0,0.2,0.5]);
yticklabels({'-0.5','-0.2','0','0.2','0.5'});
xticks(0.2:0.2:1);
xticklabels({'0.2','0.4','0.6','0.8','1.0'});
x=xlabel('$\dot\nu$');x.Interpreter='latex';x.FontSize=25;x.FontWeight='bold';
x.Units='points';
y=ylabel('$\dot\theta$');y.Interpreter='latex';y.FontSize=25;y.FontWeight='bold';
x.Units='points';
y.Rotation=0;
box off;
grid off;
h4=legend([h1 h2 h3],{'\beta = 4','\beta = 2','\beta = 1'});
h4.Box='off';
h4.FontSize=15;







