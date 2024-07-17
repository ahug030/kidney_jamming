cd('PATH TO WORKING DIRECTORY')
clear all;
close all;

pix2micron=1.3;
zextent=0.18;

stlstring=sprintf('PATH TO STL FILE DESCRIBING KIDNEY SURFACE');
tilestring0=sprintf('PATH TO CSV FILE CONTAINING TIP POSITIONS WITH X,Y COORDINATES AS COLUMNS');
%outputs image of tip positions
tilestring1=sprintf('FILENAME_tips.png');
tilestring1a=sprintf('FILENAME_tips.pdf');
%outputs image of height map
tilestring2=sprintf('FILENAME_height.png');

x=csvread(tilestring0);
x(:,2)=(1331-x(:,2));

figure;
hold on;
plot(x(:,1),x(:,2),'o','MarkerSize',12,'MarkerFaceColor',[12/255 210/255 30/255],'MarkerEdgeColor',[12/255 210/255 30/255]);
axis square
axis([0 1331 0 1331]);
set(gca,'TickDir','out');

F1 = getframe(gca);
imwrite(F1.cdata, tilestring1);
exportgraphics(gcf,tilestring1a,'ContentType','vector')

stl_tr=stlread(stlstring);
figure;
trisurf(stl_tr,'EdgeColor','none');

colormap('gray');
axis equal
%axes are in mm
%caxis([-0.2 0]);
stlmidpoint=mean([max(stl_tr.Points);min(stl_tr.Points)]);
axis([stlmidpoint(1)-1.33/2 stlmidpoint(1)+1.33/2 stlmidpoint(2)-1.33/2 stlmidpoint(2)+1.33/2 min(stl_tr.Points(:,3)) min(stl_tr.Points(:,3))+zextent]);
caxis([min(stl_tr.Points(:,3)) min(stl_tr.Points(:,3))+zextent]);
view([0 90]);
grid off
set(gca,'TickDir','out');

F2 = getframe(gca);
imwrite(F2.cdata, tilestring2);
