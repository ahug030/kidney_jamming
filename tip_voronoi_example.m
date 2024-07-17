cd('PATH TO WORKING DIRECTORY')
clear all;
close all;

%AreaLim is max allowable area of voronoi cells (for trimming edges)

cd('../2021-11-03 E16-18 kidney imaging_E17')
x=csvread('PATH TO CSV FILE CONTAINING TIP POSITIONS WITH X,Y COORDINATES AS COLUMNS');
tilestring=sprintf('FILENAME.png');
tilestring2=sprintf('FILENAME_shapeindex.png');
tilestring1a=sprintf('FILENAME_shapeindex.pdf');
tilestring3=sprintf('FILENAME_shapeindex.csv');
tilestring4=sprintf('FILENAME_edgelengths.csv');
AreaLim=0.07*10^5;

%flip y data to match fiji convention
x(:,2)=1331.2-x(:,2);
%x(:,2)=800.80-x(:,2);
%x(:,2)=1017.9-x(:,2);

%find voronoi neighbors
[ n, dim ] = size ( x );
%  V contains the Voronoi vertices,
%  C contains indices of Voronoi vertices that form the (finite sides of the)
%  Voronoi cells.
%
  [ V, C ] = voronoin ( x );
%
%  Two nodes are neighbors if they share an edge, that is, two Voronoi
%  vertices.
%
 % vn = sparse ( n, n );
 vn = zeros(n,n);
  for i = 1 : n
    for j = i + 1 : n
      s = size ( intersect ( C{i}, C{j} ) );
      if ( 1 < s(2) )
vn(i,j) = 1;
        vn(j,i) = 1;
      end
    end
  end
%find non-zero indices -- I,J are point numbers that interact
[I,J]=find(vn);
%for non-zeros, determine distances 

    bordercell_indices=[];
    bordercell_indices2=[];
    for i=1:size(x,1),
        
        %determine if cell is on the outside   
        pairs=[I(find(I==i),:),J(find(I==i),:)];
        linesvertices=[];
        for j=1:size(pairs,1),
            linesvertices=[linesvertices;[x(pairs(j,1),:),x(pairs(j,2),:)]];
        end;
        %for every pair of neighboring lines, get angles
        %determine quadrant that lines are in first and then use atan
        lines_angles=[];
        for y=1:(size(linesvertices,1)),
            if (linesvertices(y,4)-linesvertices(y,2))>0 && (linesvertices(y,3)-linesvertices(y,1))>0,
               lines_angles=[lines_angles;(atan(abs(linesvertices(y,4)-linesvertices(y,2))./abs(linesvertices(y,3)-linesvertices(y,1)))*(180/pi()))];
            elseif (linesvertices(y,4)-linesvertices(y,2))>0 && (linesvertices(y,3)-linesvertices(y,1))<0,
                lines_angles=[lines_angles;180-(atan(abs(linesvertices(y,4)-linesvertices(y,2))./abs(linesvertices(y,3)-linesvertices(y,1)))*(180/pi()))];
            elseif (linesvertices(y,4)-linesvertices(y,2))<0 && (linesvertices(y,3)-linesvertices(y,1))<0,
                lines_angles=[lines_angles;180+(atan(abs(linesvertices(y,4)-linesvertices(y,2))./abs(linesvertices(y,3)-linesvertices(y,1)))*(180/pi()))];
            else
                lines_angles=[lines_angles;360-(atan(abs(linesvertices(y,4)-linesvertices(y,2))./abs(linesvertices(y,3)-linesvertices(y,1)))*(180/pi()))];
            end;
        end;
        %sort angles with index
        b=[[1:size(linesvertices,1)]',lines_angles];
        c=sortrows(b,2);
        lines_angles_delta=[];
        %check difference in line angles for neighboring lines
        for w=1:(size(linesvertices,1)-1),
                lines_angles_delta = [lines_angles_delta;lines_angles(c(w+1,1),1)-lines_angles(c(w,1),1)];             
        end;
        %need to also wrap around the corner to get the angle between the
        %last line and first line
                lines_angles_delta=[lines_angles_delta;360-lines_angles(c(size(linesvertices,1),1),1)+lines_angles(c(1,1),1)];
        %check if angles are all below threshold
        if all(lines_angles_delta<135)==0,
            bordercell_indices=[bordercell_indices;i];
        %cell is a border cell
        else,
            %it's not
        end;         
   end;
   
outercells=x([bordercell_indices],:);
innercells=setdiff(x,outercells,'rows');

figure
plot(innercells(:,1),innercells(:,2),'o','MarkerSize',20,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
hold on
plot(outercells(:,1),outercells(:,2),'o','MarkerSize',20,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0 0 0]);
%axis([-100 100 -200 50])
axis square

% %calculate cell areas
% figure
% hold on;
% voronoi(x(:,1),x(:,2))
% A = zeros(length(C),1) ;
% for i = 1:length(C)
%     v1 = V(C{i},1) ; 
%     v2 = V(C{i},2) ;
%     patch(v1,v2,rand(1,3))
%     A(i) = polyarea(v1,v2) ;
% end

%plot convex hull of innercells
K=convhull(innercells);
plot(innercells(K,1),innercells(K,2))

%calculate cell areas
%select cells that have an non-nan area (are closed) and area below threshold
k1=figure;
hold on;
voronoi(x(:,1),x(:,2))
A = zeros(length(C),1) ;
A_true=[];
C_true=[];
v1_true=[];
v2_true=[];
original_points_true=[];
for i = 1:length(C)
    v1 = V(C{i},1) ; 
    v2 = V(C{i},2) ;
    A(i) = polyarea(v1,v2) ;
    if A(i)>0 && A(i)<AreaLim,
    %A_true=[A_true;A(i)];
    C_true=[C_true;C(i)];
    original_points_true=[original_points_true;x(i,:)];
    %v1_true=[v1_true;v1(i)];
    %v2_true=[v2_true;v2(i)];
    else,
    end;
    %patch(v1_true,v2_true,rand(1,3))
end

%create delauney data
dt = delaunayTriangulation(x(:,1),x(:,2));

coordnum=[];
shape_index=[];
%patchalpha=0.3;
patchalpha=1;
A_true = zeros(length(C_true),1);
edgel_list=[];
for i = 1:length(C_true)
    v1 = V(C_true{i},1) ; 
    v2 = V(C_true{i},2) ;
    P=polyshape(v1,v2);
    coordnum=[coordnum;numsides(P)];
    shape_index=[shape_index;perimeter(P)/sqrt(polyarea(v1,v2))];
    for k=1:(size(P.Vertices,1)-1),
    edgel=pdist([P.Vertices(k,:);P.Vertices(k+1,:)]);
    edgel_list=[edgel_list;edgel];
    end;
    %catch the last edge between last vertex and first
    edgel=pdist([P.Vertices(size(P.Vertices,1),:);P.Vertices(1,:)]);
    edgel_list=[edgel_list;edgel];
    
    A_true(i) = polyarea(v1,v2);
    
    %color patches with coord number=4 as orange
    if size(v1,1)==4,
    patch(v1,v2,[255/255 110/255 0/255],'FaceAlpha',patchalpha);
    else,
    end;
    %color patches with coord number=5 as magenta
    if size(v1,1)==5,
    patch(v1,v2,[255/255 44/255 149/255],'FaceAlpha',patchalpha);
    else,
    end;
    %color patches with coord number=7 as green
    if size(v1,1)==7,
    patch(v1,v2,[12/255 210/255 30/255],'FaceAlpha',patchalpha);
    else,
    end;
    %color patches with coord number=8 as purple
    if size(v1,1)==8,
    patch(v1,v2,[177/255 0/255 236/255],'FaceAlpha',patchalpha);
    else,
    end;
    %color all other patches gray
    if size(v1,1)~=4 && size(v1,1)~=5 && size(v1,1)~=7 && size(v1,1)~=8,
    patch(v1,v2,[0.3 0.3 0.3],'FaceAlpha',patchalpha);
    else,
    end;
end
axis equal
axis ([0 1331 0 1331]);
%axis ([0 933.40 0 800.8]);
%axis ([0 824.2 0 1017.9]);
%model cases...
% axis ([200 1700 200 1700]);
%exportgraphics(gcf,'barchart.png','Resolution',300)
print(gcf, '-dpng', tilestring);
exportgraphics(gcf,tilestring1a,'ContentType','vector')

figure
hold on;

coordnum=[];
xccyccholder=[];
xcycholder=[];
A_true = zeros(length(C_true),1);
for i = 1:length(C_true)
    v1 = V(C_true{i},1) ; 
    v2 = V(C_true{i},2) ;
    P=polyshape(v1,v2);
    coordnum=[coordnum;numsides(P)];
    A_true(i) = polyarea(v1,v2) ;
    [xc,yc]=centroid(P);
    xcyc=[xc,yc];
    xcycholder=[xcycholder;xcyc];
    %finding closest node in delaunay
    kp=dsearchn(dt.Points,[xc,yc]);
    xccycc=[dt.Points(kp,1),dt.Points(kp,2)];
    xccyccholder=[xccyccholder;xccycc];
 
    
    %color delaunay node with coord number=4 as orange
    if size(v1,1)==4,
    %patch(v1,v2,[1 1 0]);
    plot(xccycc(1),xccycc(2),'o','MarkerFaceColor',[255/255 110/255 0/255],'MarkerEdgeColor',[255/255 110/255 0/255],'MarkerSize',10);
    else,
    end;
    %color delaunay node with coord number=5 as magenta
    if size(v1,1)==5,
    %patch(v1,v2,[1 0 0]);
    plot(xccycc(1),xccycc(2),'o','MarkerFaceColor',[255/255 44/255 149/255],'MarkerEdgeColor',[255/255 44/255 149/255],'MarkerSize',10);
    else,
    end;
    %color delaunay node with coord number=7 as green
    if size(v1,1)==7,
    %patch(v1,v2,[0 1 0]);
    plot(xccycc(1),xccycc(2),'o','MarkerFaceColor',[12/255 210/255 30/255],'MarkerEdgeColor',[12/255 210/255 30/255],'MarkerSize',10);
    else,
    end;
    %color delaunay node with coord number=8 as purple
    if size(v1,1)==8,
    %patch(v1,v2,[0 1 1]);
    plot(xccycc(1),xccycc(2),'o','MarkerFaceColor',[177/255 0/255 236/255],'MarkerEdgeColor',[177/255 0/255 236/255],'MarkerSize',10);
    else,
    end;
    %color all other delaunay nodes gray
    if size(v1,1)~=4 && size(v1,1)~=5 && size(v1,1)~=7 && size(v1,1)~=8,
    %patch(v1,v2,[0.5 0.5 0.5]);
    plot(xccycc(1),xccycc(2),'o','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3],'MarkerSize',10);
    else,
    end;
end

%add delauney on top
H1=triplot(dt);
set(H1,'Color',[0.5 0.5 0.5]);




%interpolated heat map of shape index, n.b. xccycc represents at tip
%position, xcyc represents at centroid of voronoi cells

%select if plotting shape index map under voronoi cells, plotting under
%delauney otherwise
%figure(k1)

%locally averaging over window size of around 3x3 tips... or not


[xq,yq] = meshgrid(0:30:1331,0:30:1331);
%model cases...
% [xq,yq] = meshgrid(0:30:1800,0:30:1800);

%vq = griddata(xccyccholder(:,1),xccyccholder(:,2),shape_index,xq,yq,'cubic');
vq = griddata(xcycholder(:,1),xcycholder(:,2),shape_index,xq,yq,'v4');
%vq_smooth = imfilter(vq,fspecial('average',[7 7]),3.9);
vq_smooth = vq;
vqq=imresize(vq_smooth,1331/45);
%model cases...
% vqq=imresize(vq_smooth,1800/61);

%cutting colormap at shape factor of 3.8, 4
imagesc(vqq,[3.8 4]);
%imagesc(flipud(vqq),[3.81 4]);
colormap('gray');

chi=get(gca, 'Children');
%Reverse the stacking order so that the patch overlays the line
set(gca, 'Children',flipud(chi));

axis equal
axis ([0 1331 0 1331]);
%model cases...
% axis ([200 1700 200 1700]);

print(gcf, '-dpng', tilestring2);



figure
histogram(coordnum);

coordnum_mean=mean(coordnum)

%run for many samples & produce plots of metrics... need to interpret from
%Brojan, Irvine, Jimenez.

%total defect charge Q=sum(s) where s=6-Z where Z is coordnum

s=6-coordnum;
Q=sum(s)

figure;
histogram(shape_index);
shape_index_median=median(shape_index)
shape_index_mean=mean(shape_index)
csvwrite(tilestring3,shape_index);
csvwrite(tilestring4,edgel_list);

voronoi_data_out=[original_points_true coordnum shape_index A_true];
csvwrite('voronoi_data_out.csv',voronoi_data_out);

%Jimenez and Irvine plots of Q against curvature are a little different
   %Jimenez uses integrated curvature units, which for a sphere is 1/r2
   %times 4*pi*r2 = 4*pi units of curvature
   %Irvine normalizes by 3/pi = 12 units of disclination