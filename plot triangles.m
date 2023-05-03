clear all
clc
% load('x_coor.mat')
% load('y_coor.mat')
% load('sigma_x.mat')
% load('sigma_y.mat')
% load('sigma_xy.mat')
% load('x_coor_sigmay.mat')
% load('y_coor_sigmay.mat')
% load('x_coor_sigmaxy.mat')
% load('y_coor_sigmaxy.mat')
% F_sigmax=scatteredInterpolant(x_coor,y_coor,sigma_x);
% F_sigmax.Method='natural';
% F_sigmay=scatteredInterpolant(x_coor_sigmay,y_coor_sigmay,sigma_y);
% F_sigmay.Method='natural';
% F_sigmaxy=scatteredInterpolant(x_coor_sigmaxy,y_coor_sigmaxy,sigma_xy);
% F_sigmaxy.Method='natural';
%--------------------------------------
%------- find integration points and corresponding weights
%--------------------------------------
format long
num_domain=8;						%the number of domains
Points=zeros(3,2,num_domain);
% Points(:,:,1)=[0 0.4; 0.64 0.8; 0 0.8];				%from area A1
% Points(:,:,2)=[0 0.4; 0.64 0.4; 0.64 0.8];			%from area A2
% Points(:,:,3)=[0 0; 0.64 0; 0 0.4];					%from area B1
% Points(:,:,4)=[0.64 0; 0.64 0.4; 0 0.4];			%from area B2
% Points(:,:,5)=[0.64 0; 0.8 0; 0.8 0.4];				%from area C1
% Points(:,:,6)=[0.64 0; 0.8 0.4; 0.64 0.4];			%from area C2
% Points(:,:,7)=[0.8 0.4; 0.8 0.8; 0.64 0.8];			%from area D1
% Points(:,:,8)=[0.64 0.4; 0.8 0.4; 0.64 0.8];		%from area D2
Points(:,:,1)=[0 0.4; 0.5 0.8; 0 0.8];				%from area A1
Points(:,:,2)=[0 0.4; 0.5 0.4; 0.5 0.8];			%from area A2
Points(:,:,3)=[0 0; 0.5 0; 0 0.4];					%from area B1
Points(:,:,4)=[0.5 0; 0.5 0.4; 0 0.4];			%from area B2
Points(:,:,5)=[0.5 0; 0.8 0; 0.8 0.4];				%from area C1
Points(:,:,6)=[0.5 0; 0.8 0.4; 0.5 0.4];			%from area C2
Points(:,:,7)=[0.8 0.4; 0.8 0.8; 0.5 0.8];			%from area D1
Points(:,:,8)=[0.5 0.4; 0.8 0.4; 0.5 0.8];		%from area D2
ConnectivityList=[1 2 3];
N=2;						%the number of refinements    
% num_triangles=4^N;			%the number of triangular elements
% x1=zeros(num_domain,num_triangles);			%x coordinate of the nodes in triangular elements
% y1=zeros(num_domain,num_triangles);			%y coordinate of the nodes in triangular elements
% x2=zeros(num_domain,num_triangles);
% y2=zeros(num_domain,num_triangles);
% x3=zeros(num_domain,num_triangles);
% y3=zeros(num_domain,num_triangles);
% I=0;								%integration result
d=19;								%degree
% num_intpoints=73;					%the number of integartion points in each element
% intpoints=cell(num_domain,num_triangles,2);	
% for j=1:num_domain
	% for i=1:num_triangles
		% intpoints{j,i,1}=zeros(num_intpoints,2);		%storing the coordinates of integartion points in each element
		% intpoints{j,i,2}=zeros(num_intpoints,1);		%storing the weights of the above points
	% end	
% end
% intpoints_sigma=zeros(num_domain,num_triangles,num_intpoints,3);		%storing interpolated sigma
j=1;
while j<num_domain+1
	TR = triangulation(ConnectivityList,Points(:,:,j));
	figure()
	hold on; 
	box on;
	title(['Original triangulation TR',num2str(j)]);
	triplot(TR,'k','LineWidth',1)  
	[rTR{1:N,j}] = trirefine(TR,'NumberOfRefinements',N);
	figure()
	hold on; 
	box on;
	title(['refine TR',num2str(j)]);
	triplot(rTR{N,j},'k','LineWidth',1)      
	% % for i=1:num_triangles
		% % x1(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,1),1);
		% % y1(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,1),2);
		% % x2(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,2),1);
		% % y2(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,2),2);
		% % x3(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,3),1);
		% % y3(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,3),2);
	% % end
	% x1(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,1),1);
	% y1(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,1),2);
	% x2(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,2),1);
	% y2(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,2),2);
	% x3(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,3),1);
	% y3(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,3),2);
	x1(j,:)=rTR{N,j}.Points(rTR{N,j}.ConnectivityList(:,1),1);
	y1(j,:)=rTR{N,j}.Points(rTR{N,j}.ConnectivityList(:,1),2);
	x2(j,:)=rTR{N,j}.Points(rTR{N,j}.ConnectivityList(:,2),1);
	y2(j,:)=rTR{N,j}.Points(rTR{N,j}.ConnectivityList(:,2),2);
	x3(j,:)=rTR{N,j}.Points(rTR{N,j}.ConnectivityList(:,3),1);
	y3(j,:)=rTR{N,j}.Points(rTR{N,j}.ConnectivityList(:,3),2);
	% % f=@(x,y) (2.*x+y).^3;
	% tic
	% for i=1:num_triangles
		% Q = quadtriangle(d,'Domain',[x1(j,i) y1(j,i); x2(j,i) y2(j,i); x3(j,i) y3(j,i)],'Type','nonproduct');
		% intpoints{j,i,1}=Q.Points;
		% intpoints{j,i,2}=Q.Weights;
		% %--------------------------------------
		% %------- interpolation on those integration points
		% %--------------------------------------
		% intpoints_sigma(j,i,:,1)=F_sigmax(intpoints{j,i,1}(:,1),intpoints{j,i,1}(:,2));
		% intpoints_sigma(j,i,:,2)=F_sigmay(intpoints{j,i,1}(:,1),intpoints{j,i,1}(:,2));
		% intpoints_sigma(j,i,:,3)=F_sigmaxy(intpoints{j,i,1}(:,1),intpoints{j,i,1}(:,2));
		% % for j=1:size(Q.Points,1)
			% % I=I+feval(f,Q.Points(j,1),Q.Points(j,2))*Q.Weights(j);
		% % end
	% end
	% toc
	j=j+1;
end
figure
hold on
for j=1:num_domain
	triplot(rTR{N,j},'b','LineWidth',0.01) 
end
rectangle('Position',[0 0 0.8 0.8],'EdgeColor','k','LineWidth',1)
line([0.5,0.8],[0.4,0.4],'color','k','LineWidth',2)
set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
line([x3(8,8),x2(8,8)],[y3(8,8), y2(8,8)],'color','r','LineWidth',1)
line([x2(8,8),x1(8,8)],[y2(8,8), y1(8,8)],'color','r','LineWidth',1)
line([x1(8,8),x3(8,8)],[y1(8,8), y3(8,8)],'color','r','LineWidth',1)
Q = quadtriangle(d,'Domain',[x1(8,8) y1(8,8); x2(8,8) y2(8,8); x3(8,8) y3(8,8)],'Type','nonproduct');
figure
quadplot(Q,'PlotTitle','off')

	