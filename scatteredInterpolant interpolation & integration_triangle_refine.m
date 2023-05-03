clear all
clc
%--------------------------------------
%------- data from ANSYS
%--------------------------------------
x_coor=xlsread('C:\Users\xs105\Desktop\stress change vs frequency change\2021-2022 second version\manuscript 2022 version\Convergence studies\triangular\stress distribution_120by120.xlsx',1,'A1:A14400');
y_coor=xlsread('C:\Users\xs105\Desktop\stress change vs frequency change\2021-2022 second version\manuscript 2022 version\Convergence studies\triangular\stress distribution_120by120.xlsx',1,'C1:C14400');
sigma_x=xlsread('C:\Users\xs105\Desktop\stress change vs frequency change\2021-2022 second version\manuscript 2022 version\Convergence studies\triangular\stress distribution_120by120.xlsx',1,'E1:E14400');
sigma_y=xlsread('C:\Users\xs105\Desktop\stress change vs frequency change\2021-2022 second version\manuscript 2022 version\Convergence studies\triangular\stress distribution_120by120.xlsx',1,'G1:G14400');
sigma_xy=xlsread('C:\Users\xs105\Desktop\stress change vs frequency change\2021-2022 second version\manuscript 2022 version\Convergence studies\triangular\stress distribution_120by120.xlsx',1,'I1:I14400');
x_coor_sigmax=x_coor;
y_coor_sigmax=y_coor;
x_coor_sigmay=x_coor;
y_coor_sigmay=y_coor;
x_coor_sigmaxy=x_coor;
y_coor_sigmaxy=y_coor;
% % count=1;
% % for i=1:length(sigma_y)
	% % if isnan(sigma_x(i))==1 || isnan(sigma_y(i))==1 || isnan(sigma_xy(i))==1
		% % NAN_num(count,1)=i;
		% % count=count+1;
	% % end
% % end
% % for i=1:length(NAN_num)
	% % sigma_x(NAN_num(length(NAN_num)-i+1),:)=[];
	% % sigma_y(NAN_num(length(NAN_num)-i+1),:)=[];
	% % sigma_xy(NAN_num(length(NAN_num)-i+1),:)=[];
	% % x_coor(NAN_num(length(NAN_num)-i+1),:)=[];
	% % y_coor(NAN_num(length(NAN_num)-i+1),:)=[];
% % end
%%%%%%%%%%%compact sigma_xy results and coordinates
count_xy=1;
for i=1:length(sigma_xy)
	if isnan(sigma_xy(i))==1
		NAN_num_xy(count_xy,1)=i;
		count_xy=count_xy+1;
	end
end
for i=1:length(NAN_num_xy)
	sigma_xy(NAN_num_xy(length(NAN_num_xy)-i+1),:)=[];
	x_coor_sigmaxy(NAN_num_xy(length(NAN_num_xy)-i+1),:)=[];
	y_coor_sigmaxy(NAN_num_xy(length(NAN_num_xy)-i+1),:)=[];
end
%%%%%%%%%%%compact sigma_x results and coordinates
count_x=1;
for i=1:length(sigma_x)
	if isnan(sigma_x(i))==1
		NAN_num_x(count_x,1)=i;
		count_x=count_x+1;
	end
end
for i=1:length(NAN_num_x)
	sigma_x(NAN_num_x(length(NAN_num_x)-i+1),:)=[];
	x_coor_sigmax(NAN_num_x(length(NAN_num_x)-i+1),:)=[];
	y_coor_sigmax(NAN_num_x(length(NAN_num_x)-i+1),:)=[];
end
%%%%%%%%%%%compact sigma_y results and coordinates
count_y=1;
for i=1:length(sigma_y)
	if isnan(sigma_y(i))==1
		NAN_num_y(count_y,1)=i;
		count_y=count_y+1;
	end
end
for i=1:length(NAN_num_y)
	sigma_y(NAN_num_y(length(NAN_num_y)-i+1),:)=[];
	x_coor_sigmay(NAN_num_y(length(NAN_num_y)-i+1),:)=[];
	y_coor_sigmay(NAN_num_y(length(NAN_num_y)-i+1),:)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
scatter3(x_coor_sigmax,y_coor_sigmax,sigma_x)
[xi,yi]=meshgrid(0:0.01:0.8,0:0.01:0.8);
zi_x=griddata(x_coor_sigmax,y_coor_sigmax,sigma_x,xi,yi,'v4');
figure(2)
surf(xi,yi,zi_x);
zi_y=griddata(x_coor_sigmay,y_coor_sigmay,sigma_y,xi,yi,'v4');
figure(3)
surf(xi,yi,zi_y);
zi_xy=griddata(x_coor_sigmaxy,y_coor_sigmaxy,sigma_xy,xi,yi,'v4');
figure(4)
surf(xi,yi,zi_xy);
save('x_coor.mat','x_coor')
save('y_coor.mat','y_coor')
save('sigma_x.mat','sigma_x')
save('sigma_y.mat','sigma_y')
save('sigma_xy.mat','sigma_xy')
save('x_coor_sigmax.mat','x_coor_sigmax')
save('y_coor_sigmax.mat','y_coor_sigmax')
save('x_coor_sigmay.mat','x_coor_sigmay')
save('y_coor_sigmay.mat','y_coor_sigmay')
save('x_coor_sigmaxy.mat','x_coor_sigmaxy')
save('y_coor_sigmaxy.mat','y_coor_sigmaxy')
%--------------------------------------
%------- find interpolant
%--------------------------------------
clear all
clc
% load('x_coor.mat')
% load('y_coor.mat')
load('sigma_x.mat')
load('sigma_y.mat')
load('sigma_xy.mat')
load('x_coor_sigmax.mat')
load('y_coor_sigmax.mat')
load('x_coor_sigmay.mat')
load('y_coor_sigmay.mat')
load('x_coor_sigmaxy.mat')
load('y_coor_sigmaxy.mat')
F_sigmax=scatteredInterpolant(x_coor_sigmax,y_coor_sigmax,sigma_x);
F_sigmax.Method='natural';
F_sigmay=scatteredInterpolant(x_coor_sigmay,y_coor_sigmay,sigma_y);
F_sigmay.Method='natural';
F_sigmaxy=scatteredInterpolant(x_coor_sigmaxy,y_coor_sigmaxy,sigma_xy);
F_sigmaxy.Method='natural';
%--------------------------------------
%------- find integration points and corresponding weights
%--------------------------------------
format long
num_domain=8;						%the number of domains
Points=zeros(3,2,num_domain);
Points(:,:,1)=[0 0.4; 0.48 0.8; 0 0.8];				%from area A1
Points(:,:,2)=[0 0.4; 0.48 0.4; 0.48 0.8];			%from area A2
Points(:,:,3)=[0 0; 0.48 0; 0 0.4];					%from area B1
Points(:,:,4)=[0.48 0; 0.48 0.4; 0 0.4];			%from area B2
Points(:,:,5)=[0.48 0; 0.8 0; 0.8 0.4];				%from area C1
Points(:,:,6)=[0.48 0; 0.8 0.4; 0.48 0.4];			%from area C2
Points(:,:,7)=[0.8 0.4; 0.8 0.8; 0.48 0.8];			%from area D1
Points(:,:,8)=[0.48 0.4; 0.8 0.4; 0.48 0.8];		%from area D2
ConnectivityList=[1 2 3];
N=2;						%the number of refinements    
num_triangles=4^N;			%the number of triangular elements
x1=zeros(num_domain,num_triangles);			%x coordinate of the nodes in triangular elements
y1=zeros(num_domain,num_triangles);			%y coordinate of the nodes in triangular elements
x2=zeros(num_domain,num_triangles);
y2=zeros(num_domain,num_triangles);
x3=zeros(num_domain,num_triangles);
y3=zeros(num_domain,num_triangles);
I=0;								%integration result
d=19;								%degree
num_intpoints=73;					%the number of integartion points in each element
intpoints=cell(num_domain,num_triangles,2);	
for j=1:num_domain
	for i=1:num_triangles
		intpoints{j,i,1}=zeros(num_intpoints,2);		%storing the coordinates of integartion points in each element
		intpoints{j,i,2}=zeros(num_intpoints,1);		%storing the weights of the above points
	end	
end
intpoints_sigma=zeros(num_domain,num_triangles,num_intpoints,3);		%storing interpolated sigma
j=1;
while j<num_domain+1
	TR = triangulation(ConnectivityList,Points(:,:,j));
	figure()
	hold on; 
	box on;
	title(['Original triangulation TR',num2str(j)]);
	triplot(TR,'k','LineWidth',1)  
	[rTR{1:N}] = trirefine(TR,'NumberOfRefinements',N);
	figure()
	hold on; 
	box on;
	title(['refine TR',num2str(j)]);
	triplot(rTR{1,N},'k','LineWidth',1)      
	% for i=1:num_triangles
		% x1(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,1),1);
		% y1(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,1),2);
		% x2(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,2),1);
		% y2(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,2),2);
		% x3(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,3),1);
		% y3(j,i)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(i,3),2);
	% end
	x1(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,1),1);
	y1(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,1),2);
	x2(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,2),1);
	y2(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,2),2);
	x3(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,3),1);
	y3(j,:)=rTR{1,N}.Points(rTR{1,N}.ConnectivityList(:,3),2);
	% f=@(x,y) (2.*x+y).^3;
	tic
	for i=1:num_triangles
		Q = quadtriangle(d,'Domain',[x1(j,i) y1(j,i); x2(j,i) y2(j,i); x3(j,i) y3(j,i)],'Type','nonproduct');
		intpoints{j,i,1}=Q.Points;
		intpoints{j,i,2}=Q.Weights;
		%--------------------------------------
		%------- interpolation on those integration points
		%--------------------------------------
		intpoints_sigma(j,i,:,1)=F_sigmax(intpoints{j,i,1}(:,1),intpoints{j,i,1}(:,2));
		intpoints_sigma(j,i,:,2)=F_sigmay(intpoints{j,i,1}(:,1),intpoints{j,i,1}(:,2));
		intpoints_sigma(j,i,:,3)=F_sigmaxy(intpoints{j,i,1}(:,1),intpoints{j,i,1}(:,2));
		% for j=1:size(Q.Points,1)
			% I=I+feval(f,Q.Points(j,1),Q.Points(j,2))*Q.Weights(j);
		% end
	end
	toc
	j=j+1;
end
% figure()
% scatter3(intpoints{2,1,1}(:,1),intpoints{2,1,1}(:,2),intpoints_sigma(2,1,:,1))
save('intpoints_sigma.mat','intpoints_sigma')
save('intpoints.mat','intpoints')






















