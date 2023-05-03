%SSSS, rectangular, square
clear all
clc
format long
I=4;						%the number of terms in x or y direction
N=2;						%the number of terms of corner functions
a=0.8; 
b=a; 						%a, b are the plate dimensions in x and y directions
h=0.008; 					%thickness
crX=0.4;					%crack length/a
crY=0.5;					%crack location in y direction/b
c=crY*b;					%crack location in y direction
d=crX*a;					%crack length
Alpha=0;					%the angle of crack, in radians
x0=a-d*cos(Alpha);			%x coordinate of crack tip
y0=c-d*sin(Alpha);			%y coordinate of crack tip
nu=0.3; 					%the Poisson's ratio
E=73e9; 					%Young's modulus
Rho=2800; 					%density
D=E*h^3/(12*(1-nu^2)); 		%bending rigidity of the plate
% Ny=-0.1*0.9523*(pi^2*D)/a^2; 						%Ny force (unit: N/m)
y_crackline=@(x) c+(x-a)*tan(Alpha);	%the expression of the extended line of the crack
ll=1;
mm=1;
nn=1;
qq=1;						%indicators of boundary conditions(0: free, 1: simply supported, 2: clamped)
%%%%%%%%%%%%%%%%%%% define admissible functions and corner functions %%%%%%%%%%%%%%%%%%%
syms x y
WBC=x^ll*(x-a)^mm*y^nn*(y-b)^qq;
r=((x-x0)^2+(y-y0)^2)^0.5;
ThetaA=atan((y-y0)/(x-x0))-Alpha;
ThetaB=atan((y-y0)/(x-x0))-Alpha;
ThetaC=pi-Alpha+atan((y-y0)/(x-x0));
ThetaD=-pi-Alpha+atan((y-y0)/(x-x0));
for i=1:I
	PhiX(i)=x^(i-1);
	PhiY(i)=y^(i-1);
end
%corner function:
for n=1:N
	for l=0:n
		i=(3+(n-1))/2*(n-1)+(l+1);
		WSA(i)=r^((2*n+1)/2)*cos((2*l+1)/2*ThetaA);			%如果是（2*n-1）/2，得到的模态看不出discontinuity。有没有可能那种适用于面内的u、v？
		WSB(i)=r^((2*n+1)/2)*cos((2*l+1)/2*ThetaB);
		WSC(i)=r^((2*n+1)/2)*cos((2*l+1)/2*ThetaC);
		WSD(i)=r^((2*n+1)/2)*cos((2*l+1)/2*ThetaD);
		WAA(i)=r^((2*n+1)/2)*sin((2*l+1)/2*ThetaA);
		WAB(i)=r^((2*n+1)/2)*sin((2*l+1)/2*ThetaB);
		WAC(i)=r^((2*n+1)/2)*sin((2*l+1)/2*ThetaC);
		WAD(i)=r^((2*n+1)/2)*sin((2*l+1)/2*ThetaD);
	end
end
tic
%%%%%%%%%%%%%%%%%%% stiffness matrix %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% diagonal matrices %%%%%%%%%%%%%%%%%%%
K1=zeros(I^2,I^2);
K1_1=zeros(I^2,I^2);
K1_2=zeros(I^2,I^2);
K1_3=zeros(I^2,I^2);
K1_4=zeros(I^2,I^2);
K1_5=zeros(I^2,I^2);
K1_6=zeros(I^2,I^2);
K1_7=zeros(I^2,I^2);
Vmax1aa=sym(zeros(I^2,I^2));
Vmax2aa=sym(zeros(I^2,I^2));
Vmax3aa=sym(zeros(I^2,I^2));
Vmax4aa=sym(zeros(I^2,I^2));
Vmax5aa=sym(zeros(I^2,I^2));
Vmax6aa=sym(zeros(I^2,I^2));
Vmax7aa=sym(zeros(I^2,I^2));
func1_1_cell=cell(I^2);
func1_2_cell=cell(I^2);
func1_3_cell=cell(I^2);
func1_4_cell=cell(I^2);
func1_5_cell=cell(I^2);
func1_6_cell=cell(I^2);
func1_7_cell=cell(I^2);
for i=1:I
	for j=1:I
		U(i,j)=j+(i-1)*I;
	end
end
V=U;
for i=1:I^2
	[row1(i),col1(i)]=find(U==i);
	[row2(i),col2(i)]=find(V==i);
end
parfor v=1:I^2
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum strain energy Vmax (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to a(ij) and a(kl)
	Vmax1aa(:,v)=(diff(WBC,x,2).*PhiX(row2(v)).*PhiY(col2(v))+2.*diff(WBC,x,1).*diff(PhiX(row2(v)),1).*PhiY(col2(v))+WBC.*diff(PhiX(row2(v)),2).*PhiY(col2(v))).*(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1));

	Vmax2aa(:,v)=(diff(WBC,y,2).*PhiY(col2(v)).*PhiX(row2(v))+2.*diff(WBC,y,1).*diff(PhiY(col2(v)),1).*PhiX(row2(v))+WBC.*diff(PhiY(col2(v)),2).*PhiX(row2(v))).*(diff(WBC,y,2).*PhiY(col1).*PhiX(row1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1));
				
	Vmax3aa(:,v)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*PhiX(row2(v)).*PhiY(col2(v))+2.*diff(WBC,y,1).*PhiX(row2(v)).*diff(PhiY(col2(v)),1)+WBC.*PhiX(row2(v)).*diff(PhiY(col2(v)),2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*PhiX(row2(v)).*PhiY(col2(v))+2.*diff(WBC,x,1).*diff(PhiX(row2(v)),1).*PhiY(col2(v))+WBC.*diff(PhiX(row2(v)),2).*PhiY(col2(v)));
				
	Vmax4aa(:,v)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*PhiX(row2(v)).*PhiY(col2(v))+diff(WBC,x,1).*PhiX(row2(v)).*diff(PhiY(col2(v)),1)+diff(WBC,y,1).*diff(PhiX(row2(v)),1).*PhiY(col2(v))+WBC.*diff(PhiX(row2(v)),1).*diff(PhiY(col2(v)),1));
				
	Vmax5aa(:,v)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*PhiX(row2(v)).*PhiY(col2(v))+WBC.*PhiX(row2(v)).*diff(PhiY(col2(v)),1));
	
	Vmax6aa(:,v)=(diff(WBC,x,1).*PhiX(row2(v)).*PhiY(col2(v))+WBC.*diff(PhiX(row2(v)),1).*PhiY(col2(v))).*(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1));
	
	Vmax7aa(:,v)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*PhiX(row2(v)).*PhiY(col2(v))+WBC.*PhiX(row2(v)).*diff(PhiY(col2(v)),1))+(diff(WBC,x,1).*PhiX(row2(v)).*PhiY(col2(v))+WBC.*diff(PhiX(row2(v)),1).*PhiY(col2(v))).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
end
num_domain=8;						%the number of domains
refine=2;							%the number of refinements 
num_triangles=4^refine;				%the number of triangular elements in each domain
degree=19;					
num_intpoints=73;					%the number of integartion points in each element
load intpoints.mat
load intpoints_sigma.mat
for v=1:I^2
	for u=1:I^2
		%%%%%%%%%%%%%%%%%%% integration %%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax1aa(u,v)),2)==2
			func1_1_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1aa(u,v))]);
			K1_1(u,v)=D*integral2(func1_1_cell{u,v},0,a,0,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K1_1(u,v)=D*int(int(Vmax1aa(u,v),x,0,a),y,0,b);
		end
		if size(symvar(Vmax2aa(u,v)),2)==2
			func1_2_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2aa(u,v))]);
			K1_2(u,v)=D*integral2(func1_2_cell{u,v},0,a,0,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K1_2(u,v)=D*int(int(Vmax2aa(u,v),x,0,a),y,0,b);
		end
		if size(symvar(Vmax3aa(u,v)),2)==2
			func1_3_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3aa(u,v))]);
			K1_3(u,v)=D*nu*integral2(func1_3_cell{u,v},0,a,0,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K1_3(u,v)=D*nu*int(int(Vmax3aa(u,v),x,0,a),y,0,b);
		end
		if size(symvar(Vmax4aa(u,v)),2)==2
			func1_4_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4aa(u,v))]);
			K1_4(u,v)=D*2*(1-nu)*integral2(func1_4_cell{u,v},0,a,0,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K1_4(u,v)=D*2*(1-nu)*int(int(Vmax4aa(u,v),x,0,a),y,0,b);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5aa(u,v)),2)==2
			func1_5_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5aa(u,v))]);
			for j=1:num_domain
				for i=1:num_triangles
					for q=1:num_intpoints
						K1_5(u,v)=K1_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func1_5_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5aa(u,v)),2)==1
			if strcmp(char(symvar(Vmax5aa(u,v))),'x')==1			%if the symvar is x
				func1_5_cell{u,v}=eval(['@(x)',vectorize(Vmax5aa(u,v))]);
				for j=1:num_domain
					for i=1:num_triangles
						for q=1:num_intpoints
							K1_5(u,v)=K1_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func1_5_cell{u,v},intpoints{j,i,1}(q,1));
						end
					end
				end
			else													%if the symvar is y
				func1_5_cell{u,v}=eval(['@(y)',vectorize(Vmax5aa(u,v))]);
				for j=1:num_domain
					for i=1:num_triangles
						for q=1:num_intpoints
							K1_5(u,v)=K1_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func1_5_cell{u,v},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K1_5(u,v)=K1_5(u,v);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax6aa(u,v)),2)==2
			func1_6_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6aa(u,v))]);
			for j=1:num_domain
				for i=1:num_triangles
					for q=1:num_intpoints
						K1_6(u,v)=K1_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func1_6_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6aa(u,v)),2)==1
			if strcmp(char(symvar(Vmax6aa(u,v))),'x')==1			%if the symvar is x
				func1_6_cell{u,v}=eval(['@(x)',vectorize(Vmax6aa(u,v))]);
				for j=1:num_domain
					for i=1:num_triangles
						for q=1:num_intpoints
							K1_6(u,v)=K1_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func1_6_cell{u,v},intpoints{j,i,1}(q,1));
						end
					end
				end
			else													%if the symvar is y
				func1_6_cell{u,v}=eval(['@(y)',vectorize(Vmax6aa(u,v))]);
				for j=1:num_domain
					for i=1:num_triangles
						for q=1:num_intpoints
							K1_6(u,v)=K1_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func1_6_cell{u,v},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K1_6(u,v)=K1_6(u,v);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax7aa(u,v)),2)==2
			func1_7_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7aa(u,v))]);
			for j=1:num_domain
				for i=1:num_triangles
					for q=1:num_intpoints
						K1_7(u,v)=K1_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func1_7_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7aa(u,v)),2)==1
			if strcmp(char(symvar(Vmax7aa(u,v))),'x')==1			%if the symvar is x
				func1_7_cell{u,v}=eval(['@(x)',vectorize(Vmax7aa(u,v))]);
				for j=1:num_domain
					for i=1:num_triangles
						for q=1:num_intpoints
							K1_7(u,v)=K1_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func1_7_cell{u,v},intpoints{j,i,1}(q,1));
						end
					end
				end
			else													%if the symvar is y
				func1_7_cell{u,v}=eval(['@(y)',vectorize(Vmax7aa(u,v))]);
				for j=1:num_domain
					for i=1:num_triangles
						for q=1:num_intpoints
							K1_7(u,v)=K1_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func1_7_cell{u,v},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K1_7(u,v)=K1_7(u,v);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end
end
K1=K1_1+K1_2+K1_3+K1_4+K1_5+K1_6+K1_7;
toc
tic
K2=zeros(N*(N+3)/2,N*(N+3)/2);
K2_1=zeros(N*(N+3)/2,N*(N+3)/2);
K2_2=zeros(N*(N+3)/2,N*(N+3)/2);
K2_3=zeros(N*(N+3)/2,N*(N+3)/2);
K2_4=zeros(N*(N+3)/2,N*(N+3)/2);
K2_5=zeros(N*(N+3)/2,N*(N+3)/2);
K2_6=zeros(N*(N+3)/2,N*(N+3)/2);
K2_7=zeros(N*(N+3)/2,N*(N+3)/2);
K3=zeros(N*(N+3)/2,N*(N+3)/2);
K3_1=zeros(N*(N+3)/2,N*(N+3)/2);
K3_2=zeros(N*(N+3)/2,N*(N+3)/2);
K3_3=zeros(N*(N+3)/2,N*(N+3)/2);
K3_4=zeros(N*(N+3)/2,N*(N+3)/2);
K3_5=zeros(N*(N+3)/2,N*(N+3)/2);
K3_6=zeros(N*(N+3)/2,N*(N+3)/2);
K3_7=zeros(N*(N+3)/2,N*(N+3)/2);
K=zeros(I^2+N*(N+3),I^2+N*(N+3));
Vmax1bbA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1bbB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1bbC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1bbD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2bbA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2bbB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2bbC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2bbD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3bbA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3bbB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3bbC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3bbD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4bbA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4bbB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4bbC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4bbD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5bbA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5bbB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5bbC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5bbD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6bbA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6bbB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6bbC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6bbD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7bbA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7bbB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7bbC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7bbD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1ccA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1ccB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1ccC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1ccD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2ccA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2ccB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2ccC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2ccD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3ccA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3ccB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3ccC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3ccD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4ccA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4ccB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4ccC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4ccD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5ccA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5ccB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5ccC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5ccD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6ccA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6ccB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6ccC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6ccD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7ccA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7ccB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7ccC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7ccD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
func2_1A_cell=cell(N*(N+3)/2);
func2_1B_cell=cell(N*(N+3)/2);
func2_1C_cell=cell(N*(N+3)/2);
func2_1D_cell=cell(N*(N+3)/2);
func2_2A_cell=cell(N*(N+3)/2);
func2_2B_cell=cell(N*(N+3)/2);
func2_2C_cell=cell(N*(N+3)/2);
func2_2D_cell=cell(N*(N+3)/2);
func2_3A_cell=cell(N*(N+3)/2);
func2_3B_cell=cell(N*(N+3)/2);
func2_3C_cell=cell(N*(N+3)/2);
func2_3D_cell=cell(N*(N+3)/2);
func2_4A_cell=cell(N*(N+3)/2);
func2_4B_cell=cell(N*(N+3)/2);
func2_4C_cell=cell(N*(N+3)/2);
func2_4D_cell=cell(N*(N+3)/2);
func2_5A_cell=cell(N*(N+3)/2);
func2_5B_cell=cell(N*(N+3)/2);
func2_5C_cell=cell(N*(N+3)/2);
func2_5D_cell=cell(N*(N+3)/2);
func2_6A_cell=cell(N*(N+3)/2);
func2_6B_cell=cell(N*(N+3)/2);
func2_6C_cell=cell(N*(N+3)/2);
func2_6D_cell=cell(N*(N+3)/2);
func2_7A_cell=cell(N*(N+3)/2);
func2_7B_cell=cell(N*(N+3)/2);
func2_7C_cell=cell(N*(N+3)/2);
func2_7D_cell=cell(N*(N+3)/2);
func3_1A_cell=cell(N*(N+3)/2);
func3_1B_cell=cell(N*(N+3)/2);
func3_1C_cell=cell(N*(N+3)/2);
func3_1D_cell=cell(N*(N+3)/2);
func3_2A_cell=cell(N*(N+3)/2);
func3_2B_cell=cell(N*(N+3)/2);
func3_2C_cell=cell(N*(N+3)/2);
func3_2D_cell=cell(N*(N+3)/2);
func3_3A_cell=cell(N*(N+3)/2);
func3_3B_cell=cell(N*(N+3)/2);
func3_3C_cell=cell(N*(N+3)/2);
func3_3D_cell=cell(N*(N+3)/2);
func3_4A_cell=cell(N*(N+3)/2);
func3_4B_cell=cell(N*(N+3)/2);
func3_4C_cell=cell(N*(N+3)/2);
func3_4D_cell=cell(N*(N+3)/2);
func3_5A_cell=cell(N*(N+3)/2);
func3_5B_cell=cell(N*(N+3)/2);
func3_5C_cell=cell(N*(N+3)/2);
func3_5D_cell=cell(N*(N+3)/2);
func3_6A_cell=cell(N*(N+3)/2);
func3_6B_cell=cell(N*(N+3)/2);
func3_6C_cell=cell(N*(N+3)/2);
func3_6D_cell=cell(N*(N+3)/2);
func3_7A_cell=cell(N*(N+3)/2);
func3_7B_cell=cell(N*(N+3)/2);
func3_7C_cell=cell(N*(N+3)/2);
func3_7D_cell=cell(N*(N+3)/2);
for i=1:N*(N+3)/2
	row_num(i)=i;
end
parfor j=1:N*(N+3)/2
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum strain energy Vmax (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to b(i) and b(j)
	Vmax1bbA(:,j)=(diff(WBC,x,2).*WSA(row_num)+2.*diff(WBC,x,1).*diff(WSA(row_num),x)+WBC.*diff(WSA(row_num),x,2)).*(diff(WBC,x,2).*WSA(j)+2.*diff(WBC,x,1).*diff(WSA(j),x)+WBC.*diff(WSA(j),x,2));
	Vmax1bbB(:,j)=(diff(WBC,x,2).*WSB(row_num)+2.*diff(WBC,x,1).*diff(WSB(row_num),x)+WBC.*diff(WSB(row_num),x,2)).*(diff(WBC,x,2).*WSB(j)+2.*diff(WBC,x,1).*diff(WSB(j),x)+WBC.*diff(WSB(j),x,2));
	Vmax1bbC(:,j)=(diff(WBC,x,2).*WSC(row_num)+2.*diff(WBC,x,1).*diff(WSC(row_num),x)+WBC.*diff(WSC(row_num),x,2)).*(diff(WBC,x,2).*WSC(j)+2.*diff(WBC,x,1).*diff(WSC(j),x)+WBC.*diff(WSC(j),x,2));
	Vmax1bbD(:,j)=(diff(WBC,x,2).*WSD(row_num)+2.*diff(WBC,x,1).*diff(WSD(row_num),x)+WBC.*diff(WSD(row_num),x,2)).*(diff(WBC,x,2).*WSD(j)+2.*diff(WBC,x,1).*diff(WSD(j),x)+WBC.*diff(WSD(j),x,2));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax2bbA(:,j)=(diff(WBC,y,2).*WSA(row_num)+2.*diff(WBC,y,1).*diff(WSA(row_num),y)+WBC.*diff(WSA(row_num),y,2)).*(diff(WBC,y,2).*WSA(j)+2.*diff(WBC,y,1).*diff(WSA(j),y)+WBC.*diff(WSA(j),y,2));
	Vmax2bbB(:,j)=(diff(WBC,y,2).*WSB(row_num)+2.*diff(WBC,y,1).*diff(WSB(row_num),y)+WBC.*diff(WSB(row_num),y,2)).*(diff(WBC,y,2).*WSB(j)+2.*diff(WBC,y,1).*diff(WSB(j),y)+WBC.*diff(WSB(j),y,2));
	Vmax2bbC(:,j)=(diff(WBC,y,2).*WSC(row_num)+2.*diff(WBC,y,1).*diff(WSC(row_num),y)+WBC.*diff(WSC(row_num),y,2)).*(diff(WBC,y,2).*WSC(j)+2.*diff(WBC,y,1).*diff(WSC(j),y)+WBC.*diff(WSC(j),y,2));
	Vmax2bbD(:,j)=(diff(WBC,y,2).*WSD(row_num)+2.*diff(WBC,y,1).*diff(WSD(row_num),y)+WBC.*diff(WSD(row_num),y,2)).*(diff(WBC,y,2).*WSD(j)+2.*diff(WBC,y,1).*diff(WSD(j),y)+WBC.*diff(WSD(j),y,2));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax3bbA(:,j)=(diff(WBC,x,2).*WSA(row_num)+2.*diff(WBC,x).*diff(WSA(row_num),x)+WBC.*diff(WSA(row_num),x,2)).*(diff(WBC,y,2).*WSA(j)+2.*diff(WBC,y).*diff(WSA(j),y)+WBC.*diff(WSA(j),y,2))+(diff(WBC,y,2).*WSA(row_num)+2.*diff(WBC,y).*diff(WSA(row_num),y)+WBC.*diff(WSA(row_num),y,2)).*(diff(WBC,x,2).*WSA(j)+2.*diff(WBC,x).*diff(WSA(j),x)+WBC.*diff(WSA(j),x,2));
	Vmax3bbB(:,j)=(diff(WBC,x,2).*WSB(row_num)+2.*diff(WBC,x).*diff(WSB(row_num),x)+WBC.*diff(WSB(row_num),x,2)).*(diff(WBC,y,2).*WSB(j)+2.*diff(WBC,y).*diff(WSB(j),y)+WBC.*diff(WSB(j),y,2))+(diff(WBC,y,2).*WSB(row_num)+2.*diff(WBC,y).*diff(WSB(row_num),y)+WBC.*diff(WSB(row_num),y,2)).*(diff(WBC,x,2).*WSB(j)+2.*diff(WBC,x).*diff(WSB(j),x)+WBC.*diff(WSB(j),x,2));
	Vmax3bbC(:,j)=(diff(WBC,x,2).*WSC(row_num)+2.*diff(WBC,x).*diff(WSC(row_num),x)+WBC.*diff(WSC(row_num),x,2)).*(diff(WBC,y,2).*WSC(j)+2.*diff(WBC,y).*diff(WSC(j),y)+WBC.*diff(WSC(j),y,2))+(diff(WBC,y,2).*WSC(row_num)+2.*diff(WBC,y).*diff(WSC(row_num),y)+WBC.*diff(WSC(row_num),y,2)).*(diff(WBC,x,2).*WSC(j)+2.*diff(WBC,x).*diff(WSC(j),x)+WBC.*diff(WSC(j),x,2));
	Vmax3bbD(:,j)=(diff(WBC,x,2).*WSD(row_num)+2.*diff(WBC,x).*diff(WSD(row_num),x)+WBC.*diff(WSD(row_num),x,2)).*(diff(WBC,y,2).*WSD(j)+2.*diff(WBC,y).*diff(WSD(j),y)+WBC.*diff(WSD(j),y,2))+(diff(WBC,y,2).*WSD(row_num)+2.*diff(WBC,y).*diff(WSD(row_num),y)+WBC.*diff(WSD(row_num),y,2)).*(diff(WBC,x,2).*WSD(j)+2.*diff(WBC,x).*diff(WSD(j),x)+WBC.*diff(WSD(j),x,2));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax4bbA(:,j)=(diff(WBC,x,y).*WSA(row_num)+diff(WBC,x).*diff(WSA(row_num),y)+diff(WBC,y).*diff(WSA(row_num),x)+WBC.*diff(WSA(row_num),x,y)).*(diff(WBC,x,y).*WSA(j)+diff(WBC,x).*diff(WSA(j),y)+diff(WBC,y).*diff(WSA(j),x)+WBC.*diff(WSA(j),x,y));
	Vmax4bbB(:,j)=(diff(WBC,x,y).*WSB(row_num)+diff(WBC,x).*diff(WSB(row_num),y)+diff(WBC,y).*diff(WSB(row_num),x)+WBC.*diff(WSB(row_num),x,y)).*(diff(WBC,x,y).*WSB(j)+diff(WBC,x).*diff(WSB(j),y)+diff(WBC,y).*diff(WSB(j),x)+WBC.*diff(WSB(j),x,y));
	Vmax4bbC(:,j)=(diff(WBC,x,y).*WSC(row_num)+diff(WBC,x).*diff(WSC(row_num),y)+diff(WBC,y).*diff(WSC(row_num),x)+WBC.*diff(WSC(row_num),x,y)).*(diff(WBC,x,y).*WSC(j)+diff(WBC,x).*diff(WSC(j),y)+diff(WBC,y).*diff(WSC(j),x)+WBC.*diff(WSC(j),x,y));
	Vmax4bbD(:,j)=(diff(WBC,x,y).*WSD(row_num)+diff(WBC,x).*diff(WSD(row_num),y)+diff(WBC,y).*diff(WSD(row_num),x)+WBC.*diff(WSD(row_num),x,y)).*(diff(WBC,x,y).*WSD(j)+diff(WBC,x).*diff(WSD(j),y)+diff(WBC,y).*diff(WSD(j),x)+WBC.*diff(WSD(j),x,y));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax5bbA(:,j)=(diff(WBC,y).*WSA(row_num)+WBC.*diff(WSA(row_num),y)).*(diff(WBC,y).*WSA(j)+WBC.*diff(WSA(j),y));
	Vmax5bbB(:,j)=(diff(WBC,y).*WSB(row_num)+WBC.*diff(WSB(row_num),y)).*(diff(WBC,y).*WSB(j)+WBC.*diff(WSB(j),y));
	Vmax5bbC(:,j)=(diff(WBC,y).*WSC(row_num)+WBC.*diff(WSC(row_num),y)).*(diff(WBC,y).*WSC(j)+WBC.*diff(WSC(j),y));
	Vmax5bbD(:,j)=(diff(WBC,y).*WSD(row_num)+WBC.*diff(WSD(row_num),y)).*(diff(WBC,y).*WSD(j)+WBC.*diff(WSD(j),y));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax6bbA(:,j)=(diff(WBC,x).*WSA(row_num)+WBC.*diff(WSA(row_num),x)).*(diff(WBC,x).*WSA(j)+WBC.*diff(WSA(j),x));
	Vmax6bbB(:,j)=(diff(WBC,x).*WSB(row_num)+WBC.*diff(WSB(row_num),x)).*(diff(WBC,x).*WSB(j)+WBC.*diff(WSB(j),x));
	Vmax6bbC(:,j)=(diff(WBC,x).*WSC(row_num)+WBC.*diff(WSC(row_num),x)).*(diff(WBC,x).*WSC(j)+WBC.*diff(WSC(j),x));
	Vmax6bbD(:,j)=(diff(WBC,x).*WSD(row_num)+WBC.*diff(WSD(row_num),x)).*(diff(WBC,x).*WSD(j)+WBC.*diff(WSD(j),x));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax7bbA(:,j)=(diff(WBC,x).*WSA(row_num)+WBC.*diff(WSA(row_num),x)).*(diff(WBC,y).*WSA(j)+WBC.*diff(WSA(j),y))+(diff(WBC,x).*WSA(j)+WBC.*diff(WSA(j),x)).*(diff(WBC,y).*WSA(row_num)+WBC.*diff(WSA(row_num),y));
	Vmax7bbB(:,j)=(diff(WBC,x).*WSB(row_num)+WBC.*diff(WSB(row_num),x)).*(diff(WBC,y).*WSB(j)+WBC.*diff(WSB(j),y))+(diff(WBC,x).*WSB(j)+WBC.*diff(WSB(j),x)).*(diff(WBC,y).*WSB(row_num)+WBC.*diff(WSB(row_num),y));
	Vmax7bbC(:,j)=(diff(WBC,x).*WSC(row_num)+WBC.*diff(WSC(row_num),x)).*(diff(WBC,y).*WSC(j)+WBC.*diff(WSC(j),y))+(diff(WBC,x).*WSC(j)+WBC.*diff(WSC(j),x)).*(diff(WBC,y).*WSC(row_num)+WBC.*diff(WSC(row_num),y));
	Vmax7bbD(:,j)=(diff(WBC,x).*WSD(row_num)+WBC.*diff(WSD(row_num),x)).*(diff(WBC,y).*WSD(j)+WBC.*diff(WSD(j),y))+(diff(WBC,x).*WSD(j)+WBC.*diff(WSD(j),x)).*(diff(WBC,y).*WSD(row_num)+WBC.*diff(WSD(row_num),y));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum strain energy Vmax (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to c(i) and c(j)
	Vmax1ccA(:,j)=(diff(WBC,x,2).*WAA(row_num)+2.*diff(WBC,x,1).*diff(WAA(row_num),x)+WBC.*diff(WAA(row_num),x,2)).*(diff(WBC,x,2).*WAA(j)+2.*diff(WBC,x,1).*diff(WAA(j),x)+WBC.*diff(WAA(j),x,2));
	Vmax1ccB(:,j)=(diff(WBC,x,2).*WAB(row_num)+2.*diff(WBC,x,1).*diff(WAB(row_num),x)+WBC.*diff(WAB(row_num),x,2)).*(diff(WBC,x,2).*WAB(j)+2.*diff(WBC,x,1).*diff(WAB(j),x)+WBC.*diff(WAB(j),x,2));
	Vmax1ccC(:,j)=(diff(WBC,x,2).*WAC(row_num)+2.*diff(WBC,x,1).*diff(WAC(row_num),x)+WBC.*diff(WAC(row_num),x,2)).*(diff(WBC,x,2).*WAC(j)+2.*diff(WBC,x,1).*diff(WAC(j),x)+WBC.*diff(WAC(j),x,2));
	Vmax1ccD(:,j)=(diff(WBC,x,2).*WAD(row_num)+2.*diff(WBC,x,1).*diff(WAD(row_num),x)+WBC.*diff(WAD(row_num),x,2)).*(diff(WBC,x,2).*WAD(j)+2.*diff(WBC,x,1).*diff(WAD(j),x)+WBC.*diff(WAD(j),x,2));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax2ccA(:,j)=(diff(WBC,y,2).*WAA(row_num)+2.*diff(WBC,y,1).*diff(WAA(row_num),y)+WBC.*diff(WAA(row_num),y,2)).*(diff(WBC,y,2).*WAA(j)+2.*diff(WBC,y,1).*diff(WAA(j),y)+WBC.*diff(WAA(j),y,2));
	Vmax2ccB(:,j)=(diff(WBC,y,2).*WAB(row_num)+2.*diff(WBC,y,1).*diff(WAB(row_num),y)+WBC.*diff(WAB(row_num),y,2)).*(diff(WBC,y,2).*WAB(j)+2.*diff(WBC,y,1).*diff(WAB(j),y)+WBC.*diff(WAB(j),y,2));
	Vmax2ccC(:,j)=(diff(WBC,y,2).*WAC(row_num)+2.*diff(WBC,y,1).*diff(WAC(row_num),y)+WBC.*diff(WAC(row_num),y,2)).*(diff(WBC,y,2).*WAC(j)+2.*diff(WBC,y,1).*diff(WAC(j),y)+WBC.*diff(WAC(j),y,2));
	Vmax2ccD(:,j)=(diff(WBC,y,2).*WAD(row_num)+2.*diff(WBC,y,1).*diff(WAD(row_num),y)+WBC.*diff(WAD(row_num),y,2)).*(diff(WBC,y,2).*WAD(j)+2.*diff(WBC,y,1).*diff(WAD(j),y)+WBC.*diff(WAD(j),y,2));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax3ccA(:,j)=(diff(WBC,x,2).*WAA(row_num)+2.*diff(WBC,x).*diff(WAA(row_num),x)+WBC.*diff(WAA(row_num),x,2)).*(diff(WBC,y,2).*WAA(j)+2.*diff(WBC,y).*diff(WAA(j),y)+WBC.*diff(WAA(j),y,2))+(diff(WBC,y,2).*WAA(row_num)+2.*diff(WBC,y).*diff(WAA(row_num),y)+WBC.*diff(WAA(row_num),y,2)).*(diff(WBC,x,2).*WAA(j)+2.*diff(WBC,x).*diff(WAA(j),x)+WBC.*diff(WAA(j),x,2));
	Vmax3ccB(:,j)=(diff(WBC,x,2).*WAB(row_num)+2.*diff(WBC,x).*diff(WAB(row_num),x)+WBC.*diff(WAB(row_num),x,2)).*(diff(WBC,y,2).*WAB(j)+2.*diff(WBC,y).*diff(WAB(j),y)+WBC.*diff(WAB(j),y,2))+(diff(WBC,y,2).*WAB(row_num)+2.*diff(WBC,y).*diff(WAB(row_num),y)+WBC.*diff(WAB(row_num),y,2)).*(diff(WBC,x,2).*WAB(j)+2.*diff(WBC,x).*diff(WAB(j),x)+WBC.*diff(WAB(j),x,2));
	Vmax3ccC(:,j)=(diff(WBC,x,2).*WAC(row_num)+2.*diff(WBC,x).*diff(WAC(row_num),x)+WBC.*diff(WAC(row_num),x,2)).*(diff(WBC,y,2).*WAC(j)+2.*diff(WBC,y).*diff(WAC(j),y)+WBC.*diff(WAC(j),y,2))+(diff(WBC,y,2).*WAC(row_num)+2.*diff(WBC,y).*diff(WAC(row_num),y)+WBC.*diff(WAC(row_num),y,2)).*(diff(WBC,x,2).*WAC(j)+2.*diff(WBC,x).*diff(WAC(j),x)+WBC.*diff(WAC(j),x,2));
	Vmax3ccD(:,j)=(diff(WBC,x,2).*WAD(row_num)+2.*diff(WBC,x).*diff(WAD(row_num),x)+WBC.*diff(WAD(row_num),x,2)).*(diff(WBC,y,2).*WAD(j)+2.*diff(WBC,y).*diff(WAD(j),y)+WBC.*diff(WAD(j),y,2))+(diff(WBC,y,2).*WAD(row_num)+2.*diff(WBC,y).*diff(WAD(row_num),y)+WBC.*diff(WAD(row_num),y,2)).*(diff(WBC,x,2).*WAD(j)+2.*diff(WBC,x).*diff(WAD(j),x)+WBC.*diff(WAD(j),x,2));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax4ccA(:,j)=(diff(WBC,x,y).*WAA(row_num)+diff(WBC,x).*diff(WAA(row_num),y)+diff(WBC,y).*diff(WAA(row_num),x)+WBC.*diff(WAA(row_num),x,y)).*(diff(WBC,x,y).*WAA(j)+diff(WBC,x).*diff(WAA(j),y)+diff(WBC,y).*diff(WAA(j),x)+WBC.*diff(WAA(j),x,y));
	Vmax4ccB(:,j)=(diff(WBC,x,y).*WAB(row_num)+diff(WBC,x).*diff(WAB(row_num),y)+diff(WBC,y).*diff(WAB(row_num),x)+WBC.*diff(WAB(row_num),x,y)).*(diff(WBC,x,y).*WAB(j)+diff(WBC,x).*diff(WAB(j),y)+diff(WBC,y).*diff(WAB(j),x)+WBC.*diff(WAB(j),x,y));
	Vmax4ccC(:,j)=(diff(WBC,x,y).*WAC(row_num)+diff(WBC,x).*diff(WAC(row_num),y)+diff(WBC,y).*diff(WAC(row_num),x)+WBC.*diff(WAC(row_num),x,y)).*(diff(WBC,x,y).*WAC(j)+diff(WBC,x).*diff(WAC(j),y)+diff(WBC,y).*diff(WAC(j),x)+WBC.*diff(WAC(j),x,y));
	Vmax4ccD(:,j)=(diff(WBC,x,y).*WAD(row_num)+diff(WBC,x).*diff(WAD(row_num),y)+diff(WBC,y).*diff(WAD(row_num),x)+WBC.*diff(WAD(row_num),x,y)).*(diff(WBC,x,y).*WAD(j)+diff(WBC,x).*diff(WAD(j),y)+diff(WBC,y).*diff(WAD(j),x)+WBC.*diff(WAD(j),x,y));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax5ccA(:,j)=(diff(WBC,y).*WAA(row_num)+WBC.*diff(WAA(row_num),y)).*(diff(WBC,y).*WAA(j)+WBC.*diff(WAA(j),y));
	Vmax5ccB(:,j)=(diff(WBC,y).*WAB(row_num)+WBC.*diff(WAB(row_num),y)).*(diff(WBC,y).*WAB(j)+WBC.*diff(WAB(j),y));
	Vmax5ccC(:,j)=(diff(WBC,y).*WAC(row_num)+WBC.*diff(WAC(row_num),y)).*(diff(WBC,y).*WAC(j)+WBC.*diff(WAC(j),y));
	Vmax5ccD(:,j)=(diff(WBC,y).*WAD(row_num)+WBC.*diff(WAD(row_num),y)).*(diff(WBC,y).*WAD(j)+WBC.*diff(WAD(j),y));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax6ccA(:,j)=(diff(WBC,x).*WAA(row_num)+WBC.*diff(WAA(row_num),x)).*(diff(WBC,x).*WAA(j)+WBC.*diff(WAA(j),x));
	Vmax6ccB(:,j)=(diff(WBC,x).*WAB(row_num)+WBC.*diff(WAB(row_num),x)).*(diff(WBC,x).*WAB(j)+WBC.*diff(WAB(j),x));
	Vmax6ccC(:,j)=(diff(WBC,x).*WAC(row_num)+WBC.*diff(WAC(row_num),x)).*(diff(WBC,x).*WAC(j)+WBC.*diff(WAC(j),x));
	Vmax6ccD(:,j)=(diff(WBC,x).*WAD(row_num)+WBC.*diff(WAD(row_num),x)).*(diff(WBC,x).*WAD(j)+WBC.*diff(WAD(j),x));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax7ccA(:,j)=(diff(WBC,x).*WAA(row_num)+WBC.*diff(WAA(row_num),x)).*(diff(WBC,y).*WAA(j)+WBC.*diff(WAA(j),y))+(diff(WBC,x).*WAA(j)+WBC.*diff(WAA(j),x)).*(diff(WBC,y).*WAA(row_num)+WBC.*diff(WAA(row_num),y));
	Vmax7ccB(:,j)=(diff(WBC,x).*WAB(row_num)+WBC.*diff(WAB(row_num),x)).*(diff(WBC,y).*WAB(j)+WBC.*diff(WAB(j),y))+(diff(WBC,x).*WAB(j)+WBC.*diff(WAB(j),x)).*(diff(WBC,y).*WAB(row_num)+WBC.*diff(WAB(row_num),y));
	Vmax7ccC(:,j)=(diff(WBC,x).*WAC(row_num)+WBC.*diff(WAC(row_num),x)).*(diff(WBC,y).*WAC(j)+WBC.*diff(WAC(j),y))+(diff(WBC,x).*WAC(j)+WBC.*diff(WAC(j),x)).*(diff(WBC,y).*WAC(row_num)+WBC.*diff(WAC(row_num),y));
	Vmax7ccD(:,j)=(diff(WBC,x).*WAD(row_num)+WBC.*diff(WAD(row_num),x)).*(diff(WBC,y).*WAD(j)+WBC.*diff(WAD(j),y))+(diff(WBC,x).*WAD(j)+WBC.*diff(WAD(j),x)).*(diff(WBC,y).*WAD(row_num)+WBC.*diff(WAD(row_num),y));
end
for v=1:N*(N+3)/2
	for u=1:N*(N+3)/2
		%%%%%%%%%%%%%%%%%%% integration %%%%%%%%%%%%%%%%%%%
		func2_1A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1bbA(u,v))]);
		func2_1B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1bbB(u,v))]);
		func2_1C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1bbC(u,v))]);
		func2_1D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1bbD(u,v))]);
		K2_1(u,v)=D*(integral2(func2_1A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_1B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_1C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_1D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func2_2A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2bbA(u,v))]);
		func2_2B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2bbB(u,v))]);
		func2_2C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2bbC(u,v))]);
		func2_2D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2bbD(u,v))]);
		K2_2(u,v)=D*(integral2(func2_2A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_2B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_2C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_2D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func2_3A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3bbA(u,v))]);
		func2_3B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3bbB(u,v))]);
		func2_3C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3bbC(u,v))]);
		func2_3D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3bbD(u,v))]);
		K2_3(u,v)=D*nu*(integral2(func2_3A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_3B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_3C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_3D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func2_4A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4bbA(u,v))]);
		func2_4B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4bbB(u,v))]);
		func2_4C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4bbC(u,v))]);
		func2_4D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4bbD(u,v))]);
		K2_4(u,v)=D*2*(1-nu)*(integral2(func2_4A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_4B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_4C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func2_4D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func2_5A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5bbA(u,v))]);
		func2_6A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6bbA(u,v))]);
		func2_7A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7bbA(u,v))]);
		func3_5A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5ccA(u,v))]);
		func3_6A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6ccA(u,v))]);
		func3_7A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7ccA(u,v))]);
		func2_5B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5bbB(u,v))]);
		func2_6B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6bbB(u,v))]);
		func2_7B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7bbB(u,v))]);
		func3_5B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5ccB(u,v))]);
		func3_6B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6ccB(u,v))]);
		func3_7B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7ccB(u,v))]);
		func2_5C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5bbC(u,v))]);
		func2_6C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6bbC(u,v))]);
		func2_7C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7bbC(u,v))]);
		func3_5C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5ccC(u,v))]);
		func3_6C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6ccC(u,v))]);
		func3_7C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7ccC(u,v))]);
		func2_5D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5bbD(u,v))]);
		func2_6D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6bbD(u,v))]);
		func2_7D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7bbD(u,v))]);
		func3_5D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5ccD(u,v))]);
		func3_6D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6ccD(u,v))]);
		func3_7D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7ccD(u,v))]);
		for j=1:num_domain
			if j==1 || j==2								%for area A1 and A2
				for i=1:num_triangles
					for q=1:num_intpoints
						K2_5(u,v)=K2_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func2_5A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K2_6(u,v)=K2_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func2_6A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K2_7(u,v)=K2_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func2_7A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_5(u,v)=K3_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func3_5A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_6(u,v)=K3_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func3_6A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_7(u,v)=K3_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func3_7A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			elseif j==3 || j==4							%for area B1 and B2
				for i=1:num_triangles
					for q=1:num_intpoints
						K2_5(u,v)=K2_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func2_5B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K2_6(u,v)=K2_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func2_6B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K2_7(u,v)=K2_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func2_7B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_5(u,v)=K3_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func3_5B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_6(u,v)=K3_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func3_6B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_7(u,v)=K3_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func3_7B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			elseif j==5 || j==6							%for area C1 and C2
				for i=1:num_triangles
					for q=1:num_intpoints
						K2_5(u,v)=K2_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func2_5C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K2_6(u,v)=K2_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func2_6C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K2_7(u,v)=K2_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func2_7C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_5(u,v)=K3_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func3_5C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_6(u,v)=K3_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func3_6C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_7(u,v)=K3_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func3_7C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			else										%for area D1 and D2
				for i=1:num_triangles
					for q=1:num_intpoints
						K2_5(u,v)=K2_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func2_5D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K2_6(u,v)=K2_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func2_6D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K2_7(u,v)=K2_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func2_7D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_5(u,v)=K3_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func3_5D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_6(u,v)=K3_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func3_6D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K3_7(u,v)=K3_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func3_7D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func3_1A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1ccA(u,v))]);
		func3_1B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1ccB(u,v))]);
		func3_1C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1ccC(u,v))]);
		func3_1D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1ccD(u,v))]);
		K3_1(u,v)=D*(integral2(func3_1A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_1B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_1C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_1D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func3_2A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2ccA(u,v))]);
		func3_2B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2ccB(u,v))]);
		func3_2C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2ccC(u,v))]);
		func3_2D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2ccD(u,v))]);
		K3_2(u,v)=D*(integral2(func3_2A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_2B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_2C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_2D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func3_3A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3ccA(u,v))]);
		func3_3B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3ccB(u,v))]);
		func3_3C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3ccC(u,v))]);
		func3_3D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3ccD(u,v))]);
		K3_3(u,v)=D*nu*(integral2(func3_3A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_3B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_3C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_3D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func3_4A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4ccA(u,v))]);
		func3_4B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4ccB(u,v))]);
		func3_4C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4ccC(u,v))]);
		func3_4D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4ccD(u,v))]);
		K3_4(u,v)=D*2*(1-nu)*(integral2(func3_4A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_4B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_4C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func3_4D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
	end
end
K2=K2_1+K2_2+K2_3+K2_4+K2_5+K2_6+K2_7;
K3=K3_1+K3_2+K3_3+K3_4+K3_5+K3_6+K3_7;
K(1:I^2,1:I^2)=K1;
K(I^2+1:I^2+N*(N+3)/2,I^2+1:I^2+N*(N+3)/2)=K2;
K(I^2+N*(N+3)/2+1:I^2+N*(N+3),I^2+N*(N+3)/2+1:I^2+N*(N+3))=K3;
toc
tic
%%%%%%%%%%%%%%%%%%% stiffness matrix %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% off-diagonal matrices %%%%%%%%%%%%%%%%%%%
K12=zeros(I^2,N*(N+3)/2);
K12_1=zeros(I^2,N*(N+3)/2);
K12_1A=zeros(I^2,N*(N+3)/2);
K12_1B=zeros(I^2,N*(N+3)/2);
K12_1C=zeros(I^2,N*(N+3)/2);
K12_1D=zeros(I^2,N*(N+3)/2);
K12_2=zeros(I^2,N*(N+3)/2);
K12_2A=zeros(I^2,N*(N+3)/2);
K12_2B=zeros(I^2,N*(N+3)/2);
K12_2C=zeros(I^2,N*(N+3)/2);
K12_2D=zeros(I^2,N*(N+3)/2);
K12_3=zeros(I^2,N*(N+3)/2);
K12_3A=zeros(I^2,N*(N+3)/2);
K12_3B=zeros(I^2,N*(N+3)/2);
K12_3C=zeros(I^2,N*(N+3)/2);
K12_3D=zeros(I^2,N*(N+3)/2);
K12_4=zeros(I^2,N*(N+3)/2);
K12_4A=zeros(I^2,N*(N+3)/2);
K12_4B=zeros(I^2,N*(N+3)/2);
K12_4C=zeros(I^2,N*(N+3)/2);
K12_4D=zeros(I^2,N*(N+3)/2);
K12_5=zeros(I^2,N*(N+3)/2);
K12_5A=zeros(I^2,N*(N+3)/2);
K12_5B=zeros(I^2,N*(N+3)/2);
K12_5C=zeros(I^2,N*(N+3)/2);
K12_5D=zeros(I^2,N*(N+3)/2);
K12_6=zeros(I^2,N*(N+3)/2);
K12_6A=zeros(I^2,N*(N+3)/2);
K12_6B=zeros(I^2,N*(N+3)/2);
K12_6C=zeros(I^2,N*(N+3)/2);
K12_6D=zeros(I^2,N*(N+3)/2);
K12_7=zeros(I^2,N*(N+3)/2);
K12_7A=zeros(I^2,N*(N+3)/2);
K12_7B=zeros(I^2,N*(N+3)/2);
K12_7C=zeros(I^2,N*(N+3)/2);
K12_7D=zeros(I^2,N*(N+3)/2);
K13=zeros(I^2,N*(N+3)/2);
K13_1=zeros(I^2,N*(N+3)/2);
K13_1A=zeros(I^2,N*(N+3)/2);
K13_1B=zeros(I^2,N*(N+3)/2);
K13_1C=zeros(I^2,N*(N+3)/2);
K13_1D=zeros(I^2,N*(N+3)/2);
K13_2=zeros(I^2,N*(N+3)/2);
K13_2A=zeros(I^2,N*(N+3)/2);
K13_2B=zeros(I^2,N*(N+3)/2);
K13_2C=zeros(I^2,N*(N+3)/2);
K13_2D=zeros(I^2,N*(N+3)/2);
K13_3=zeros(I^2,N*(N+3)/2);
K13_3A=zeros(I^2,N*(N+3)/2);
K13_3B=zeros(I^2,N*(N+3)/2);
K13_3C=zeros(I^2,N*(N+3)/2);
K13_3D=zeros(I^2,N*(N+3)/2);
K13_4=zeros(I^2,N*(N+3)/2);
K13_4A=zeros(I^2,N*(N+3)/2);
K13_4B=zeros(I^2,N*(N+3)/2);
K13_4C=zeros(I^2,N*(N+3)/2);
K13_4D=zeros(I^2,N*(N+3)/2);
K13_5=zeros(I^2,N*(N+3)/2);
K13_5A=zeros(I^2,N*(N+3)/2);
K13_5B=zeros(I^2,N*(N+3)/2);
K13_5C=zeros(I^2,N*(N+3)/2);
K13_5D=zeros(I^2,N*(N+3)/2);
K13_6=zeros(I^2,N*(N+3)/2);
K13_6A=zeros(I^2,N*(N+3)/2);
K13_6B=zeros(I^2,N*(N+3)/2);
K13_6C=zeros(I^2,N*(N+3)/2);
K13_6D=zeros(I^2,N*(N+3)/2);
K13_7=zeros(I^2,N*(N+3)/2);
K13_7A=zeros(I^2,N*(N+3)/2);
K13_7B=zeros(I^2,N*(N+3)/2);
K13_7C=zeros(I^2,N*(N+3)/2);
K13_7D=zeros(I^2,N*(N+3)/2);
Vmax1abA=sym(zeros(I^2,N*(N+3)/2));
Vmax1abB=sym(zeros(I^2,N*(N+3)/2));
Vmax1abC=sym(zeros(I^2,N*(N+3)/2));
Vmax1abD=sym(zeros(I^2,N*(N+3)/2));
Vmax2abA=sym(zeros(I^2,N*(N+3)/2));
Vmax2abB=sym(zeros(I^2,N*(N+3)/2));
Vmax2abC=sym(zeros(I^2,N*(N+3)/2));
Vmax2abD=sym(zeros(I^2,N*(N+3)/2));
Vmax3abA=sym(zeros(I^2,N*(N+3)/2));
Vmax3abB=sym(zeros(I^2,N*(N+3)/2));
Vmax3abC=sym(zeros(I^2,N*(N+3)/2));
Vmax3abD=sym(zeros(I^2,N*(N+3)/2));
Vmax4abA=sym(zeros(I^2,N*(N+3)/2));
Vmax4abB=sym(zeros(I^2,N*(N+3)/2));
Vmax4abC=sym(zeros(I^2,N*(N+3)/2));
Vmax4abD=sym(zeros(I^2,N*(N+3)/2));
Vmax5abA=sym(zeros(I^2,N*(N+3)/2));
Vmax5abB=sym(zeros(I^2,N*(N+3)/2));
Vmax5abC=sym(zeros(I^2,N*(N+3)/2));
Vmax5abD=sym(zeros(I^2,N*(N+3)/2));
Vmax6abA=sym(zeros(I^2,N*(N+3)/2));
Vmax6abB=sym(zeros(I^2,N*(N+3)/2));
Vmax6abC=sym(zeros(I^2,N*(N+3)/2));
Vmax6abD=sym(zeros(I^2,N*(N+3)/2));
Vmax7abA=sym(zeros(I^2,N*(N+3)/2));
Vmax7abB=sym(zeros(I^2,N*(N+3)/2));
Vmax7abC=sym(zeros(I^2,N*(N+3)/2));
Vmax7abD=sym(zeros(I^2,N*(N+3)/2));
Vmax1acA=sym(zeros(I^2,N*(N+3)/2));
Vmax1acB=sym(zeros(I^2,N*(N+3)/2));
Vmax1acC=sym(zeros(I^2,N*(N+3)/2));
Vmax1acD=sym(zeros(I^2,N*(N+3)/2));
Vmax2acA=sym(zeros(I^2,N*(N+3)/2));
Vmax2acB=sym(zeros(I^2,N*(N+3)/2));
Vmax2acC=sym(zeros(I^2,N*(N+3)/2));
Vmax2acD=sym(zeros(I^2,N*(N+3)/2));
Vmax3acA=sym(zeros(I^2,N*(N+3)/2));
Vmax3acB=sym(zeros(I^2,N*(N+3)/2));
Vmax3acC=sym(zeros(I^2,N*(N+3)/2));
Vmax3acD=sym(zeros(I^2,N*(N+3)/2));
Vmax4acA=sym(zeros(I^2,N*(N+3)/2));
Vmax4acB=sym(zeros(I^2,N*(N+3)/2));
Vmax4acC=sym(zeros(I^2,N*(N+3)/2));
Vmax4acD=sym(zeros(I^2,N*(N+3)/2));
Vmax5acA=sym(zeros(I^2,N*(N+3)/2));
Vmax5acB=sym(zeros(I^2,N*(N+3)/2));
Vmax5acC=sym(zeros(I^2,N*(N+3)/2));
Vmax5acD=sym(zeros(I^2,N*(N+3)/2));
Vmax6acA=sym(zeros(I^2,N*(N+3)/2));
Vmax6acB=sym(zeros(I^2,N*(N+3)/2));
Vmax6acC=sym(zeros(I^2,N*(N+3)/2));
Vmax6acD=sym(zeros(I^2,N*(N+3)/2));
Vmax7acA=sym(zeros(I^2,N*(N+3)/2));
Vmax7acB=sym(zeros(I^2,N*(N+3)/2));
Vmax7acC=sym(zeros(I^2,N*(N+3)/2));
Vmax7acD=sym(zeros(I^2,N*(N+3)/2));
func12_1A_cell=cell(I^2,N*(N+3)/2);
func12_1B_cell=cell(I^2,N*(N+3)/2);
func12_1C_cell=cell(I^2,N*(N+3)/2);
func12_1D_cell=cell(I^2,N*(N+3)/2);
func12_2A_cell=cell(I^2,N*(N+3)/2);
func12_2B_cell=cell(I^2,N*(N+3)/2);
func12_2C_cell=cell(I^2,N*(N+3)/2);
func12_2D_cell=cell(I^2,N*(N+3)/2);
func12_3A_cell=cell(I^2,N*(N+3)/2);
func12_3B_cell=cell(I^2,N*(N+3)/2);
func12_3C_cell=cell(I^2,N*(N+3)/2);
func12_3D_cell=cell(I^2,N*(N+3)/2);
func12_4A_cell=cell(I^2,N*(N+3)/2);
func12_4B_cell=cell(I^2,N*(N+3)/2);
func12_4C_cell=cell(I^2,N*(N+3)/2);
func12_4D_cell=cell(I^2,N*(N+3)/2);
func12_5A_cell=cell(I^2,N*(N+3)/2);
func12_5B_cell=cell(I^2,N*(N+3)/2);
func12_5C_cell=cell(I^2,N*(N+3)/2);
func12_5D_cell=cell(I^2,N*(N+3)/2);
func12_6A_cell=cell(I^2,N*(N+3)/2);
func12_6B_cell=cell(I^2,N*(N+3)/2);
func12_6C_cell=cell(I^2,N*(N+3)/2);
func12_6D_cell=cell(I^2,N*(N+3)/2);
func12_7A_cell=cell(I^2,N*(N+3)/2);
func12_7B_cell=cell(I^2,N*(N+3)/2);
func12_7C_cell=cell(I^2,N*(N+3)/2);
func12_7D_cell=cell(I^2,N*(N+3)/2);
func13_1A_cell=cell(I^2,N*(N+3)/2);
func13_1B_cell=cell(I^2,N*(N+3)/2);
func13_1C_cell=cell(I^2,N*(N+3)/2);
func13_1D_cell=cell(I^2,N*(N+3)/2);
func13_2A_cell=cell(I^2,N*(N+3)/2);
func13_2B_cell=cell(I^2,N*(N+3)/2);
func13_2C_cell=cell(I^2,N*(N+3)/2);
func13_2D_cell=cell(I^2,N*(N+3)/2);
func13_3A_cell=cell(I^2,N*(N+3)/2);
func13_3B_cell=cell(I^2,N*(N+3)/2);
func13_3C_cell=cell(I^2,N*(N+3)/2);
func13_3D_cell=cell(I^2,N*(N+3)/2);
func13_4A_cell=cell(I^2,N*(N+3)/2);
func13_4B_cell=cell(I^2,N*(N+3)/2);
func13_4C_cell=cell(I^2,N*(N+3)/2);
func13_4D_cell=cell(I^2,N*(N+3)/2);
func13_5A_cell=cell(I^2,N*(N+3)/2);
func13_5B_cell=cell(I^2,N*(N+3)/2);
func13_5C_cell=cell(I^2,N*(N+3)/2);
func13_5D_cell=cell(I^2,N*(N+3)/2);
func13_6A_cell=cell(I^2,N*(N+3)/2);
func13_6B_cell=cell(I^2,N*(N+3)/2);
func13_6C_cell=cell(I^2,N*(N+3)/2);
func13_6D_cell=cell(I^2,N*(N+3)/2);
func13_7A_cell=cell(I^2,N*(N+3)/2);
func13_7B_cell=cell(I^2,N*(N+3)/2);
func13_7C_cell=cell(I^2,N*(N+3)/2);
func13_7D_cell=cell(I^2,N*(N+3)/2);
parfor k=1:N*(N+3)/2
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum strain energy Vmax (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to a(ij) and b(k)
	Vmax1abA(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,x,2).*WSA(k)+2.*diff(WBC,x,1).*diff(WSA(k),x,1)+WBC.*diff(WSA(k),x,2));

	Vmax1abB(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,x,2).*WSB(k)+2.*diff(WBC,x,1).*diff(WSB(k),x,1)+WBC.*diff(WSB(k),x,2));

	Vmax1abC(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,x,2).*WSC(k)+2.*diff(WBC,x,1).*diff(WSC(k),x,1)+WBC.*diff(WSC(k),x,2));
	
	Vmax1abD(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,x,2).*WSD(k)+2.*diff(WBC,x,1).*diff(WSD(k),x,1)+WBC.*diff(WSD(k),x,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax2abA(:,k)=(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1)).*(diff(WBC,y,2).*WSA(k)+2.*diff(WBC,y,1).*diff(WSA(k),y,1)+WBC.*diff(WSA(k),y,2));
	
	Vmax2abB(:,k)=(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1)).*(diff(WBC,y,2).*WSB(k)+2.*diff(WBC,y,1).*diff(WSB(k),y,1)+WBC.*diff(WSB(k),y,2));
	
	Vmax2abC(:,k)=(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1)).*(diff(WBC,y,2).*WSC(k)+2.*diff(WBC,y,1).*diff(WSC(k),y,1)+WBC.*diff(WSC(k),y,2));
	
	Vmax2abD(:,k)=(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1)).*(diff(WBC,y,2).*WSD(k)+2.*diff(WBC,y,1).*diff(WSD(k),y,1)+WBC.*diff(WSD(k),y,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax3abA(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*WSA(k)+2.*diff(WBC,y,1).*diff(WSA(k),y,1)+WBC.*diff(WSA(k),y,2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*WSA(k)+2.*diff(WBC,x,1).*diff(WSA(k),x,1)+WBC.*diff(WSA(k),x,2));
	
	Vmax3abB(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*WSB(k)+2.*diff(WBC,y,1).*diff(WSB(k),y,1)+WBC.*diff(WSB(k),y,2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*WSB(k)+2.*diff(WBC,x,1).*diff(WSB(k),x,1)+WBC.*diff(WSB(k),x,2));
	
	Vmax3abC(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*WSC(k)+2.*diff(WBC,y,1).*diff(WSC(k),y,1)+WBC.*diff(WSC(k),y,2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*WSC(k)+2.*diff(WBC,x,1).*diff(WSC(k),x,1)+WBC.*diff(WSC(k),x,2));
	
	Vmax3abD(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*WSD(k)+2.*diff(WBC,y,1).*diff(WSD(k),y,1)+WBC.*diff(WSD(k),y,2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*WSD(k)+2.*diff(WBC,x,1).*diff(WSD(k),x,1)+WBC.*diff(WSD(k),x,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax4abA(:,k)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*WSA(k)+diff(WBC,x,1).*diff(WSA(k),y,1)+diff(WBC,y,1).*diff(WSA(k),x,1)+WBC.*diff(WSA(k),x,y));
	
	Vmax4abB(:,k)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*WSB(k)+diff(WBC,x,1).*diff(WSB(k),y,1)+diff(WBC,y,1).*diff(WSB(k),x,1)+WBC.*diff(WSB(k),x,y));
	
	Vmax4abC(:,k)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*WSC(k)+diff(WBC,x,1).*diff(WSC(k),y,1)+diff(WBC,y,1).*diff(WSC(k),x,1)+WBC.*diff(WSC(k),x,y));
	
	Vmax4abD(:,k)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*WSD(k)+diff(WBC,x,1).*diff(WSD(k),y,1)+diff(WBC,y,1).*diff(WSD(k),x,1)+WBC.*diff(WSD(k),x,y));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax5abA(:,k)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*WSA(k)+WBC.*diff(WSA(k),y,1));
	
	Vmax5abB(:,k)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*WSB(k)+WBC.*diff(WSB(k),y,1));
	
	Vmax5abC(:,k)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*WSC(k)+WBC.*diff(WSC(k),y,1));
	
	Vmax5abD(:,k)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*WSD(k)+WBC.*diff(WSD(k),y,1));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax6abA(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,x,1).*WSA(k)+WBC.*diff(WSA(k),x,1));
	
	Vmax6abB(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,x,1).*WSB(k)+WBC.*diff(WSB(k),x,1));
	
	Vmax6abC(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,x,1).*WSC(k)+WBC.*diff(WSC(k),x,1));
	
	Vmax6abD(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,x,1).*WSD(k)+WBC.*diff(WSD(k),x,1));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax7abA(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*WSA(k)+WBC.*diff(WSA(k),y,1))+(diff(WBC,x,1).*WSA(k)+WBC.*diff(WSA(k),x,1)).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
	
	Vmax7abB(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*WSB(k)+WBC.*diff(WSB(k),y,1))+(diff(WBC,x,1).*WSB(k)+WBC.*diff(WSB(k),x,1)).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
	
	Vmax7abC(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*WSC(k)+WBC.*diff(WSC(k),y,1))+(diff(WBC,x,1).*WSC(k)+WBC.*diff(WSC(k),x,1)).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
	
	Vmax7abD(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*WSD(k)+WBC.*diff(WSD(k),y,1))+(diff(WBC,x,1).*WSD(k)+WBC.*diff(WSD(k),x,1)).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
	
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum strain energy Vmax (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to a(ij) and c(k)
	Vmax1acA(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,x,2).*WAA(k)+2.*diff(WBC,x,1).*diff(WAA(k),x,1)+WBC.*diff(WAA(k),x,2));
	
	Vmax1acB(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,x,2).*WAB(k)+2.*diff(WBC,x,1).*diff(WAB(k),x,1)+WBC.*diff(WAB(k),x,2));
	
	Vmax1acC(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,x,2).*WAC(k)+2.*diff(WBC,x,1).*diff(WAC(k),x,1)+WBC.*diff(WAC(k),x,2));
	
	Vmax1acD(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,x,2).*WAD(k)+2.*diff(WBC,x,1).*diff(WAD(k),x,1)+WBC.*diff(WAD(k),x,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax2acA(:,k)=(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1)).*(diff(WBC,y,2).*WAA(k)+2.*diff(WBC,y,1).*diff(WAA(k),y,1)+WBC.*diff(WAA(k),y,2));
	
	Vmax2acB(:,k)=(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1)).*(diff(WBC,y,2).*WAB(k)+2.*diff(WBC,y,1).*diff(WAB(k),y,1)+WBC.*diff(WAB(k),y,2));
	
	Vmax2acC(:,k)=(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1)).*(diff(WBC,y,2).*WAC(k)+2.*diff(WBC,y,1).*diff(WAC(k),y,1)+WBC.*diff(WAC(k),y,2));
	
	Vmax2acD(:,k)=(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*diff(PhiY(col1),1).*PhiX(row1)+WBC.*diff(PhiY(col1),2).*PhiX(row1)).*(diff(WBC,y,2).*WAD(k)+2.*diff(WBC,y,1).*diff(WAD(k),y,1)+WBC.*diff(WAD(k),y,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax3acA(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*WAA(k)+2.*diff(WBC,y,1).*diff(WAA(k),y,1)+WBC.*diff(WAA(k),y,2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*WAA(k)+2.*diff(WBC,x,1).*diff(WAA(k),x,1)+WBC.*diff(WAA(k),x,2));
	
	Vmax3acB(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*WAB(k)+2.*diff(WBC,y,1).*diff(WAB(k),y,1)+WBC.*diff(WAB(k),y,2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*WAB(k)+2.*diff(WBC,x,1).*diff(WAB(k),x,1)+WBC.*diff(WAB(k),x,2));
	
	Vmax3acC(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*WAC(k)+2.*diff(WBC,y,1).*diff(WAC(k),y,1)+WBC.*diff(WAC(k),y,2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*WAC(k)+2.*diff(WBC,x,1).*diff(WAC(k),x,1)+WBC.*diff(WAC(k),x,2));
	
	Vmax3acD(:,k)=(diff(WBC,x,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,x,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),2).*PhiY(col1)).*(diff(WBC,y,2).*WAD(k)+2.*diff(WBC,y,1).*diff(WAD(k),y,1)+WBC.*diff(WAD(k),y,2))+(diff(WBC,y,2).*PhiX(row1).*PhiY(col1)+2.*diff(WBC,y,1).*PhiX(row1).*diff(PhiY(col1),1)+WBC.*PhiX(row1).*diff(PhiY(col1),2)).*(diff(WBC,x,2).*WAD(k)+2.*diff(WBC,x,1).*diff(WAD(k),x,1)+WBC.*diff(WAD(k),x,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax4acA(:,k)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*WAA(k)+diff(WBC,x,1).*diff(WAA(k),y,1)+diff(WBC,y,1).*diff(WAA(k),x,1)+WBC.*diff(WAA(k),x,y));
	
	Vmax4acB(:,k)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*WAB(k)+diff(WBC,x,1).*diff(WAB(k),y,1)+diff(WBC,y,1).*diff(WAB(k),x,1)+WBC.*diff(WAB(k),x,y));
	
	Vmax4acC(:,k)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*WAC(k)+diff(WBC,x,1).*diff(WAC(k),y,1)+diff(WBC,y,1).*diff(WAC(k),x,1)+WBC.*diff(WAC(k),x,y));
	
	Vmax4acD(:,k)=(diff(WBC,x,y).*PhiX(row1).*PhiY(col1)+diff(WBC,x,1).*PhiX(row1).*diff(PhiY(col1),1)+diff(WBC,y,1).*diff(PhiX(row1),1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*diff(PhiY(col1),1)).*(diff(WBC,x,y).*WAD(k)+diff(WBC,x,1).*diff(WAD(k),y,1)+diff(WBC,y,1).*diff(WAD(k),x,1)+WBC.*diff(WAD(k),x,y));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax5acA(:,k)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*WAA(k)+WBC.*diff(WAA(k),y,1));
	
	Vmax5acB(:,k)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*WAB(k)+WBC.*diff(WAB(k),y,1));
	
	Vmax5acC(:,k)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*WAC(k)+WBC.*diff(WAC(k),y,1));
	
	Vmax5acD(:,k)=(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1)).*(diff(WBC,y,1).*WAD(k)+WBC.*diff(WAD(k),y,1));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax6acA(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,x,1).*WAA(k)+WBC.*diff(WAA(k),x,1));
	
	Vmax6acB(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,x,1).*WAB(k)+WBC.*diff(WAB(k),x,1));
	
	Vmax6acC(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,x,1).*WAC(k)+WBC.*diff(WAC(k),x,1));
	
	Vmax6acD(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,x,1).*WAD(k)+WBC.*diff(WAD(k),x,1));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax7acA(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*WAA(k)+WBC.*diff(WAA(k),y,1))+(diff(WBC,x,1).*WAA(k)+WBC.*diff(WAA(k),x,1)).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
	
	Vmax7acB(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*WAB(k)+WBC.*diff(WAB(k),y,1))+(diff(WBC,x,1).*WAB(k)+WBC.*diff(WAB(k),x,1)).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
	
	Vmax7acC(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*WAC(k)+WBC.*diff(WAC(k),y,1))+(diff(WBC,x,1).*WAC(k)+WBC.*diff(WAC(k),x,1)).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
	
	Vmax7acD(:,k)=(diff(WBC,x,1).*PhiX(row1).*PhiY(col1)+WBC.*diff(PhiX(row1),1).*PhiY(col1)).*(diff(WBC,y,1).*WAD(k)+WBC.*diff(WAD(k),y,1))+(diff(WBC,x,1).*WAD(k)+WBC.*diff(WAD(k),x,1)).*(diff(WBC,y,1).*PhiX(row1).*PhiY(col1)+WBC.*PhiX(row1).*diff(PhiY(col1),1));
end	
for k=1:N*(N+3)/2
	for u=1:I^2
		%%%%%%%%%%%%%%%%%%% integration %%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%
		if size(symvar(Vmax1abA(u,k)),2)==2
			func12_1A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax1abA(u,k))]);
			K12_1A(u,k)=D*integral2(func12_1A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_1A(u,k)=D*int(int(Vmax1abA(u,k),x,0,x0),y,y_crackline,b);
		end
		if size(symvar(Vmax1abB(u,k)),2)==2
			func12_1B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax1abB(u,k))]);
			K12_1B(u,k)=D*integral2(func12_1B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_1B(u,k)=D*int(int(Vmax1abB(u,k),x,0,x0),y,0,y_crackline);
		end
		if size(symvar(Vmax1abC(u,k)),2)==2
			func12_1C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax1abC(u,k))]);
			K12_1C(u,k)=D*integral2(func12_1C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_1C(u,k)=D*int(int(Vmax1abC(u,k),x,x0,a),y,0,y_crackline);
		end
		if size(symvar(Vmax1abD(u,k)),2)==2
			func12_1D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax1abD(u,k))]);
			K12_1D(u,k)=D*integral2(func12_1D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_1D(u,k)=D*int(int(Vmax1abD(u,k),x,x0,a),y,y_crackline,b);
		end
		%%%%%%%%%%%
		if size(symvar(Vmax2abA(u,k)),2)==2
			func12_2A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax2abA(u,k))]);
			K12_2A(u,k)=D*integral2(func12_2A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_2A(u,k)=D*int(int(Vmax2abA(u,k),x,0,x0),y,y_crackline,b);
		end
		if size(symvar(Vmax2abB(u,k)),2)==2
			func12_2B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax2abB(u,k))]);
			K12_2B(u,k)=D*integral2(func12_2B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_2B(u,k)=D*int(int(Vmax2abB(u,k),x,0,x0),y,0,y_crackline);
		end
		if size(symvar(Vmax2abC(u,k)),2)==2
			func12_2C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax2abC(u,k))]);
			K12_2C(u,k)=D*integral2(func12_2C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_2C(u,k)=D*int(int(Vmax2abC(u,k),x,x0,a),y,0,y_crackline);
		end
		if size(symvar(Vmax2abD(u,k)),2)==2
			func12_2D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax2abD(u,k))]);
			K12_2D(u,k)=D*integral2(func12_2D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_2D(u,k)=D*int(int(Vmax2abD(u,k),x,x0,a),y,y_crackline,b);
		end
		%%%%%%%%%%%
		if size(symvar(Vmax3abA(u,k)),2)==2
			func12_3A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax3abA(u,k))]);
			K12_3A(u,k)=D*nu*integral2(func12_3A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_3A(u,k)=D*nu*int(int(Vmax3abA(u,k),x,0,x0),y,y_crackline,b);
		end
		if size(symvar(Vmax3abB(u,k)),2)==2
			func12_3B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax3abB(u,k))]);
			K12_3B(u,k)=D*nu*integral2(func12_3B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_3B(u,k)=D*nu*int(int(Vmax3abB(u,k),x,0,x0),y,0,y_crackline);
		end
		if size(symvar(Vmax3abC(u,k)),2)==2
			func12_3C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax3abC(u,k))]);
			K12_3C(u,k)=D*nu*integral2(func12_3C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_3C(u,k)=D*nu*int(int(Vmax3abC(u,k),x,x0,a),y,0,y_crackline);
		end
		if size(symvar(Vmax3abD(u,k)),2)==2
			func12_3D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax3abD(u,k))]);
			K12_3D(u,k)=D*nu*integral2(func12_3D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_3D(u,k)=D*nu*int(int(Vmax3abD(u,k),x,x0,a),y,y_crackline,b);
		end
		%%%%%%%%%%%
		if size(symvar(Vmax4abA(u,k)),2)==2
			func12_4A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax4abA(u,k))]);
			K12_4A(u,k)=D*2*(1-nu)*integral2(func12_4A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_4A(u,k)=D*2*(1-nu)*int(int(Vmax4abA(u,k),x,0,x0),y,y_crackline,b);
		end
		if size(symvar(Vmax4abB(u,k)),2)==2
			func12_4B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax4abB(u,k))]);
			K12_4B(u,k)=D*2*(1-nu)*integral2(func12_4B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_4B(u,k)=D*2*(1-nu)*int(int(Vmax4abB(u,k),x,0,x0),y,0,y_crackline);
		end
		if size(symvar(Vmax4abC(u,k)),2)==2
			func12_4C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax4abC(u,k))]);
			K12_4C(u,k)=D*2*(1-nu)*integral2(func12_4C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_4C(u,k)=D*2*(1-nu)*int(int(Vmax4abC(u,k),x,x0,a),y,0,y_crackline);
		end
		if size(symvar(Vmax4abD(u,k)),2)==2
			func12_4D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax4abD(u,k))]);
			K12_4D(u,k)=D*2*(1-nu)*integral2(func12_4D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K12_4D(u,k)=D*2*(1-nu)*int(int(Vmax4abD(u,k),x,x0,a),y,y_crackline,b);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5abA(u,k)),2)==2
			func12_5A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax5abA(u,k))]);
			for j=1:2							%for area A1 and A2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_5A(u,k)=K12_5A(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5A_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5abA(u,k)),2)==1
			if strcmp(char(symvar(Vmax5abA(u,k))),'x')==1			%if the symvar is x
				func12_5A_cell{u,k}=eval(['@(x)',vectorize(Vmax5abA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_5A(u,k)=K12_5A(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5A_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_5A_cell{u,k}=eval(['@(y)',vectorize(Vmax5abA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_5A(u,k)=K12_5A(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5A_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_5A(u,k)=K12_5A(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5abB(u,k)),2)==2
			func12_5B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax5abB(u,k))]);
			for j=3:4							%for area B1 and B2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_5B(u,k)=K12_5B(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5B_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5abB(u,k)),2)==1
			if strcmp(char(symvar(Vmax5abB(u,k))),'x')==1			%if the symvar is x
				func12_5B_cell{u,k}=eval(['@(x)',vectorize(Vmax5abB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_5B(u,k)=K12_5B(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5B_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_5B_cell{u,k}=eval(['@(y)',vectorize(Vmax5abB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_5B(u,k)=K12_5B(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5B_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_5B(u,k)=K12_5B(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5abC(u,k)),2)==2
			func12_5C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax5abC(u,k))]);
			for j=5:6							%for area C1 and C2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_5C(u,k)=K12_5C(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5C_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5abC(u,k)),2)==1
			if strcmp(char(symvar(Vmax5abC(u,k))),'x')==1			%if the symvar is x
				func12_5C_cell{u,k}=eval(['@(x)',vectorize(Vmax5abC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_5C(u,k)=K12_5C(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5C_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_5C_cell{u,k}=eval(['@(y)',vectorize(Vmax5abC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_5C(u,k)=K12_5C(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5C_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_5C(u,k)=K12_5C(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5abD(u,k)),2)==2
			func12_5D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax5abD(u,k))]);
			for j=7:8							%for area D1 and D2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_5D(u,k)=K12_5D(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5D_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5abD(u,k)),2)==1
			if strcmp(char(symvar(Vmax5abD(u,k))),'x')==1			%if the symvar is x
				func12_5D_cell{u,k}=eval(['@(x)',vectorize(Vmax5abD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_5D(u,k)=K12_5D(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5D_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_5D_cell{u,k}=eval(['@(y)',vectorize(Vmax5abD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_5D(u,k)=K12_5D(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func12_5D_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_5D(u,k)=K12_5D(u,k);
		end
		%%%%%%%%%%%
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax6abA(u,k)),2)==2
			func12_6A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax6abA(u,k))]);
			for j=1:2							%for area A1 and A2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_6A(u,k)=K12_6A(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6A_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6abA(u,k)),2)==1
			if strcmp(char(symvar(Vmax6abA(u,k))),'x')==1			%if the symvar is x
				func12_6A_cell{u,k}=eval(['@(x)',vectorize(Vmax6abA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_6A(u,k)=K12_6A(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6A_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_6A_cell{u,k}=eval(['@(y)',vectorize(Vmax6abA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_6A(u,k)=K12_6A(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6A_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_6A(u,k)=K12_6A(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax6abB(u,k)),2)==2
			func12_6B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax6abB(u,k))]);
			for j=3:4							%for area B1 and B2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_6B(u,k)=K12_6B(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6B_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6abB(u,k)),2)==1
			if strcmp(char(symvar(Vmax6abB(u,k))),'x')==1			%if the symvar is x
				func12_6B_cell{u,k}=eval(['@(x)',vectorize(Vmax6abB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_6B(u,k)=K12_6B(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6B_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_6B_cell{u,k}=eval(['@(y)',vectorize(Vmax6abB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_6B(u,k)=K12_6B(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6B_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_6B(u,k)=K12_6B(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax6abC(u,k)),2)==2
			func12_6C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax6abC(u,k))]);
			for j=5:6							%for area C1 and C2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_6C(u,k)=K12_6C(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6C_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6abC(u,k)),2)==1
			if strcmp(char(symvar(Vmax6abC(u,k))),'x')==1			%if the symvar is x
				func12_6C_cell{u,k}=eval(['@(x)',vectorize(Vmax6abC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_6C(u,k)=K12_6C(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6C_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_6C_cell{u,k}=eval(['@(y)',vectorize(Vmax6abC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_6C(u,k)=K12_6C(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6C_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_6C(u,k)=K12_6C(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax6abD(u,k)),2)==2
			func12_6D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax6abD(u,k))]);
			for j=7:8							%for area D1 and D2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_6D(u,k)=K12_6D(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6D_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6abD(u,k)),2)==1
			if strcmp(char(symvar(Vmax6abD(u,k))),'x')==1			%if the symvar is x
				func12_6D_cell{u,k}=eval(['@(x)',vectorize(Vmax6abD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_6D(u,k)=K12_6D(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6D_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_6D_cell{u,k}=eval(['@(y)',vectorize(Vmax6abD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_6D(u,k)=K12_6D(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func12_6D_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_6D(u,k)=K12_6D(u,k);
		end
		% %%%%%%%%%%%
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax7abA(u,k)),2)==2
			func12_7A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax7abA(u,k))]);
			for j=1:2							%for area A1 and A2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_7A(u,k)=K12_7A(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7A_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7abA(u,k)),2)==1
			if strcmp(char(symvar(Vmax7abA(u,k))),'x')==1			%if the symvar is x
				func12_7A_cell{u,k}=eval(['@(x)',vectorize(Vmax7abA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_7A(u,k)=K12_7A(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7A_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_7A_cell{u,k}=eval(['@(y)',vectorize(Vmax7abA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_7A(u,k)=K12_7A(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7A_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_7A(u,k)=K12_7A(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax7abB(u,k)),2)==2
			func12_7B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax7abB(u,k))]);
			for j=3:4							%for area B1 and B2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_7B(u,k)=K12_7B(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7B_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7abB(u,k)),2)==1
			if strcmp(char(symvar(Vmax7abB(u,k))),'x')==1			%if the symvar is x
				func12_7B_cell{u,k}=eval(['@(x)',vectorize(Vmax7abB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_7B(u,k)=K12_7B(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7B_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_7B_cell{u,k}=eval(['@(y)',vectorize(Vmax7abB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_7B(u,k)=K12_7B(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7B_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_7B(u,k)=K12_7B(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax7abC(u,k)),2)==2
			func12_7C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax7abC(u,k))]);
			for j=5:6							%for area C1 and C2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_7C(u,k)=K12_7C(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7C_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7abC(u,k)),2)==1
			if strcmp(char(symvar(Vmax7abC(u,k))),'x')==1			%if the symvar is x
				func12_7C_cell{u,k}=eval(['@(x)',vectorize(Vmax7abC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_7C(u,k)=K12_7C(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7C_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_7C_cell{u,k}=eval(['@(y)',vectorize(Vmax7abC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_7C(u,k)=K12_7C(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7C_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_7C(u,k)=K12_7C(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax7abD(u,k)),2)==2
			func12_7D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax7abD(u,k))]);
			for j=7:8							%for area D1 and D2
				for i=1:num_triangles
					for q=1:num_intpoints
						K12_7D(u,k)=K12_7D(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7D_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7abD(u,k)),2)==1
			if strcmp(char(symvar(Vmax7abD(u,k))),'x')==1			%if the symvar is x
				func12_7D_cell{u,k}=eval(['@(x)',vectorize(Vmax7abD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_7D(u,k)=K12_7D(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7D_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func12_7D_cell{u,k}=eval(['@(y)',vectorize(Vmax7abD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K12_7D(u,k)=K12_7D(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func12_7D_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K12_7D(u,k)=K12_7D(u,k);
		end
		% %%%%%%%%%%%
		if size(symvar(Vmax1acA(u,k)),2)==2
			func13_1A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax1acA(u,k))]);
			K13_1A(u,k)=D*integral2(func13_1A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_1A(u,k)=D*int(int(Vmax1acA(u,k),x,0,x0),y,y_crackline,b);
		end
		if size(symvar(Vmax1acB(u,k)),2)==2
			func13_1B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax1acB(u,k))]);
			K13_1B(u,k)=D*integral2(func13_1B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_1B(u,k)=D*int(int(Vmax1acB(u,k),x,0,x0),y,0,y_crackline);
		end
		if size(symvar(Vmax1acC(u,k)),2)==2
			func13_1C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax1acC(u,k))]);
			K13_1C(u,k)=D*integral2(func13_1C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_1C(u,k)=D*int(int(Vmax1acC(u,k),x,x0,a),y,0,y_crackline);
		end
		if size(symvar(Vmax1acD(u,k)),2)==2
			func13_1D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax1acD(u,k))]);
			K13_1D(u,k)=D*integral2(func13_1D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_1D(u,k)=D*int(int(Vmax1acD(u,k),x,x0,a),y,y_crackline,b);
		end
		%%%%%%%%%%%
		if size(symvar(Vmax2acA(u,k)),2)==2
			func13_2A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax2acA(u,k))]);
			K13_2A(u,k)=D*integral2(func13_2A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_2A(u,k)=D*int(int(Vmax2acA(u,k),x,0,x0),y,y_crackline,b);
		end
		if size(symvar(Vmax2acB(u,k)),2)==2
			func13_2B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax2acB(u,k))]);
			K13_2B(u,k)=D*integral2(func13_2B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_2B(u,k)=D*int(int(Vmax2acB(u,k),x,0,x0),y,0,y_crackline);
		end
		if size(symvar(Vmax2acC(u,k)),2)==2
			func13_2C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax2acC(u,k))]);
			K13_2C(u,k)=D*integral2(func13_2C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_2C(u,k)=D*int(int(Vmax2acC(u,k),x,x0,a),y,0,y_crackline);
		end
		if size(symvar(Vmax2acD(u,k)),2)==2
			func13_2D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax2acD(u,k))]);
			K13_2D(u,k)=D*integral2(func13_2D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_2D(u,k)=D*int(int(Vmax2acD(u,k),x,x0,a),y,y_crackline,b);
		end
		%%%%%%%%%%%
		if size(symvar(Vmax3acA(u,k)),2)==2
			func13_3A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax3acA(u,k))]);
			K13_3A(u,k)=D*nu*integral2(func13_3A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_3A(u,k)=D*nu*int(int(Vmax3acA(u,k),x,0,x0),y,y_crackline,b);
		end
		if size(symvar(Vmax3acB(u,k)),2)==2
			func13_3B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax3acB(u,k))]);
			K13_3B(u,k)=D*nu*integral2(func13_3B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_3B(u,k)=D*nu*int(int(Vmax3acB(u,k),x,0,x0),y,0,y_crackline);
		end
		if size(symvar(Vmax3acC(u,k)),2)==2
			func13_3C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax3acC(u,k))]);
			K13_3C(u,k)=D*nu*integral2(func13_3C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_3C(u,k)=D*nu*int(int(Vmax3acC(u,k),x,x0,a),y,0,y_crackline);
		end
		if size(symvar(Vmax3acD(u,k)),2)==2
			func13_3D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax3acD(u,k))]);
			K13_3D(u,k)=D*nu*integral2(func13_3D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_3D(u,k)=D*nu*int(int(Vmax3acD(u,k),x,x0,a),y,y_crackline,b);
		end
		%%%%%%%%%%%
		if size(symvar(Vmax4acA(u,k)),2)==2
			func13_4A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax4acA(u,k))]);
			K13_4A(u,k)=D*2*(1-nu)*integral2(func13_4A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_4A(u,k)=D*2*(1-nu)*int(int(Vmax4acA(u,k),x,0,x0),y,y_crackline,b);
		end
		if size(symvar(Vmax4acB(u,k)),2)==2
			func13_4B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax4acB(u,k))]);
			K13_4B(u,k)=D*2*(1-nu)*integral2(func13_4B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_4B(u,k)=D*2*(1-nu)*int(int(Vmax4acB(u,k),x,0,x0),y,0,y_crackline);
		end
		if size(symvar(Vmax4acC(u,k)),2)==2
			func13_4C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax4acC(u,k))]);
			K13_4C(u,k)=D*2*(1-nu)*integral2(func13_4C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_4C(u,k)=D*2*(1-nu)*int(int(Vmax4acC(u,k),x,x0,a),y,0,y_crackline);
		end
		if size(symvar(Vmax4acD(u,k)),2)==2
			func13_4D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax4acD(u,k))]);
			K13_4D(u,k)=D*2*(1-nu)*integral2(func13_4D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			K13_4D(u,k)=D*2*(1-nu)*int(int(Vmax4acD(u,k),x,x0,a),y,y_crackline,b);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5acA(u,k)),2)==2
			func13_5A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax5acA(u,k))]);
			for j=1:2							%for area A1 and A2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_5A(u,k)=K13_5A(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5A_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5acA(u,k)),2)==1
			if strcmp(char(symvar(Vmax5acA(u,k))),'x')==1			%if the symvar is x
				func13_5A_cell{u,k}=eval(['@(x)',vectorize(Vmax5acA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_5A(u,k)=K13_5A(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5A_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_5A_cell{u,k}=eval(['@(y)',vectorize(Vmax5acA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_5A(u,k)=K13_5A(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5A_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_5A(u,k)=K13_5A(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5acB(u,k)),2)==2
			func13_5B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax5acB(u,k))]);
			for j=3:4							%for area B1 and B2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_5B(u,k)=K13_5B(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5B_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5acB(u,k)),2)==1
			if strcmp(char(symvar(Vmax5acB(u,k))),'x')==1			%if the symvar is x
				func13_5B_cell{u,k}=eval(['@(x)',vectorize(Vmax5acB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_5B(u,k)=K13_5B(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5B_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_5B_cell{u,k}=eval(['@(y)',vectorize(Vmax5acB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_5B(u,k)=K13_5B(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5B_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_5B(u,k)=K13_5B(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5acC(u,k)),2)==2
			func13_5C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax5acC(u,k))]);
			for j=5:6							%for area C1 and C2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_5C(u,k)=K13_5C(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5C_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5acC(u,k)),2)==1
			if strcmp(char(symvar(Vmax5acC(u,k))),'x')==1			%if the symvar is x
				func13_5C_cell{u,k}=eval(['@(x)',vectorize(Vmax5acC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_5C(u,k)=K13_5C(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5C_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_5C_cell{u,k}=eval(['@(y)',vectorize(Vmax5acC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_5C(u,k)=K13_5C(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5C_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_5C(u,k)=K13_5C(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax5acD(u,k)),2)==2
			func13_5D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax5acD(u,k))]);
			for j=7:8							%for area D1 and D2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_5D(u,k)=K13_5D(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5D_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax5acD(u,k)),2)==1
			if strcmp(char(symvar(Vmax5acD(u,k))),'x')==1			%if the symvar is x
				func13_5D_cell{u,k}=eval(['@(x)',vectorize(Vmax5acD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_5D(u,k)=K13_5D(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5D_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_5D_cell{u,k}=eval(['@(y)',vectorize(Vmax5acD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_5D(u,k)=K13_5D(u,k)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func13_5D_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_5D(u,k)=K13_5D(u,k);
		end
		%%%%%%%%%%%
		if size(symvar(Vmax6acA(u,k)),2)==2
			func13_6A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax6acA(u,k))]);
			for j=1:2							%for area A1 and A2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_6A(u,k)=K13_6A(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6A_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6acA(u,k)),2)==1
			if strcmp(char(symvar(Vmax6acA(u,k))),'x')==1			%if the symvar is x
				func13_6A_cell{u,k}=eval(['@(x)',vectorize(Vmax6acA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_6A(u,k)=K13_6A(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6A_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_6A_cell{u,k}=eval(['@(y)',vectorize(Vmax6acA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_6A(u,k)=K13_6A(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6A_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_6A(u,k)=K13_6A(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax6acB(u,k)),2)==2
			func13_6B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax6acB(u,k))]);
			for j=3:4							%for area B1 and B2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_6B(u,k)=K13_6B(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6B_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6acB(u,k)),2)==1
			if strcmp(char(symvar(Vmax6acB(u,k))),'x')==1			%if the symvar is x
				func13_6B_cell{u,k}=eval(['@(x)',vectorize(Vmax6acB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_6B(u,k)=K13_6B(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6B_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_6B_cell{u,k}=eval(['@(y)',vectorize(Vmax6acB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_6B(u,k)=K13_6B(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6B_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_6B(u,k)=K13_6B(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax6acC(u,k)),2)==2
			func13_6C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax6acC(u,k))]);
			for j=5:6							%for area C1 and C2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_6C(u,k)=K13_6C(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6C_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6acC(u,k)),2)==1
			if strcmp(char(symvar(Vmax6acC(u,k))),'x')==1			%if the symvar is x
				func13_6C_cell{u,k}=eval(['@(x)',vectorize(Vmax6acC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_6C(u,k)=K13_6C(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6C_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_6C_cell{u,k}=eval(['@(y)',vectorize(Vmax6acC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_6C(u,k)=K13_6C(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6C_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_6C(u,k)=K13_6C(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax6acD(u,k)),2)==2
			func13_6D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax6acD(u,k))]);
			for j=7:8							%for area D1 and D2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_6D(u,k)=K13_6D(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6D_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax6acD(u,k)),2)==1
			if strcmp(char(symvar(Vmax6acD(u,k))),'x')==1			%if the symvar is x
				func13_6D_cell{u,k}=eval(['@(x)',vectorize(Vmax6acD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_6D(u,k)=K13_6D(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6D_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_6D_cell{u,k}=eval(['@(y)',vectorize(Vmax6acD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_6D(u,k)=K13_6D(u,k)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func13_6D_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_6D(u,k)=K13_6D(u,k);
		end
		% %%%%%%%%%%%
		if size(symvar(Vmax7acA(u,k)),2)==2
			func13_7A_cell{u,k}=eval(['@(x,y)',vectorize(Vmax7acA(u,k))]);
			for j=1:2							%for area A1 and A2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_7A(u,k)=K13_7A(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7A_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7acA(u,k)),2)==1
			if strcmp(char(symvar(Vmax7acA(u,k))),'x')==1			%if the symvar is x
				func13_7A_cell{u,k}=eval(['@(x)',vectorize(Vmax7acA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_7A(u,k)=K13_7A(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7A_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_7A_cell{u,k}=eval(['@(y)',vectorize(Vmax7acA(u,k))]);
				for j=1:2							%for area A1 and A2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_7A(u,k)=K13_7A(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7A_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_7A(u,k)=K13_7A(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax7acB(u,k)),2)==2
			func13_7B_cell{u,k}=eval(['@(x,y)',vectorize(Vmax7acB(u,k))]);
			for j=3:4							%for area B1 and B2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_7B(u,k)=K13_7B(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7B_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7acB(u,k)),2)==1
			if strcmp(char(symvar(Vmax7acB(u,k))),'x')==1			%if the symvar is x
				func13_7B_cell{u,k}=eval(['@(x)',vectorize(Vmax7acB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_7B(u,k)=K13_7B(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7B_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_7B_cell{u,k}=eval(['@(y)',vectorize(Vmax7acB(u,k))]);
				for j=3:4							%for area B1 and B2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_7B(u,k)=K13_7B(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7B_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_7B(u,k)=K13_7B(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax7acC(u,k)),2)==2
			func13_7C_cell{u,k}=eval(['@(x,y)',vectorize(Vmax7acC(u,k))]);
			for j=5:6							%for area C1 and C2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_7C(u,k)=K13_7C(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7C_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7acC(u,k)),2)==1
			if strcmp(char(symvar(Vmax7acC(u,k))),'x')==1			%if the symvar is x
				func13_7C_cell{u,k}=eval(['@(x)',vectorize(Vmax7acC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_7C(u,k)=K13_7C(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7C_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_7C_cell{u,k}=eval(['@(y)',vectorize(Vmax7acC(u,k))]);
				for j=5:6							%for area C1 and C2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_7C(u,k)=K13_7C(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7C_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_7C(u,k)=K13_7C(u,k);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if size(symvar(Vmax7acD(u,k)),2)==2
			func13_7D_cell{u,k}=eval(['@(x,y)',vectorize(Vmax7acD(u,k))]);
			for j=7:8							%for area D1 and D2
				for i=1:num_triangles
					for q=1:num_intpoints
						K13_7D(u,k)=K13_7D(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7D_cell{u,k},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		elseif size(symvar(Vmax7acD(u,k)),2)==1
			if strcmp(char(symvar(Vmax7acD(u,k))),'x')==1			%if the symvar is x
				func13_7D_cell{u,k}=eval(['@(x)',vectorize(Vmax7acD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_7D(u,k)=K13_7D(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7D_cell{u,k},intpoints{j,i,1}(q,1));
						end
					end
				end
			
			else													%if the symvar is y				
				func13_7D_cell{u,k}=eval(['@(y)',vectorize(Vmax7acD(u,k))]);
				for j=7:8							%for area D1 and D2
					for i=1:num_triangles
						for q=1:num_intpoints
							K13_7D(u,k)=K13_7D(u,k)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func13_7D_cell{u,k},intpoints{j,i,1}(q,2));
						end
					end
				end
			end
		else
			K13_7D(u,k)=K13_7D(u,k);
		end
	end
end
K12_1=K12_1A+K12_1B+K12_1C+K12_1D;
K12_2=K12_2A+K12_2B+K12_2C+K12_2D;
K12_3=K12_3A+K12_3B+K12_3C+K12_3D;
K12_4=K12_4A+K12_4B+K12_4C+K12_4D;
K12_5=K12_5A+K12_5B+K12_5C+K12_5D;
K12_6=K12_6A+K12_6B+K12_6C+K12_6D;
K12_7=K12_7A+K12_7B+K12_7C+K12_7D;
K13_1=K13_1A+K13_1B+K13_1C+K13_1D;
K13_2=K13_2A+K13_2B+K13_2C+K13_2D;
K13_3=K13_3A+K13_3B+K13_3C+K13_3D;
K13_4=K13_4A+K13_4B+K13_4C+K13_4D;
K13_5=K13_5A+K13_5B+K13_5C+K13_5D;
K13_6=K13_6A+K13_6B+K13_6C+K13_6D;
K13_7=K13_7A+K13_7B+K13_7C+K13_7D;
K12=K12_1+K12_2+K12_3+K12_4+K12_5+K12_6+K12_7;
K13=K13_1+K13_2+K13_3+K13_4+K13_5+K13_6+K13_7;
toc
tic
K23=zeros(N*(N+3)/2,N*(N+3)/2);
K23_1=zeros(N*(N+3)/2,N*(N+3)/2);
K23_2=zeros(N*(N+3)/2,N*(N+3)/2);
K23_3=zeros(N*(N+3)/2,N*(N+3)/2);
K23_4=zeros(N*(N+3)/2,N*(N+3)/2);
K23_5=zeros(N*(N+3)/2,N*(N+3)/2);
K23_6=zeros(N*(N+3)/2,N*(N+3)/2);
K23_7=zeros(N*(N+3)/2,N*(N+3)/2);
Vmax1bcA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1bcB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1bcC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax1bcD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2bcA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2bcB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2bcC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax2bcD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3bcA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3bcB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3bcC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax3bcD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4bcA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4bcB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4bcC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax4bcD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5bcA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5bcB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5bcC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax5bcD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6bcA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6bcB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6bcC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax6bcD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7bcA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7bcB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7bcC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
Vmax7bcD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
func23_1A_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_1B_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_1C_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_1D_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_2A_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_2B_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_2C_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_2D_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_3A_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_3B_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_3C_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_3D_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_4A_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_4B_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_4C_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_4D_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_5A_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_5B_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_5C_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_5D_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_6A_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_6B_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_6C_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_6D_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_7A_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_7B_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_7C_cell=cell(N*(N+3)/2,N*(N+3)/2);
func23_7D_cell=cell(N*(N+3)/2,N*(N+3)/2);
parfor j=1:N*(N+3)/2
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum strain energy Vmax (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to b(i) and c(j)
	Vmax1bcA(:,j)=(diff(WBC,x,2).*WSA(row_num)+2.*diff(WBC,x,1).*diff(WSA(row_num),x,1)+WBC.*diff(WSA(row_num),x,2)).*(diff(WBC,x,2).*WAA(j)+2.*diff(WBC,x,1).*diff(WAA(j),x,1)+WBC.*diff(WAA(j),x,2));
	
	Vmax1bcB(:,j)=(diff(WBC,x,2).*WSB(row_num)+2.*diff(WBC,x,1).*diff(WSB(row_num),x,1)+WBC.*diff(WSB(row_num),x,2)).*(diff(WBC,x,2).*WAB(j)+2.*diff(WBC,x,1).*diff(WAB(j),x,1)+WBC.*diff(WAB(j),x,2));
	
	Vmax1bcC(:,j)=(diff(WBC,x,2).*WSC(row_num)+2.*diff(WBC,x,1).*diff(WSC(row_num),x,1)+WBC.*diff(WSC(row_num),x,2)).*(diff(WBC,x,2).*WAC(j)+2.*diff(WBC,x,1).*diff(WAC(j),x,1)+WBC.*diff(WAC(j),x,2));
	
	Vmax1bcD(:,j)=(diff(WBC,x,2).*WSD(row_num)+2.*diff(WBC,x,1).*diff(WSD(row_num),x,1)+WBC.*diff(WSD(row_num),x,2)).*(diff(WBC,x,2).*WAD(j)+2.*diff(WBC,x,1).*diff(WAD(j),x,1)+WBC.*diff(WAD(j),x,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax2bcA(:,j)=(diff(WBC,y,2).*WSA(row_num)+2.*diff(WBC,y,1).*diff(WSA(row_num),y,1)+WBC.*diff(WSA(row_num),y,2)).*(diff(WBC,y,2).*WAA(j)+2.*diff(WBC,y,1).*diff(WAA(j),y,1)+WBC.*diff(WAA(j),y,2));
	
	Vmax2bcB(:,j)=(diff(WBC,y,2).*WSB(row_num)+2.*diff(WBC,y,1).*diff(WSB(row_num),y,1)+WBC.*diff(WSB(row_num),y,2)).*(diff(WBC,y,2).*WAB(j)+2.*diff(WBC,y,1).*diff(WAB(j),y,1)+WBC.*diff(WAB(j),y,2));
	
	Vmax2bcC(:,j)=(diff(WBC,y,2).*WSC(row_num)+2.*diff(WBC,y,1).*diff(WSC(row_num),y,1)+WBC.*diff(WSC(row_num),y,2)).*(diff(WBC,y,2).*WAC(j)+2.*diff(WBC,y,1).*diff(WAC(j),y,1)+WBC.*diff(WAC(j),y,2));
	
	Vmax2bcD(:,j)=(diff(WBC,y,2).*WSD(row_num)+2.*diff(WBC,y,1).*diff(WSD(row_num),y,1)+WBC.*diff(WSD(row_num),y,2)).*(diff(WBC,y,2).*WAD(j)+2.*diff(WBC,y,1).*diff(WAD(j),y,1)+WBC.*diff(WAD(j),y,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax3bcA(:,j)=(diff(WBC,x,2).*WSA(row_num)+2.*diff(WBC,x,1).*diff(WSA(row_num),x,1)+WBC.*diff(WSA(row_num),x,2)).*(diff(WBC,y,2).*WAA(j)+2.*diff(WBC,y,1).*diff(WAA(j),y,1)+WBC.*diff(WAA(j),y,2))+(diff(WBC,y,2).*WSA(row_num)+2.*diff(WBC,y,1).*diff(WSA(row_num),y,1)+WBC.*diff(WSA(row_num),y,2)).*(diff(WBC,x,2).*WAA(j)+2.*diff(WBC,x,1).*diff(WAA(j),x,1)+WBC.*diff(WAA(j),x,2));
	
	Vmax3bcB(:,j)=(diff(WBC,x,2).*WSB(row_num)+2.*diff(WBC,x,1).*diff(WSB(row_num),x,1)+WBC.*diff(WSB(row_num),x,2)).*(diff(WBC,y,2).*WAB(j)+2.*diff(WBC,y,1).*diff(WAB(j),y,1)+WBC.*diff(WAB(j),y,2))+(diff(WBC,y,2).*WSB(row_num)+2.*diff(WBC,y,1).*diff(WSB(row_num),y,1)+WBC.*diff(WSB(row_num),y,2)).*(diff(WBC,x,2).*WAB(j)+2.*diff(WBC,x,1).*diff(WAB(j),x,1)+WBC.*diff(WAB(j),x,2));
	
	Vmax3bcC(:,j)=(diff(WBC,x,2).*WSC(row_num)+2.*diff(WBC,x,1).*diff(WSC(row_num),x,1)+WBC.*diff(WSC(row_num),x,2)).*(diff(WBC,y,2).*WAC(j)+2.*diff(WBC,y,1).*diff(WAC(j),y,1)+WBC.*diff(WAC(j),y,2))+(diff(WBC,y,2).*WSC(row_num)+2.*diff(WBC,y,1).*diff(WSC(row_num),y,1)+WBC.*diff(WSC(row_num),y,2)).*(diff(WBC,x,2).*WAC(j)+2.*diff(WBC,x,1).*diff(WAC(j),x,1)+WBC.*diff(WAC(j),x,2));
	
	Vmax3bcD(:,j)=(diff(WBC,x,2).*WSD(row_num)+2.*diff(WBC,x,1).*diff(WSD(row_num),x,1)+WBC.*diff(WSD(row_num),x,2)).*(diff(WBC,y,2).*WAD(j)+2.*diff(WBC,y,1).*diff(WAD(j),y,1)+WBC.*diff(WAD(j),y,2))+(diff(WBC,y,2).*WSD(row_num)+2.*diff(WBC,y,1).*diff(WSD(row_num),y,1)+WBC.*diff(WSD(row_num),y,2)).*(diff(WBC,x,2).*WAD(j)+2.*diff(WBC,x,1).*diff(WAD(j),x,1)+WBC.*diff(WAD(j),x,2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax4bcA(:,j)=(diff(WBC,x,y).*WSA(row_num)+diff(WBC,x,1).*diff(WSA(row_num),y,1)+diff(WBC,y,1).*diff(WSA(row_num),x,1)+WBC.*diff(WSA(row_num),x,y)).*(diff(WBC,x,y).*WAA(j)+diff(WBC,x,1).*diff(WAA(j),y,1)+diff(WBC,y,1).*diff(WAA(j),x,1)+WBC.*diff(WAA(j),x,y));
	
	Vmax4bcB(:,j)=(diff(WBC,x,y).*WSB(row_num)+diff(WBC,x,1).*diff(WSB(row_num),y,1)+diff(WBC,y,1).*diff(WSB(row_num),x,1)+WBC.*diff(WSB(row_num),x,y)).*(diff(WBC,x,y).*WAB(j)+diff(WBC,x,1).*diff(WAB(j),y,1)+diff(WBC,y,1).*diff(WAB(j),x,1)+WBC.*diff(WAB(j),x,y));
	
	Vmax4bcC(:,j)=(diff(WBC,x,y).*WSC(row_num)+diff(WBC,x,1).*diff(WSC(row_num),y,1)+diff(WBC,y,1).*diff(WSC(row_num),x,1)+WBC.*diff(WSC(row_num),x,y)).*(diff(WBC,x,y).*WAC(j)+diff(WBC,x,1).*diff(WAC(j),y,1)+diff(WBC,y,1).*diff(WAC(j),x,1)+WBC.*diff(WAC(j),x,y));
	
	Vmax4bcD(:,j)=(diff(WBC,x,y).*WSD(row_num)+diff(WBC,x,1).*diff(WSD(row_num),y,1)+diff(WBC,y,1).*diff(WSD(row_num),x,1)+WBC.*diff(WSD(row_num),x,y)).*(diff(WBC,x,y).*WAD(j)+diff(WBC,x,1).*diff(WAD(j),y,1)+diff(WBC,y,1).*diff(WAD(j),x,1)+WBC.*diff(WAD(j),x,y));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax5bcA(:,j)=(diff(WBC,y,1).*WSA(row_num)+WBC.*diff(WSA(row_num),y,1)).*(diff(WBC,y,1).*WAA(j)+WBC.*diff(WAA(j),y,1));
	
	Vmax5bcB(:,j)=(diff(WBC,y,1).*WSB(row_num)+WBC.*diff(WSB(row_num),y,1)).*(diff(WBC,y,1).*WAB(j)+WBC.*diff(WAB(j),y,1));
	
	Vmax5bcC(:,j)=(diff(WBC,y,1).*WSC(row_num)+WBC.*diff(WSC(row_num),y,1)).*(diff(WBC,y,1).*WAC(j)+WBC.*diff(WAC(j),y,1));
	
	Vmax5bcD(:,j)=(diff(WBC,y,1).*WSD(row_num)+WBC.*diff(WSD(row_num),y,1)).*(diff(WBC,y,1).*WAD(j)+WBC.*diff(WAD(j),y,1));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax6bcA(:,j)=(diff(WBC,x,1).*WSA(row_num)+WBC.*diff(WSA(row_num),x,1)).*(diff(WBC,x,1).*WAA(j)+WBC.*diff(WAA(j),x,1));
	
	Vmax6bcB(:,j)=(diff(WBC,x,1).*WSB(row_num)+WBC.*diff(WSB(row_num),x,1)).*(diff(WBC,x,1).*WAB(j)+WBC.*diff(WAB(j),x,1));
	
	Vmax6bcC(:,j)=(diff(WBC,x,1).*WSC(row_num)+WBC.*diff(WSC(row_num),x,1)).*(diff(WBC,x,1).*WAC(j)+WBC.*diff(WAC(j),x,1));
	
	Vmax6bcD(:,j)=(diff(WBC,x,1).*WSD(row_num)+WBC.*diff(WSD(row_num),x,1)).*(diff(WBC,x,1).*WAD(j)+WBC.*diff(WAD(j),x,1));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Vmax7bcA(:,j)=(diff(WBC,x,1).*WSA(row_num)+WBC.*diff(WSA(row_num),x,1)).*(diff(WBC,y,1).*WAA(j)+WBC.*diff(WAA(j),y,1))+(diff(WBC,x,1).*WAA(j)+WBC.*diff(WAA(j),x,1)).*(diff(WBC,y,1).*WSA(row_num)+WBC.*diff(WSA(row_num),y,1));
	
	Vmax7bcB(:,j)=(diff(WBC,x,1).*WSB(row_num)+WBC.*diff(WSB(row_num),x,1)).*(diff(WBC,y,1).*WAB(j)+WBC.*diff(WAB(j),y,1))+(diff(WBC,x,1).*WAB(j)+WBC.*diff(WAB(j),x,1)).*(diff(WBC,y,1).*WSB(row_num)+WBC.*diff(WSB(row_num),y,1));
	
	Vmax7bcC(:,j)=(diff(WBC,x,1).*WSC(row_num)+WBC.*diff(WSC(row_num),x,1)).*(diff(WBC,y,1).*WAC(j)+WBC.*diff(WAC(j),y,1))+(diff(WBC,x,1).*WAC(j)+WBC.*diff(WAC(j),x,1)).*(diff(WBC,y,1).*WSC(row_num)+WBC.*diff(WSC(row_num),y,1));
	
	Vmax7bcD(:,j)=(diff(WBC,x,1).*WSD(row_num)+WBC.*diff(WSD(row_num),x,1)).*(diff(WBC,y,1).*WAD(j)+WBC.*diff(WAD(j),y,1))+(diff(WBC,x,1).*WAD(j)+WBC.*diff(WAD(j),x,1)).*(diff(WBC,y,1).*WSD(row_num)+WBC.*diff(WSD(row_num),y,1));
end	
for v=1:N*(N+3)/2
	for u=1:N*(N+3)/2
		%%%%%%%%%%%%%%%%%%% integration %%%%%%%%%%%%%%%%%%%
		func23_1A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1bcA(u,v))]);
		func23_1B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1bcB(u,v))]);
		func23_1C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1bcC(u,v))]);
		func23_1D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax1bcD(u,v))]);
		K23_1(u,v)=D*(integral2(func23_1A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_1B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_1C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_1D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		func23_2A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2bcA(u,v))]);
		func23_2B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2bcB(u,v))]);
		func23_2C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2bcC(u,v))]);
		func23_2D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax2bcD(u,v))]);
		K23_2(u,v)=D*(integral2(func23_2A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_2B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_2C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_2D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		func23_3A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3bcA(u,v))]);
		func23_3B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3bcB(u,v))]);
		func23_3C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3bcC(u,v))]);
		func23_3D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax3bcD(u,v))]);
		K23_3(u,v)=D*nu*(integral2(func23_3A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_3B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_3C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_3D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		func23_4A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4bcA(u,v))]);
		func23_4B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4bcB(u,v))]);
		func23_4C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4bcC(u,v))]);
		func23_4D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax4bcD(u,v))]);
		K23_4(u,v)=D*2*(1-nu)*(integral2(func23_4A_cell{u,v},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_4B_cell{u,v},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_4C_cell{u,v},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func23_4D_cell{u,v},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		func23_5A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5bcA(u,v))]);
		func23_5B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5bcB(u,v))]);
		func23_5C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5bcC(u,v))]);
		func23_5D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax5bcD(u,v))]);
		func23_6A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6bcA(u,v))]);
		func23_6B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6bcB(u,v))]);
		func23_6C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6bcC(u,v))]);
		func23_6D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax6bcD(u,v))]);
		func23_7A_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7bcA(u,v))]);
		func23_7B_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7bcB(u,v))]);
		func23_7C_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7bcC(u,v))]);
		func23_7D_cell{u,v}=eval(['@(x,y)',vectorize(Vmax7bcD(u,v))]);
		for j=1:num_domain
			if j==1 || j==2								%for area A1 and A2
				for i=1:num_triangles
					for q=1:num_intpoints
						K23_5(u,v)=K23_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func23_5A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K23_6(u,v)=K23_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func23_6A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K23_7(u,v)=K23_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func23_7A_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			elseif j==3 || j==4							%for area B1 and B2
				for i=1:num_triangles
					for q=1:num_intpoints
						K23_5(u,v)=K23_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func23_5B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K23_6(u,v)=K23_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func23_6B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K23_7(u,v)=K23_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func23_7B_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			elseif j==5 || j==6							%for area C1 and C2
				for i=1:num_triangles
					for q=1:num_intpoints
						K23_5(u,v)=K23_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func23_5C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K23_6(u,v)=K23_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func23_6C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K23_7(u,v)=K23_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func23_7C_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			else										%for area D1 and D2
				for i=1:num_triangles
					for q=1:num_intpoints
						K23_5(u,v)=K23_5(u,v)+h*intpoints_sigma(j,i,q,2)*intpoints{j,i,2}(q,1)*feval(func23_5D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K23_6(u,v)=K23_6(u,v)+h*intpoints_sigma(j,i,q,1)*intpoints{j,i,2}(q,1)*feval(func23_6D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
						K23_7(u,v)=K23_7(u,v)+h*intpoints_sigma(j,i,q,3)*intpoints{j,i,2}(q,1)*feval(func23_7D_cell{u,v},intpoints{j,i,1}(q,1),intpoints{j,i,1}(q,2));
					end
				end
			end
		end
	end
end
K23=K23_1+K23_2+K23_3+K23_4+K23_5+K23_6+K23_7;	
K(1:I^2,I^2+1:I^2+N*(N+3)/2)=K12;
K(1:I^2,I^2+N*(N+3)/2+1:I^2+N*(N+3))=K13;
K(I^2+1:I^2+N*(N+3)/2,I^2+N*(N+3)/2+1:I^2+N*(N+3))=K23;	
for i=I^2+1:I^2+N*(N+3)/2
	for j=1:I^2
		K(i,j)=K(j,i);
	end
end
for i=I^2+N*(N+3)/2+1:I^2+N*(N+3)
	for j=1:I^2+N*(N+3)/2
		K(i,j)=K(j,i);
	end
end
toc
tic
%%%%%%%%%%%%%%%%%%% mass matrix %%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%% diagonal matrices %%%%%%%%%%%%%%%%%%%		
M1=zeros(I^2,I^2);	
Tmaxaa=sym(zeros(I^2,I^2));
func_M1_cell=cell(I^2);
parfor v=1:I^2
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum kinetic energy function Tmax/Omega^2 (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to a(ij) and a(kl)
	Tmaxaa(:,v)=WBC.^2.*PhiX(row1).*PhiY(col1).*PhiX(row2(v)).*PhiY(col2(v));
end			
for v=1:I^2
	for u=1:I^2
		%%%%%%%%%%%%%%%%%%% integration %%%%%%%%%%%%%%%%%%%
		if size(symvar(Tmaxaa(u,v)),2)==2
			func_M1_cell{u,v}=eval(['@(x,y)',vectorize(Tmaxaa(u,v))]);
			M1(u,v)=Rho*h*integral2(func_M1_cell{u,v},0,a,0,b,'AbsTol',1e-12,'RelTol',1e-12);
		else
			M1(u,v)=Rho*h*int(int(Tmaxaa(u,v),x,0,a),y,0,b);
		end
	end
end
toc
tic
M2=zeros(N*(N+3)/2,N*(N+3)/2);
M3=zeros(N*(N+3)/2,N*(N+3)/2);
M=zeros(I^2+N*(N+3),I^2+N*(N+3));	
TmaxbbA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxbbB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxbbC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxbbD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxccA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxccB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxccC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxccD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
func_M2A_cell=cell(N*(N+3)/2);
func_M2B_cell=cell(N*(N+3)/2);
func_M2C_cell=cell(N*(N+3)/2);
func_M2D_cell=cell(N*(N+3)/2);
func_M3A_cell=cell(N*(N+3)/2);
func_M3B_cell=cell(N*(N+3)/2);
func_M3C_cell=cell(N*(N+3)/2);
func_M3D_cell=cell(N*(N+3)/2);
for j=1:N*(N+3)/2
	%%%%%%%%%%%%%%%%%%% differentiation with regard to b(i) and b(j)
	TmaxbbA(:,j)=WBC.^2.*WSA(row_num).*WSA(j);
	TmaxbbB(:,j)=WBC.^2.*WSB(row_num).*WSB(j);
	TmaxbbC(:,j)=WBC.^2.*WSC(row_num).*WSC(j);
	TmaxbbD(:,j)=WBC.^2.*WSD(row_num).*WSD(j);
	%%%%%%%%%%%%%%%%%%% differentiation with regard to c(i) and c(j)
	TmaxccA(:,j)=WBC.^2.*WAA(row_num).*WAA(j);
	TmaxccB(:,j)=WBC.^2.*WAB(row_num).*WAB(j);
	TmaxccC(:,j)=WBC.^2.*WAC(row_num).*WAC(j);
	TmaxccD(:,j)=WBC.^2.*WAD(row_num).*WAD(j);
end
for j=1:N*(N+3)/2
	for i=1:N*(N+3)/2		
		%%%%%%%%%%%%%%%%%%% integration %%%%%%%%%%%%%%%%%%%
		func_M2A_cell{i,j}=eval(['@(x,y)',vectorize(TmaxbbA(i,j))]);
		func_M2B_cell{i,j}=eval(['@(x,y)',vectorize(TmaxbbB(i,j))]);
		func_M2C_cell{i,j}=eval(['@(x,y)',vectorize(TmaxbbC(i,j))]);
		func_M2D_cell{i,j}=eval(['@(x,y)',vectorize(TmaxbbD(i,j))]);
		M2(i,j)=Rho*h*(integral2(func_M2A_cell{i,j},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M2B_cell{i,j},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M2C_cell{i,j},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M2D_cell{i,j},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		func_M3A_cell{i,j}=eval(['@(x,y)',vectorize(TmaxccA(i,j))]);
		func_M3B_cell{i,j}=eval(['@(x,y)',vectorize(TmaxccB(i,j))]);
		func_M3C_cell{i,j}=eval(['@(x,y)',vectorize(TmaxccC(i,j))]);
		func_M3D_cell{i,j}=eval(['@(x,y)',vectorize(TmaxccD(i,j))]);
		M3(i,j)=Rho*h*(integral2(func_M3A_cell{i,j},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M3B_cell{i,j},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M3C_cell{i,j},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M3D_cell{i,j},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
	end
end
M(1:I^2,1:I^2)=M1;
M(I^2+1:I^2+N*(N+3)/2,I^2+1:I^2+N*(N+3)/2)=M2;
M(I^2+N*(N+3)/2+1:I^2+N*(N+3),I^2+N*(N+3)/2+1:I^2+N*(N+3))=M3;
toc
tic
%%%%%%%%%%%%%%%%%%% mass matrix %%%%%%%%%%%%%%%%%%%		
%%%%%%%%%%%%%%%%%%% off-diagonal matrices %%%%%%%%%%%%%%%%%%%	
M12=zeros(I^2,N*(N+3)/2);
M13=zeros(I^2,N*(N+3)/2);
TmaxabA=sym(zeros(I^2,N*(N+3)/2));
TmaxabB=sym(zeros(I^2,N*(N+3)/2));
TmaxabC=sym(zeros(I^2,N*(N+3)/2));
TmaxabD=sym(zeros(I^2,N*(N+3)/2));
TmaxacA=sym(zeros(I^2,N*(N+3)/2));
TmaxacB=sym(zeros(I^2,N*(N+3)/2));
TmaxacC=sym(zeros(I^2,N*(N+3)/2));
TmaxacD=sym(zeros(I^2,N*(N+3)/2));
func_M12A_cell=cell(I^2,N*(N+3)/2);
func_M12B_cell=cell(I^2,N*(N+3)/2);
func_M12C_cell=cell(I^2,N*(N+3)/2);
func_M12D_cell=cell(I^2,N*(N+3)/2);
func_M13A_cell=cell(I^2,N*(N+3)/2);
func_M13B_cell=cell(I^2,N*(N+3)/2);
func_M13C_cell=cell(I^2,N*(N+3)/2);
func_M13D_cell=cell(I^2,N*(N+3)/2);
parfor k=1:N*(N+3)/2
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum kinetic energy function Tmax/Omega.^2 (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to a(ij) and b(k)
	TmaxabA(:,k)=WBC.^2.*PhiX(row1).*PhiY(col1).*WSA(k);
	
	TmaxabB(:,k)=WBC.^2.*PhiX(row1).*PhiY(col1).*WSB(k);
	
	TmaxabC(:,k)=WBC.^2.*PhiX(row1).*PhiY(col1).*WSC(k);
	
	TmaxabD(:,k)=WBC.^2.*PhiX(row1).*PhiY(col1).*WSD(k);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	TmaxacA(:,k)=WBC.^2.*PhiX(row1).*PhiY(col1).*WAA(k);
	
	TmaxacB(:,k)=WBC.^2.*PhiX(row1).*PhiY(col1).*WAB(k);
	
	TmaxacC(:,k)=WBC.^2.*PhiX(row1).*PhiY(col1).*WAC(k);
	
	TmaxacD(:,k)=WBC.^2.*PhiX(row1).*PhiY(col1).*WAD(k);
end
for k=1:N*(N+3)/2
	for u=1:I^2
		%%%%%%%%%%%%%%%%%%% integration %%%%%%%%%%%%%%%%%%%
		func_M12A_cell{u,k}=eval(['@(x,y)',vectorize(TmaxabA(u,k))]);
		func_M12B_cell{u,k}=eval(['@(x,y)',vectorize(TmaxabB(u,k))]);
		func_M12C_cell{u,k}=eval(['@(x,y)',vectorize(TmaxabC(u,k))]);
		func_M12D_cell{u,k}=eval(['@(x,y)',vectorize(TmaxabD(u,k))]);
		M12(u,k)=Rho*h*(integral2(func_M12A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M12B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M12C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M12D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
		func_M13A_cell{u,k}=eval(['@(x,y)',vectorize(TmaxacA(u,k))]);
		func_M13B_cell{u,k}=eval(['@(x,y)',vectorize(TmaxacB(u,k))]);
		func_M13C_cell{u,k}=eval(['@(x,y)',vectorize(TmaxacC(u,k))]);
		func_M13D_cell{u,k}=eval(['@(x,y)',vectorize(TmaxacD(u,k))]);
		M13(u,k)=Rho*h*(integral2(func_M13A_cell{u,k},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M13B_cell{u,k},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M13C_cell{u,k},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M13D_cell{u,k},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
	end
end
toc
tic
M23=zeros(N*(N+3)/2,N*(N+3)/2);
TmaxbcA=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxbcB=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxbcC=sym(zeros(N*(N+3)/2,N*(N+3)/2));
TmaxbcD=sym(zeros(N*(N+3)/2,N*(N+3)/2));
func_M23A_cell=cell(N*(N+3)/2);
func_M23B_cell=cell(N*(N+3)/2);
func_M23C_cell=cell(N*(N+3)/2);
func_M23D_cell=cell(N*(N+3)/2);
parfor j=1:N*(N+3)/2
	%%%%%%%%%%%%%%%%%%% the expression of the terms of maximum kinetic energy function Tmax/Omega.^2 (integrand) 
	%%%%%%%%%%%%%%%%%%% differentiation with regard to b(i) and c(j)
	TmaxbcA(:,j)=WBC.^2.*WSA(row_num).*WAA(j);
	
	TmaxbcB(:,j)=WBC.^2.*WSB(row_num).*WAB(j);
	
	TmaxbcC(:,j)=WBC.^2.*WSC(row_num).*WAC(j);
	
	TmaxbcD(:,j)=WBC.^2.*WSD(row_num).*WAD(j);
end
for j=1:N*(N+3)/2
	for i=1:N*(N+3)/2		
		%%%%%%%%%%%%%%%%%%% integration %%%%%%%%%%%%%%%%%%%
		func_M23A_cell{i,j}=eval(['@(x,y)',vectorize(TmaxbcA(i,j))]);
		func_M23B_cell{i,j}=eval(['@(x,y)',vectorize(TmaxbcB(i,j))]);
		func_M23C_cell{i,j}=eval(['@(x,y)',vectorize(TmaxbcC(i,j))]);
		func_M23D_cell{i,j}=eval(['@(x,y)',vectorize(TmaxbcD(i,j))]);
		M23(i,j)=Rho*h*(integral2(func_M23A_cell{i,j},0,x0,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M23B_cell{i,j},0,x0,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M23C_cell{i,j},x0,a,0,y_crackline,'AbsTol',1e-12,'RelTol',1e-12)+integral2(func_M23D_cell{i,j},x0,a,y_crackline,b,'AbsTol',1e-12,'RelTol',1e-12));
	end
end
M(1:I^2,I^2+1:I^2+N*(N+3)/2)=M12;
M(1:I^2,I^2+N*(N+3)/2+1:I^2+N*(N+3))=M13;
M(I^2+1:I^2+N*(N+3)/2,I^2+N*(N+3)/2+1:I^2+N*(N+3))=M23;	
for i=I^2+1:I^2+N*(N+3)/2
	for j=1:I^2
		M(i,j)=M(j,i);
	end
end
for i=I^2+N*(N+3)/2+1:I^2+N*(N+3)
	for j=1:I^2+N*(N+3)/2
		M(i,j)=M(j,i);
	end
end
toc
tic
X1=zeros(I^2+N*(N+3),I^2+N*(N+3));	
eigenvalues=zeros(I^2+N*(N+3),I^2+N*(N+3));
[X,eigenvalues]=eig(K,M);
Eigenvalues=diag(eigenvalues);
Hz=Eigenvalues.^0.5/(2*pi);
Hz_sort=sort(real(Hz));	
dimensionlessFreq_sort=Hz_sort*2*pi*a^2*(Rho*h/D)^0.5	%sorted dimensionless frequencies
dimensionlessFreq_no_sort=Hz*2*pi*a^2*(Rho*h/D)^0.5		%unsorted dimensionless frequencies
toc
% mailme('1441252861@qq.com','MATLAB calculation finished','MATLAB calculation finished'); %		
%%%%%%%%%%%%%%%%%%%%%%%%% mode shape calculation %%%%%%%%%%%%%%%%%%%%%%%%%%
tic
x_nodes=linspace(0,a,52);
y_nodes=linspace(0,b,52);
for MODE_NUMBER=1:I^2+N*(N+3)
	EIGENVECTOR=real(X(:,MODE_NUMBER))';
	for i=1:size(x_nodes,2)
		xx=x_nodes(i);	%x coordinate of the node 
		for j=1:size(y_nodes,2)
			yy=y_nodes(j);	%y coordinate of the node 
			for k=1:I
				SET_X(k)=xx^ll*(xx-a)^mm*xx^(k-1);
				SET_Y(k)=yy^nn*(yy-b)^qq*yy^(k-1);
			end
			counter_1=0;
			for n=1:I
				for m=1:I
					counter_1=counter_1+1;
					EVALUATED_FUNCTIONS_1(counter_1)=SET_X(n)*SET_Y(m);
				end
			end
			if xx<=x0 && yy>=(c+(xx-a)*tan(Alpha))	%if the node is in area A
				counter_2=0;
				for n=1:N
					for l=0:n
						counter_2=counter_2+1;
						EVALUATED_FUNCTIONS_2(counter_2)=(xx^ll*(xx-a)^mm*yy^nn*(yy-b)^qq)*((((xx-x0)^2+(yy-y0)^2)^0.5)^((2*n+1)/2)*cos((2*l+1)/2*(atan((yy-y0)/(xx-x0))-Alpha)));	%WBC*WSA(p)
					end
				end
				counter_3=0;
				for n=1:N
					for l=0:n
						counter_3=counter_3+1;
						EVALUATED_FUNCTIONS_3(counter_3)=(xx^ll*(xx-a)^mm*yy^nn*(yy-b)^qq)*((((xx-x0)^2+(yy-y0)^2)^0.5)^((2*n+1)/2)*sin((2*l+1)/2*(atan((yy-y0)/(xx-x0))-Alpha)));	%WBC*WAA(p)
					end
				end
			elseif xx<=x0 && yy<(c+(xx-a)*tan(Alpha))	%if the node is in area B
				counter_2=0;
				for n=1:N
					for l=0:n
						counter_2=counter_2+1;
						EVALUATED_FUNCTIONS_2(counter_2)=(xx^ll*(xx-a)^mm*yy^nn*(yy-b)^qq)*((((xx-x0)^2+(yy-y0)^2)^0.5)^((2*n+1)/2)*cos((2*l+1)/2*(atan((yy-y0)/(xx-x0))-Alpha)));
					end
				end
				counter_3=0;
				for n=1:N
					for l=0:n
						counter_3=counter_3+1;
						EVALUATED_FUNCTIONS_3(counter_3)=(xx^ll*(xx-a)^mm*yy^nn*(yy-b)^qq)*((((xx-x0)^2+(yy-y0)^2)^0.5)^((2*n+1)/2)*sin((2*l+1)/2*(atan((yy-y0)/(xx-x0))-Alpha)));
					end
				end
			elseif xx>x0 && yy<=(c+(xx-a)*tan(Alpha))	%if the node is in area C
				counter_2=0;
				for n=1:N
					for l=0:n
						counter_2=counter_2+1;
						EVALUATED_FUNCTIONS_2(counter_2)=(xx^ll*(xx-a)^mm*yy^nn*(yy-b)^qq)*((((xx-x0)^2+(yy-y0)^2)^0.5)^((2*n+1)/2)*cos((2*l+1)/2*(pi-Alpha+atan((yy-y0)/(xx-x0)))));
					end
				end
				counter_3=0;
				for n=1:N
					for l=0:n
						counter_3=counter_3+1;
						EVALUATED_FUNCTIONS_3(counter_3)=(xx^ll*(xx-a)^mm*yy^nn*(yy-b)^qq)*((((xx-x0)^2+(yy-y0)^2)^0.5)^((2*n+1)/2)*sin((2*l+1)/2*(pi-Alpha+atan((yy-y0)/(xx-x0)))));
					end
				end
			else										%if the node is in area D
				counter_2=0;
				for n=1:N
					for l=0:n
						counter_2=counter_2+1;
						EVALUATED_FUNCTIONS_2(counter_2)=(xx^ll*(xx-a)^mm*yy^nn*(yy-b)^qq)*((((xx-x0)^2+(yy-y0)^2)^0.5)^((2*n+1)/2)*cos((2*l+1)/2*(-pi-Alpha+atan((yy-y0)/(xx-x0)))));
					end
				end
				counter_3=0;
				for n=1:N
					for l=0:n
						counter_3=counter_3+1;
						EVALUATED_FUNCTIONS_3(counter_3)=(xx^ll*(xx-a)^mm*yy^nn*(yy-b)^qq)*((((xx-x0)^2+(yy-y0)^2)^0.5)^((2*n+1)/2)*sin((2*l+1)/2*(-pi-Alpha+atan((yy-y0)/(xx-x0)))));
					end
				end
			end
			EVALUATED_FUNCTIONS=[EVALUATED_FUNCTIONS_1, double(EVALUATED_FUNCTIONS_2), double(EVALUATED_FUNCTIONS_3)];
			EVALUATED_MODE(i,j)=-dot(EVALUATED_FUNCTIONS,EIGENVECTOR);
		end
	end
	figure(MODE_NUMBER)
	surf(x_nodes,y_nodes,EVALUATED_MODE)
	view(60,30);
end
toc




























