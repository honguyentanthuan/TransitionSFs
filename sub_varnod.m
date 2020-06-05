function [D_node,D_elem]=sub_varnod(var_nod,ite)
% =========================================================================
% Coded by : Thuan Ho-Nguyen-Tan
% Email    : honguyentanthuan@seoultech.ac.kr
%
% Computational Solid Mechanics Laboratory,
% Department of Mechanical and Automotive Engineering,
% Seoul National University of Science and Technology (SeoulTech),Korea.
% =========================================================================
x1=[-1 -1];
x2=[1 -1];
x3=[1 1];
x4=[-1 1];

n1=var_nod(1);
n2=var_nod(2);
n3=var_nod(3);
n4=var_nod(4);

if n1==0
    coord1=[];
else
    k1=n1+1;
    for i=1:n1
        t=i/k1;
        coord1(i,:)=x1*(1-t)+x2*t;
    end
end

if n2==0
    coord2=[];
else
    k2=n2+1;
    for i=1:n2
        t=i/k2;
        coord2(i,:)=x2*(1-t)+x3*t;
    end
end

if n3==0
    coord3=[];
else
    k3=n3+1;
    for i=1:n3
        t=i/k3;
        coord3(i,:)=x3*(1-t)+x4*t;
    end
end

if n4==0
    coord4=[];
else
    k4=n4+1;
    for i=1:n4
        t=i/k4;
        coord4(i,:)=x4*(1-t)+x1*t;
    end
end

% create sub-domains
if n1>n3
    temp1=2/(n1+1);
    numx=n1+1;
elseif n3>n1
    temp1=2/(n3+1);
    numx=n3+1;
elseif n3==n1  && n1~=0
    temp1=2/(n3+1);
    numx=n1+1;
else
    temp1=2;
    numx=1;
end

if n2>n4
    temp2=2/(n2+1);
    numy=n2+1;
elseif n4>n2
    temp2=2/(n4+1);
    numy=n4+1;
elseif n4==n2  && n4~=0
    temp2=2/(n4+1);
    numy=n4+1;
else
    temp2=2;
    numy=1;
end

[xvar,yvar]=meshgrid(-1:temp1:1,-1:temp2:1);
xvar1=reshape(xvar',size(xvar,1)*size(yvar,2),1);
yvar1=reshape(yvar',size(xvar,1)*size(yvar,2),1);

temp = 0;

for i = 1: numx
    for j = 1:numy
        
        sw = i + ( j - 1 ) * ( numx + 1 );
        se = i + 1 + ( j - 1 ) * ( numx + 1 );
        nw = i + j * ( numx + 1 );
        ne = i + 1 + j * ( numx + 1 );
        
        temp = temp + 1;
        
        D_elem1(1) = sw;
        D_elem1(2) = se;
        D_elem1(3) = ne;
        D_elem1(4) = nw;
        
        D_elem{temp,:}=D_elem1;
    end
end
D_node=[xvar1 yvar1];

 % --- refine mesh
 for i=1:ite
     [D_node,D_elem]=create_rectangle(D_node,D_elem);
 end


