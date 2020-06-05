function [coord1,coord2,coord3,coord4] = making_variable_node_elem(lkmn,nodes)
% =========================================================================
% Coded by : Thuan Ho-Nguyen-Tan
% Email    : honguyentanthuan@seoultech.ac.kr
% 
% Computational Solid Mechanics Laboratory,
% Department of Mechanical and Automotive Engineering,
% Seoul National University of Science and Technology (SeoulTech),Korea.
% =========================================================================
% --- initialize
x1 = nodes(1,:);
x2 = nodes(2,:);
x3 = nodes(3,:);
x4 = nodes(4,:);

n1 = lkmn(1);
n2 = lkmn(2);
n3 = lkmn(3);
n4 = lkmn(4);

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

