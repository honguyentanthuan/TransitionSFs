% =========================================================================
% 
%                 FEM TRANSITION SHAPE FUNCTIONS
%                     based on rectagle element
%                       Finite Element Program
% 
% Coded by : Thuan Ho-Nguyen-Tan
% Email    : honguyentanthuan@seoultech.ac.kr
% 
% Computational Solid Mechanics Laboratory,
% Department of Mechanical and Automotive Engineering,
% Seoul National University of Science and Technology (SeoulTech),Korea.
% =========================================================================
%   Last modify: 2016.Dec.23 
% -------------------------------------------------------------------------

clc
clear all
close all

% --- rectangle master element 
nodes = [-1 -1;1 -1;1 1;-1 1];

% --- input number of variable nodes 
lkmn = [1,1,1,1];
n = 4+ lkmn(1) +lkmn(2) +lkmn(3) +lkmn(4) ;

% --- gauss points
ite=2;
[pospg,nodes_sub]=sub_varnod(lkmn,ite);

% --- construct MLS shape functions
xyzGauss =zeros(size(pospg,1),2);
psi =zeros(size(pospg,1),n);
dpsi_xi =zeros(size(pospg,1),n);
dpsi_eta =zeros(size(pospg,1),n);
  

for igauss=1:size(pospg,1)
    xyzGauss(igauss,:)=pospg(igauss,:);
    
    [psi(igauss,:),dpsi_xi(igauss,:),dpsi_eta(igauss,:)]= shape_transition(...
       xyzGauss(igauss,:),lkmn);
   
end

% --- plot shape functions
gcoord_sub=xyzGauss;
for i=1:n
    figure()
    gcoord_sub(:,3)= psi(:,i);
    plot_mesh(gcoord_sub,nodes_sub,0);
end

% -------------------------------------------------------------------------
