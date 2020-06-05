function [psi,dpsi_xi,dpsi_eta] = shape_transition(pospg,lkmn)
%==========================================================================
%     determine the values of the 2D shape functions on the Gauss points
%                         for rectangle element
%==========================================================================
% Coded by: Thuan Ho-Nguyen-Tan
% Email: honguyentanthuan@gmail.com
%
% Computational Solid Mechanics Laboratory,
% Department of Mechanical and Automotive Engineering,
% Seoul National University of Science and Technology (SeoulTech), Korea.
% Last modified: 7 July 2016
%--------------------------------------------------------------------------
n1 = lkmn(1);
n2 = lkmn(2);
n3 = lkmn(3);
n4 = lkmn(4);
nod = 4+n1+n2+n3+n4;

% --- initialize
psi = zeros(1,nod);
dpsi_xi = zeros(1,nod);
dpsi_eta = zeros(1,nod);
xi = pospg(:,1);
eta = pospg(:,2); 


% --- making variable-node element
nodes = [-1 -1;1 -1;1 1;-1 1];
[coord1,coord2,coord3,coord4]=making_variable_node_elem(lkmn,nodes);

%---- shape function and derivatives
psi(:,1) = (1-xi).*(1-eta)/4;
psi(:,2) = (1+xi).*(1-eta)/4;
psi(:,3) = (1+xi).*(1+eta)/4;
psi(:,4) = (1-xi).*(1+eta)/4;

dpsi_xi(:,1) = -(1-eta)/4;
dpsi_xi(:,2) = (1-eta)/4;
dpsi_xi(:,3) = (1+eta)/4;
dpsi_xi(:,4) = -(1+eta)/4;

dpsi_eta(:,1) = -(1-xi)/4;
dpsi_eta(:,2) = -(1+xi)/4;
dpsi_eta(:,3) = (1+xi)/4;
dpsi_eta(:,4) = (1-xi)/4;

% --- bottom side
 if n1>0 
     co = [-1 -1;coord1;1 -1];
     for i=2:size(co,1)-1
         x1 = co(i-1,1);
         x2 = co(i+1,1);
         xi_i = co(i,1);
         h = abs(xi_i-x1);
       
         if xi>=x1 && xi<= x2  
             psi(:,4+i-1) = (1-abs(xi-xi_i)/h).*(1-eta)/2;
             if xi~=xi_i
                 dpsi_xi(:,4+i-1) = -(xi-xi_i)/(abs(xi-xi_i)*h).*(1-eta)/2;
                 dpsi_eta(:,4+i-1) = -(1-abs(xi-xi_i)/h)/2;
             end
         end
     end
     for i=1:n1
         psi(:,1) = psi(:,1)-psi(:,4+i)*(1-(i/(n1+1)));
         psi(:,2) = psi(:,2)-psi(:,4+i)*i/(n1+1);
         
         dpsi_xi(:,1) = dpsi_xi(:,1)- dpsi_xi(:,4+i)*(1-(i/(n1+1)));
         dpsi_eta(:,1) = dpsi_eta(:,1) - dpsi_eta(:,4+i)*(1-(i/(n1+1)));
         dpsi_xi(:,2) = dpsi_xi(:,2)- dpsi_xi(:,4+i)*i/(n1+1);
         dpsi_eta(:,2) = dpsi_eta(:,2) - dpsi_eta(:,4+i)*i/(n1+1);
     end
 end  

 % --- right side
 if n2>0 
     co = [1 -1;coord2;1 1];
     for i=2:size(co,1)-1
         y1 = co(i-1,2);
         y2 = co(i+1,2);
         yi_i = co(i,2);
         h = abs(yi_i-y1);
       
         if eta>=y1 && eta<= y2 
             psi(:,4+n1+i-1) = (1-abs(eta-yi_i)/h).*(1+xi)/2;
             if  eta~=yi_i
                 dpsi_xi(:,4+n1+i-1) = (1-abs(eta-yi_i)/h)/2;
                 dpsi_eta(:,4+n1+i-1) = -(eta-yi_i)/(abs(eta-yi_i)*h).*(1+xi)/2;
             end
         end
     end
     for i=1:n2
         psi(:,2) = psi(:,2)-psi(:,4+n1+i)*(1-(i/(n2+1)));
         psi(:,3) = psi(:,3)-psi(:,4+n1+i)*i/(n2+1);
         
         dpsi_xi(:,2) = dpsi_xi(:,2)-dpsi_xi(:,4+n1+i)*(1-(i/(n2+1)));
         dpsi_eta(:,2) = dpsi_eta(:,2)-dpsi_eta(:,4+n1+i)*(1-(i/(n2+1)));
         dpsi_xi(:,3) = dpsi_xi(:,3)-dpsi_xi(:,4+n1+i)*i/(n2+1);
         dpsi_eta(:,3) = dpsi_eta(:,3)-dpsi_eta(:,4+n1+i)*i/(n2+1);
     end
 end  

% --- top side
 if n3>0 
     co = [1 1;coord3;-1 1];
     for i=2:size(co,1)-1
         x1 = co(i-1,1);
         x2 = co(i+1,1);
         xi_i = co(i,1);
         h = abs(xi_i-x1);
       
         if xi<=x1 && xi>= x2 
             psi(:,4+n1+n2+i-1) = (1-abs(xi-xi_i)/h).*(1+eta)/2;
             if xi~=xi_i
                 dpsi_xi(:,4+n1+n2+i-1) = -(xi-xi_i)/(abs(xi-xi_i)*h).*(1+eta)/2;
                 dpsi_eta(:,4+n1+n2+i-1) = (1-abs(xi-xi_i)/h)/2;
             end
         end
     end
   
     for i=1:n3
         psi(:,3) = psi(:,3)-psi(:,4+n1+n2+i)*(1-(i/(n3+1)));
         psi(:,4) = psi(:,4)-psi(:,4+n1+n2+i)*i/(n3+1);
         
         dpsi_xi(:,3) = dpsi_xi(:,3)-dpsi_xi(:,4+n1+n2+i)*(1-(i/(n3+1)));
         dpsi_eta(:,3) = dpsi_eta(:,3)-dpsi_eta(:,4+n1+n2+i)*(1-(i/(n3+1)));
         dpsi_xi(:,4) = dpsi_xi(:,4)-dpsi_xi(:,4+n1+n2+i)*i/(n3+1);
         dpsi_eta(:,4) = dpsi_eta(:,4)-dpsi_eta(:,4+n1+n2+i)*i/(n3+1);
     end
 end
 
  % --- left side
 if n4>0 
     co = [-1 1;coord4;-1 -1];
     for i=2:size(co,1)-1
         y1 = co(i-1,2);
         y2 = co(i+1,2);
         yi_i = co(i,2);
         h = abs(yi_i-y1);
       
         if eta<=y1 && eta>= y2
             psi(:,4+n1+n2+n3+i-1) = (1-abs(eta-yi_i)/h).*(1-xi)/2;
             if eta~=yi_i
                 dpsi_xi(:,4+n1+n2+n3+i-1) = -(1-abs(eta-yi_i)/h)/2;
                 dpsi_eta(:,4+n1+n2+n3+i-1) = -(eta-yi_i)/(abs(eta-yi_i)*h).*(1-xi)/2;
             end
         end
     end
     for i=1:n4
         psi(:,4) = psi(:,4)-psi(:,4+n1+n2+n3+i)*(1-(i/(n4+1)));
         psi(:,1) = psi(:,1)-psi(:,4+n1+n2+n3+i)*i/(n4+1);
         
         dpsi_xi(:,4) = dpsi_xi(:,4)-dpsi_xi(:,4+n1+n2+n3+i)*(1-(i/(n4+1)));
         dpsi_eta(:,4) = dpsi_eta(:,4)-dpsi_eta(:,4+n1+n2+n3+i)*(1-(i/(n4+1)));
         dpsi_xi(:,1) = dpsi_xi(:,1)-dpsi_xi(:,4+n1+n2+n3+i)*i/(n4+1);
         dpsi_eta(:,1) = dpsi_eta(:,1)-dpsi_eta(:,4+n1+n2+n3+i)*i/(n4+1);
     end
 end
 
 % --- arrange number order
 psi_n1 =zeros(1,n1);
 dpsi_xi_n1 =zeros(1,n1);
 dpsi_eta_n1 =zeros(1,n1);
 for i=1:n1
     psi_n1(:,i) = psi(:,4+i);
     dpsi_xi_n1(:,i) = dpsi_xi(:,4+i);
     dpsi_eta_n1(:,i) = dpsi_eta(:,4+i);
 end

 psi_n2 =zeros(1,n2);
 dpsi_xi_n2 =zeros(1,n2);
 dpsi_eta_n2 =zeros(1,n2);
 for i=1:n2
     psi_n2(:,i) = psi(:,4+n1+i);
     dpsi_xi_n2(:,i) = dpsi_xi(:,4+n1+i);
     dpsi_eta_n2(:,i) = dpsi_eta(:,4+n1+i);
 end
 
 psi_n3 =zeros(1,n3);
 dpsi_xi_n3 =zeros(1,n3);
 dpsi_eta_n3 =zeros(1,n3);
for i=1:n3
     psi_n3(:,i) = psi(:,4+n1+n2+i); 
     dpsi_xi_n3(:,i) = dpsi_xi(:,4+n1+n2+i);
     dpsi_eta_n3(:,i) = dpsi_eta(:,4+n1+n2+i);
end 

 psi_n4 =zeros(1,n4);
 dpsi_xi_n4 =zeros(1,n4);
 dpsi_eta_n4 =zeros(1,n4);
for i=1:n4
     psi_n4(:,i) = psi(:,4+n1+n2+n3+i);
     dpsi_xi_n4(:,i) = dpsi_xi(:,4+n1+n2+n3+i);
     dpsi_eta_n4(:,i) = dpsi_eta(:,4+n1+n2+n3+i);
end

psi = [psi(:,3) psi_n3 psi(:,4) psi_n4 psi(:,1) psi_n1 psi(:,2) psi_n2];
dpsi_xi = [dpsi_xi(:,3) dpsi_xi_n3 dpsi_xi(:,4) dpsi_xi_n4 dpsi_xi(:,1) dpsi_xi_n1 dpsi_xi(:,2) dpsi_xi_n2];
dpsi_eta = [dpsi_eta(:,3) dpsi_eta_n3 dpsi_eta(:,4) dpsi_eta_n4 dpsi_eta(:,1) dpsi_eta_n1 dpsi_eta(:,2) dpsi_eta_n2];



 
 