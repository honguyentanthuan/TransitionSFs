function [gcoord,elem]=create_rectangle(coord,ele)

gcoord = coord;
elem = ele;
for i=1:size(elem,1);
    
    poly = gcoord(elem{i,:},:);
    centroid = sum(poly,1)/4;
    
    m=zeros(5,2);
    m(1,:)=poly(1,:);
    m(2,:)=poly(2,:);
    m(3,:)=poly(3,:);
    m(4,:)=poly(4,:);
    m(5,:)=centroid;
    m(6,:)=(m(1,:)+m(2,:))/2;
    m(7,:)=(m(2,:)+m(3,:))/2;
    m(8,:)=(m(3,:)+m(4,:))/2;
    m(9,:)=(m(1,:)+m(4,:))/2;
    
    for h=1:size(m,1)
        temp=0;
        for e=1:size(gcoord,1)
            if m(h,:)==gcoord(e,:)
                id1(h,:)=e;
                temp=1;
            end
        end
        
        if temp==0
            gcoord=[gcoord;m(h,:)];
            id1(h,:)=size(gcoord,1);
        end
    end
    
    k=size(elem,1);
    
    elem{i,:}  =[id1(5) id1(9) id1(1) id1(6)];
    elem{k+1,:}=[id1(5) id1(6) id1(2) id1(7)];
    elem{k+2,:}=[id1(5) id1(7) id1(3) id1(8)];
    elem{k+3,:}=[id1(5) id1(8) id1(4) id1(9)];
    
end

