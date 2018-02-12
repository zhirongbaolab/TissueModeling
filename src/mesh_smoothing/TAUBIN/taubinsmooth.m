function [ Vnieuw ] = taubinsmooth( F,V,iterations,lambda,mu )


%TAUBINSMOOTH performs a back and forward Laplacian smoothing "without shrinking" of a triangulated mesh,as described
%by Gabriel Taubin (ICCV '95)

%INPUT
% - F m*3 faces list
% - V n*3 vertices list
% - iterations (default 10)
% - lamba and mu smoothing variables 0 <lambda <mu <1 (default 0.5 and
% 0.53)

%OUTPUT
% Vnieuw : smoothnend vertices

%EXAMPLE
%load example.mat
%[ Vnieuw ] = taubinsmooth( F,VT,30);
% VT(:,1)=VT(:,1)+30;
%clf
%trisurf(F,VT(:,1),VT(:,2),VT(:,3),'FaceColor','b');
%hold
%daspect([1 1 1]);
%light
%trisurf(F,Vnieuw(:,1),Vnieuw(:,2),Vnieuw(:,3),'FaceColor','y');


if(nargin<3), iterations=10; end
if(nargin<4), lambda=0.5; end
if(nargin<5), mu=0.53; end


[Ne]=neigboursTRI2(V,F);
% Ne=vertex_neighbours_double(F,V);
Voriginal=V;
[Indices_edges]=detectedges(V,F);

Vnieuw=V;

for k=1:iterations
    V=Vnieuw;
for i=1:size(V,1)
       indices=[Ne{i}]';    
    distances=sqrt(sum((V(indices,:)-repmat(V(i,:),size(indices,1),1)).^2,2));
    weights=distances.^-1;   
    vectoren=[weights.*V(indices,1) weights.*V(indices,2) weights.*V(indices,3)];
    totaldist=sum(weights);  
    Vnieuw(i,:)=V(i,:)+lambda*(sum(vectoren)./totaldist-V(i,:)); 
end

Vnieuw(Indices_edges,:)=V(Indices_edges,:);
V=Vnieuw;
for i=1:size(V,1)
    indices=[Ne{i}]';
    distances=sqrt(sum((V(indices,:)-repmat(V(i,:),size(indices,1),1)).^2,2));
    weights=distances.^-1;
    vectoren=[weights.*V(indices,1) weights.*V(indices,2) weights.*V(indices,3)];
    totaldist=sum(weights);
    Vnieuw(i,:)=V(i,:)-mu*(sum(vectoren)./totaldist-V(i,:));
end
Vnieuw(Indices_edges,:)=V(Indices_edges,:);
V=Vnieuw;

end

TFind = find(isnan(Vnieuw(:,1)));
Vnieuw(TFind,:)=Voriginal(TFind,:);


function[Indices_edges]=detectedges(V,F)

fk1 = F(:,1);
fk2 = F(:,2);
fk3 = F(:,3);

ed1=sort([fk1 fk2 ]')';
ed2=sort([fk1 fk3 ]')';
ed3=sort([fk2 fk3 ]')';

%single edges
ed=[ed1 ;ed2 ;ed3];
[etemp1,ia,ic]=unique(ed,'rows','stable');
esingle=ed(ia,:);

%dubbles
edouble=removerows(ed,ia);

C = setdiff(esingle,edouble,'rows');

Indices_edges=reshape(C,size(C,1)*2,1);



function [Ne]=neigboursTRI2(V,F)

Ne=cell(1,size(V,1));

for i=1:length(F)
    Ne{F(i,1)}=[Ne{F(i,1)} [F(i,2) F(i,3)]];
    Ne{F(i,2)}=[Ne{F(i,2)} [F(i,3) F(i,1)]];
    Ne{F(i,3)}=[Ne{F(i,3)} [F(i,1) F(i,2)]];
end

for i=1:size(V,1)
    Ne{i}=unique(Ne{i});
    if isempty (Ne{i})==1
        Ne{i}=[];
    end
    
end




