function [ uvp ,omega, streamline] = BoundaryCondition( X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode)

U=zeros(QnodeY,QnodeX);
V=zeros(QnodeY,QnodeX);
P=ones(LnodeY,LnodeX);
omega=zeros(QnodeY,QnodeX);
streamline=zeros(QnodeY,QnodeX);
%P(:,1)=5;
 U(QnodeY,:)=1; 
% U(1,:)=1;
u=U';
u=u(:);
v=V';
v=V(:);
p=P';
p=p(:);
omega=omega';
omega=omega(:);
streamline=streamline';
streamline=streamline(:);
uvp=[u;v;p];
end

