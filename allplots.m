function [ Ufinal,Vfinal ] = allplots( uvp,QnodeX,QnodeY ,tQnode,tLnode,nX,nY,X,Y,omega,streamline)

U=uvp(1:tQnode,1);
U=reshape(U,QnodeX,QnodeY);
U=U';
Ufinal=U;

omega=reshape(omega,QnodeX,QnodeY);
omega=omega';

streamline=reshape(streamline,QnodeX,QnodeY);
streamline=streamline';



V=uvp(tQnode+1:2*tQnode,1);
V=reshape(V,QnodeX,QnodeY);
V=V';
Vfinal=V;

P=uvp(2*tQnode+1:2*tQnode+tLnode,1);
P=reshape(P,nX+1,nY+1);
P=P';


figure(1)%%%%%%%%%%%%%%%%
quiver(X,Y,U,V);

figure(2)
contour(P,150);

figure(3)%%%%%%%%%%%%%%%%
contour(X,Y,omega,75)



figure(4)%%%%%%%%%%%%%%%%

contour(X,Y,streamline,75)


figure


surf(X,Y,streamline);

figure%%%%%%%%%%%%%%%%%
surf(P);




figure
contour(U,150);

figure
contour(V,150);
end

