function [ X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode,Qnnod,Lnnod] =notify( nX,nY,initialX,initialY,finalX,finalY )
QnodeX=2*nX+1;
QnodeY=2*nY+1;
LnodeX=nX+1;
LnodeY=nY+1;
tQnode=QnodeX*QnodeY;
tLnode=LnodeX*LnodeY;


X=linspace(initialX,finalX,QnodeX);
Y=linspace(initialY,finalY,QnodeY);
Qnnod=linspace(1,tQnode,tQnode);
Qnnod=reshape(Qnnod,QnodeX,QnodeY);
Qnnod=Qnnod';
Lnnod=linspace(1,tLnode,tLnode);
Lnnod=reshape(Lnnod,LnodeX,LnodeY);
Lnnod=Lnnod';
        
end

