function [uvp]=Driver(uvp,nX,nY,initialX,initialY,NRe,finalX,finalY,tol,X,Y,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod)

error=inf;
E1=error;E2=error;

itr=1;
while E1>tol || E2>tol
    Residu=zeros(2*tQnode+tLnode,1);
    jacobi=zeros(2*tQnode+tLnode);
    
    [ Residu, jacobi ] =ContinuityEq( Residu, jacobi,uvp,nX,nY,initialX,initialY,NRe,finalX,finalY,tol,X,Y,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod ,Lnnod );
    
    [ Residu, jacobi ] =MomentumEq( Residu, jacobi,uvp,nX,nY,initialX,initialY,NRe,finalX,finalY,tol,X,Y,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod );
   
[ jacobi ] = jacobiwrtbdry( Qnnod, jacobi,tQnode,QnodeX,QnodeY,nX,nY );
[ Residu] = residuwrtbdry( Residu,Qnnod ,QnodeX,QnodeY,nX,nY ,tQnode );




a=jacobi;
b=Residu;


[ duvp ] = a\b; 

E1=0;
E2=0;
dr(1:2*tQnode,1)=(Residu(1:2*tQnode,1)).*(Residu(1:2*tQnode,1));
du(1:2*tQnode,1)=(duvp(1:2*tQnode,1)).*(duvp(1:2*tQnode,1));
error1=sum(dr);
error2=sum(du);
E1=(error1^0.5)/(2*tQnode);
E2=(error2^0.5)/(2*tQnode);
%*********************************************************************

uvp=uvp+duvp;
itr=itr+1;
E1
E2

iteration =itr-1;
iteration
end
fprintf('final errors \n major calculation over wait few second');
end
