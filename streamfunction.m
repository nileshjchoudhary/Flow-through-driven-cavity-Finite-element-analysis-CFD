function [ streamline ] =streamfunction(  uvp ,nX,nY,initialX,initialY,NRe,finalX,finalY,tol, X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode,Qnnod,Lnnod ,omega,streamline )
[ w ,gp ] = gausspoints;



Res=zeros(tQnode,1);
jac=zeros(tQnode);
uy=1;
Luy=1;
for iY=1:nY
    ux=1;
    Lux=1;
    for iX=1:nX
        [ TrkNode ] = nodetrack( ux,uy,nX,nY,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod,Lux,Luy);
    
  for  j=1:9
for k=1:9

     for r=1:3 
     for s=1:3 
 zai=gp(s);
 eta=gp(r);

[ phi ,dphie ,dphiz ] = biquadBF( zai ,eta);
 [ oneDphiz,oneDdphizz, oneDphie,oneDdphiee] = oneDbasis( zai,eta );

  [ sai  ] = bilinearBF( zai ,eta);
  
  



x_z=X(ux)*oneDdphizz(1)+X(ux+1)*oneDdphizz(2)+X(ux+2)*oneDdphizz(3);
y_e=Y(uy)*oneDdphiee(1)+Y(uy+1)*oneDdphiee(2)+Y(uy+2)*oneDdphiee(3);
 [  U_z,V_e,U_e,Ui,V_z,Vi,P] = veloVec( TrkNode, dphiz,dphie,phi,sai,uvp );

          
                    det=x_z*y_e;
                    
                    if k==1
                        
                        Res(TrkNode(1,j))=Res(TrkNode(1,j))+(-(omega(TrkNode(1,j))*phi(j))*det)*w(r)*w(s);
                        
                    end
                    
                    jac(TrkNode(1,j),TrkNode(1,k))=jac(TrkNode(1,j),TrkNode(1,k))...
                        +((dphiz(j)*dphiz(k)*(1/x_z^2)+dphie(j)*dphie(k)*(1/y_e^2))*det)*w(s)*w(r);
                    
                 
                 


               end
     end
    
end
 end
ux=ux+2;
Lux=Lux+1;
    end
uy=uy+2;  
Luy=Luy+1;
end
temp6=zeros(2*QnodeX+2*QnodeY-4,1);

temp6(1:QnodeX,1)=Qnnod(1,1:QnodeX);
 temp6(QnodeX+1:2*QnodeX,1)=Qnnod(QnodeY,1:QnodeX);  
  temp6(2*QnodeX+1:2*QnodeX+QnodeY-2,1)=Qnnod(2:QnodeY-1,1); 
  temp6(2*QnodeX+QnodeY-1:2*QnodeX+2*QnodeY-4,1)=Qnnod(2:QnodeY-1,QnodeX); 
  for t6=1:2*QnodeX+2*QnodeY-4
      jac(temp6(t6),:)=0;
      jac(temp6(t6),temp6(t6))=1;
      Res(temp6(t6),1)=0;
  end
  streamline=jac\Res;
end



