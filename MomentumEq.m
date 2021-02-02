function  [ Residu, jacobi ] =MomentumEq( Residu, jacobi,uvp,nX,nY,initialX,initialY,NRe,finalX,finalY,tol,X,Y,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod )
[ w ,gp ] = gausspoints;

uy=1;
Luy=1;
for iY=1:nY
    ux=1;
    Lux=1;
    for iX=1:nX
        [ TrkNode ] = nodetrack( ux,uy,nX,nY,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod,Lux,Luy);
    
  for  j=1:9
 for  k=1:9;

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

                 if k==1; 
Residu(TrkNode(1,j))=Residu(TrkNode(1,j))...
                    +((phi(j)*NRe*(Ui*(1/x_z)*U_z+Vi*(1/y_e)*U_e)...
                    +(-P+2*(1/x_z)*U_z)*(1/x_z)*dphiz(j)...
                    +((1/y_e)*U_e+(1/x_z)*V_z)*(1/y_e)*dphie(j))*x_z*y_e)*w(r)*w(s);
                       
Residu(TrkNode(2,j))=Residu(TrkNode(2,j))...
                    +((phi(j)*NRe*(Ui*(1/x_z)*V_z+Vi*(1/y_e)*V_e)...
                    +(-P+2*(1/y_e)*V_e)*(1/y_e)*dphie(j)...
                    +((1/y_e)*U_e+(1/x_z)*V_z)*(1/x_z)*dphiz(j))*x_z*y_e)*w(r)*w(s);      

                                                  
 
                 end

 jacobi(TrkNode(1,j),TrkNode(1,k))=jacobi(TrkNode(1,j),TrkNode(1,k))...
                    +((NRe*phi(j)*(phi(k)*(1/x_z)*U_z+Ui*(1/x_z)*dphiz(k)...
                    + Vi*(1/y_e)*dphie(k))+2*(1/x_z^2)*dphiz(j)*dphiz(k)...
                    +(1/y_e^2)*dphie(j)*dphie(k))*x_z*y_e)*w(r)*w(s);
                             
 jacobi(TrkNode(1,j),TrkNode(2,k))=jacobi(TrkNode(1,j),TrkNode(2,k))...
     +((NRe*(1/y_e)*U_e*phi(j)*phi(k)+(1/y_e)*dphie(j)*(1/x_z)*dphiz(k))*x_z*y_e)*w(r)*w(s);

                 
                 
 jacobi(TrkNode(2,j),TrkNode(1,k))=jacobi(TrkNode(2,j),TrkNode(1,k))...
     +((NRe*(1/x_z)*V_z*phi(j)*phi(k)+(1/y_e)*dphie(k)*(1/x_z)*dphiz(j))*x_z*y_e)*w(r)*w(s);
                 
 jacobi(TrkNode(2,j),TrkNode(2,k))=jacobi(TrkNode(2,j),TrkNode(2,k))...
                     +((NRe*phi(j)*(dphiz(k)*(1/x_z)*Ui+(1/y_e)*phi(k)*V_e+Vi*(1/y_e)*dphie(k))...
                     +2*(1/y_e^2)*dphie(j)*dphie(k)+(1/x_z^2)*dphiz(j)*dphiz(k))*x_z*y_e)*w(r)*w(s); 



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
   
end
