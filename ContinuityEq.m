function  [ Residu, jacobi ] =ContinuityEq( Residu, jacobi,uvp,nX,nY,initialX,initialY,NRe,finalX,finalY,tol,X,Y,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod )
[ w ,gp ] = gausspoints;

uy=1;
Luy=1;
for iY=1:nY
    ux=1;
    Lux=1;
    for iX=1:nX
        [ TrkNode ] = nodetrack( ux,uy,nX,nY,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod,Lux,Luy);
    
 for  j=1:4
 for  k=1:9;

 for r=1:3 
 for s=1:3 

 zai=gp(s);
 eta=gp(r);

[ phi ,dphie ,dphiz ] = biquadBF( zai ,eta);
 [ oneDphiz,oneDdphizz, oneDphie,oneDdphiee] = oneDbasis( zai,eta );

  [ sai ] = bilinearBF( zai ,eta);
  
x_z=X(ux)*oneDdphizz(1)+X(ux+1)*oneDdphizz(2)+X(ux+2)*oneDdphizz(3);
y_e=Y(uy)*oneDdphiee(1)+Y(uy+1)*oneDdphiee(2)+Y(uy+2)*oneDdphiee(3);
det=x_z*y_e;
[ U_z,V_e ] = velocity( TrkNode, dphiz,dphie,uvp);


                        if k==1; 

                                                  
   Residu(TrkNode(3,j))=Residu(TrkNode(3,j))+sai(j)*(U_z/x_z+V_e/y_e)*x_z*y_e*w(r)*w(s);  
    
                        end

 jacobi(TrkNode(3,j),TrkNode(1,k))=jacobi(TrkNode(3,j),TrkNode(1,k))...
                              +((sai(j)*(1/x_z)*dphiz(k))*det)*w(r)*w(s); 
      
      
      jacobi(TrkNode(3,j),TrkNode(2,k))=jacobi(TrkNode(3,j),TrkNode(2,k))...
          +((sai(j)*(1/y_e)*dphie(k))*det)*w(r)*w(s); 
      
                   
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
jacobi=jacobi-jacobi';  % by formulation other part of jacobi are earlier transpose 
end

