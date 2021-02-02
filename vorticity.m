function [ omega ] =vorticity(  uvp ,nX,nY,initialX,initialY,NRe,finalX,finalY,tol, X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode,Qnnod,Lnnod,omega  )

% [ w ,gp ] = gausspoints;
[pt,ptt ] = points;
ptx=[pt,pt,pt];
pty=[pt',pt',pt'];
pty=pty';
pty=pty(:);
pty=pty';
uy=1;
Luy=1;
for iY=1:nY
    ux=1;
    Lux=1;
    for iX=1:nX
        [ TrkNode ] = nodetrack( ux,uy,nX,nY,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod,Lux,Luy);
    
 for  j=1:9



 zai=ptx(j);
 eta=pty(j);


 [ phi ,dphie ,dphiz ] = biquadBF( zai ,eta);
 [ oneDphiz,oneDdphizz, oneDphie,oneDdphiee] = oneDbasis( zai,eta );

  [ sai ] = bilinearBF( zai ,eta);
  x_z=X(ux)*oneDdphizz(1)+X(ux+1)*oneDdphizz(2)+X(ux+2)*oneDdphizz(3);
y_e=Y(uy)*oneDdphiee(1)+Y(uy+1)*oneDdphiee(2)+Y(uy+2)*oneDdphiee(3);

[  U_z,V_e,U_e,Ui,V_z,Vi,P] = veloVec( TrkNode, dphiz,dphie,phi,sai,uvp );


                       
                                                  
   
     omega(TrkNode(1,j))=(-U_e/y_e+V_z/x_z);
                    


  
  
  end
ux=ux+2;
Lux=Lux+1;
    end
uy=uy+2;  
Luy=Luy+1;

  
end



end

