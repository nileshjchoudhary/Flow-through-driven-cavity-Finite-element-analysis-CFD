function [  U_z,V_e,U_e,Ui,V_z,Vi,P] = veloVec( TrkNode, dphiz,dphie,phi,sai,uvp )



U_z=0;
V_e=0;
U_e=0;
Ui=0;
V_z=0;
Vi=0;
P=0;

for vi=1:9
    U_z=U_z+dphiz(1,vi)*uvp(TrkNode(1,vi));
    Ui=Ui+phi(1,vi)*uvp(TrkNode(1,vi));
    U_e=U_e+dphie(1,vi)*uvp(TrkNode(1,vi));
   
    Vi=Vi+phi(1,vi)*uvp(TrkNode(2,vi));
    V_z=V_z+dphiz(1,vi)*uvp(TrkNode(2,vi));
    V_e=V_e+dphie(1,vi)*uvp(TrkNode(2,vi));
end

for pi=1:4
    P=P+sai(1,pi)*uvp(TrkNode(3,pi));
end


end

