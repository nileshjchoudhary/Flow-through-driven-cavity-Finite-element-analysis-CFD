function [ U_z,V_e ] = velocity( TrkNode, dphiz,dphie,uvp)

U_z=0;
V_e=0;
for vi=1:9;
    U_z=U_z+dphiz(1,vi)*uvp(TrkNode(1,vi));
    V_e=V_e+dphie(1,vi)*uvp(TrkNode(2,vi));
end

end

