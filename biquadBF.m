function [ phi ,dphie ,dphiz ] = biquadBF( zai ,eta)


%*************** Bi_Quadratic Basis Functions************%
phi1=(1-3*zai+2*zai^2)*(1-3*eta+2*eta^2);
phi4=(1-3*zai+2*zai^2)*(4*(eta-eta^2));
phi7=(1-3*zai+2*zai^2)*(-eta+2*eta^2);
phi2=(4*(zai-zai^2))*(1-3*eta+2*eta^2);
phi5=(4*(zai-zai^2))*(4*(eta-eta^2));
phi8=(4*(zai-zai^2))*(-eta+2*eta^2);
phi3=(-zai+2*zai^2)*(1-3*eta+2*eta^2);
phi6=(-zai+2*zai^2)*(4*(eta-eta^2));
phi9=(-zai+2*zai^2)*(-eta+2*eta^2);

%**********  Derivative with respect to zai  ***************%

dphi1z=(-3+4*zai)*(1-3*eta+2*eta^2);
dphi4z=(-3+4*zai)*(4*(eta-eta^2));
dphi7z=(-3+4*zai)*(-eta+2*eta^2);
dphi2z=4*(1-2*zai)*(1-3*eta+2*eta^2);
dphi5z=4*(1-2*zai)*(4*(eta-eta^2));
dphi8z=4*(1-2*zai)*(-eta+2*eta^2);
dphi3z=(-1+4*zai)*(1-3*eta+2*eta^2);
dphi6z=(-1+4*zai)*(4*(eta-eta^2));
dphi9z=(-1+4*zai)*(-eta+2*eta^2);

%**********  Derivative with respect to eta  ***************%
dphi1e=(1-3*zai+2*zai^2)*(-3+4*eta);
dphi4e=(1-3*zai+2*zai^2)*4*(1-2*eta);
dphi7e=(1-3*zai+2*zai^2)*(-1+4*eta);
dphi2e=(4*(zai-zai^2))*(-3+4*eta);
dphi5e=(4*(zai-zai^2))*4*(1-2*eta);
dphi8e=(4*(zai-zai^2))*(-1+4*eta);
dphi3e=(-zai+2*zai^2)*(-3+4*eta);
dphi6e=(-zai+2*zai^2)*4*(1-2*eta);
dphi9e=(-zai+2*zai^2)*(-1+4*eta);
phi=[phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9];

dphie=[dphi1e,dphi2e,dphi3e,dphi4e,dphi5e,dphi6e,dphi7e,dphi8e,dphi9e];

dphiz=[dphi1z,dphi2z,dphi3z,dphi4z,dphi5z,dphi6z,dphi7z,dphi8z,dphi9z];

end

