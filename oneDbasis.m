function [ oneDphiz,oneDdphizz, oneDphie,oneDdphiee] = oneDbasis( zai,eta )
oneDphiz(1)=(1-3*zai+2*zai^2);
oneDphiz(2)=(4*(zai-zai^2));
oneDphiz(3)=(-zai+2*zai^2);

oneDdphizz(1)=(-3+4*zai);
oneDdphizz(2)=(4*(1-2*zai));
oneDdphizz(3)=(-1+4*zai);

oneDphie(1)=(1-3*eta+2*eta^2);
oneDphie(2)=(4*(eta-eta^2));
oneDphie(3)=(-eta+2*eta^2);

oneDdphiee(1)=(-3+4*eta);
oneDdphiee(2)=(4*(1-2*eta));
oneDdphiee(3)=(-1+4*eta);
end

