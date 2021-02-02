function [ sai  ] = bilinearBF( zai ,eta)

sai1= (1-zai)*(1-eta);    
sai2= zai*(1-eta);
sai3= (1-zai)*eta;
sai4= zai*eta;
sai=[sai1,sai2,sai3,sai4];
end

