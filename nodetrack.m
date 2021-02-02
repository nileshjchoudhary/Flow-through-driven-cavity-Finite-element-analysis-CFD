function [ TrkNode ] = nodetrack( ux,uy,nX,nY,QnodeX,QnodeY,LnodeX,LnodeY,tLnode,tQnode,Qnnod,Lnnod,Lux,Luy)


temp2=tQnode*ones(1,9);
temp1=zeros(1,9);
skipy=0;
skip=1;
for ty=1:3;
    skipx=0;
for tx= 1:3;
        temp1(1,skip)=Qnnod((uy+skipy),(ux+skipx));
skip=skip+1;
skipx=skipx+1;
end
skipy=skipy+1;
end
temp2=temp2+temp1;
temp3=zeros(1,9);
temp3(1,1)=Lnnod(Luy,Lux)+2*tQnode;
temp3(1,2)=Lnnod(Luy,Lux+1)+2*tQnode;
temp3(1,3)=Lnnod(Luy+1,Lux)+2*tQnode;
temp3(1,4)=Lnnod(Luy+1,Lux+1)+2*tQnode;

TrkNode=[temp1;temp2;temp3];

end

