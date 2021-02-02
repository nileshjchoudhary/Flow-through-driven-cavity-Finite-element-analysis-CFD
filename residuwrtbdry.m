function [ Residu] = residuwrtbdry( Residu,Qnnod ,QnodeX,QnodeY,nX,nY ,tQnode )

temp4(1:QnodeX,1)=Qnnod(1,1:QnodeX);
temp4(QnodeX+1:2*QnodeX,1)=Qnnod(QnodeY,1:QnodeX);
temp4(2*QnodeX+1:2*QnodeX+QnodeY-2,1)=Qnnod(2:QnodeY-1,1);
temp4(2*QnodeX+QnodeY-1:2*QnodeX+2*(QnodeY-2),1)=Qnnod(2:QnodeY-1,QnodeX);
temp5=temp4+tQnode*ones(2*QnodeX+2*(QnodeY-2),1);
track=[temp4;temp5;];
for tr=1:2*(2*QnodeX+2*(QnodeY-2))
    Residu(track(tr))=0;
end
Residu=-Residu;


end

