function [ jacobi ] = jacobiwrtbdry( Qnnod, jacobi,tQnode,QnodeX,QnodeY,nX,nY )



  temp1(1:QnodeX,1)=Qnnod(1,1:QnodeX);
  temp1(QnodeX+1:2*QnodeX,1)=Qnnod(QnodeY,1:QnodeX); 
  temp1(2*QnodeX+1:2*QnodeX+QnodeY-2,1)=Qnnod(2:QnodeY-1,1); 
  temp1(2*QnodeX+QnodeY-2+1:2*QnodeX+2*(QnodeY-2),1)=Qnnod(2:QnodeY-1,QnodeX);
  temp2=tQnode*ones((2*QnodeX+2*QnodeY-4),1);
  temp2=temp2+temp1;
  temp=[temp1;temp2];
  for ii=1:2*(2*QnodeX+2*QnodeY-4)
  jacobi(temp(ii),1:2*tQnode+(nX+1)*(nY+1))=0;
  jacobi(temp(ii),temp(ii))=1;
  end
  

end
