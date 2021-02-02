function fluidflow
% .M file for change as per problem 
% 1.inputdata.m for no of element ,Reynolds no. tolerance
% 2.boundaryCondition.m 
% ********************************************************
%  
% Program illustration for 2X2 element  
%  
%          Quadratic mesh U ,V           Linear mesh for P
%     Qnnod is below                    Lnnod is below
%        21__22__23__24__25             7__________8________9
%        |   |   |    |   |             |          |        |
%        16__17__18__19__20             |    (3)   |  (4)   |
%        |(3)|   |(4) |   |             |          |        |
%        11__12__13__14__15             4__________5________6
%        |   |   |    |   |             |          |        |
%   ^    6___7___8___9___10        ^    |    (1)   |  (2)   |
%   |    |(1)|   |(2) |   |        |    |          |        |
% QnodeY 1___2___3___4____5    LnodeY   1__________2________3
%  QnodeX -->                   LnodeX-->
%  
%  
% for each element=>> let for 4th element 
% in program we work in each loop from element to element to x direction 
% each element we use for quadratic mesh uX,uY initial both are one 
% uX=uX+2 for next element in X direction uY =uY+2 in y direction 
% similarly for linear mesh we us LuX,LuY which increase by unity in x and y
% direction 
%  
% for 4th element uX,uY=3,3      LuX ,LuY=2,2
% for node tracking 
%  
%  uY+2  23____24__25     LuY+1  8________9
%        |     |   |             |        |
%  uY+1  8____ 19__20            |  (4)   |
%        |(4)  |   |             |        |
%    uY  13____14__15       LuY  5________6
%        uX  uX+1  uX+2          LuX     LuX+1
%  
% notification
% Reynolds number         NRe
%  tolerance            tol
%                      quadratic                Linear 
% element in x           nX                    nX
% element in Y           nY                    nY
% mesh point in x,y      X,Y                     
% node in x,y           QnodeX ,QnodeY        LnodeX ,LnodeY
% total node in x,y     tQnodeX ,tQnodeY       tLnodeX ,tLnodeY
% tracked node number     Qnnod                 Lnnod 
%  
%  vorticity               omega
% stream function value     streamline 
% Pressure                     P
% velocity in x,y            U,Y and lower and upper case 
% combined     for Gauss elimination  uvp   which are [U;V;P]%
%

tic,
[ nX,nY,initialX,initialY,NRe,finalX,finalY,tol ] = Inputdata;

[ X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode,Qnnod,Lnnod] =notify( nX,nY,initialX,initialY,finalX,finalY );

[ uvp,omega ,streamline] = BoundaryCondition( X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode);

[ uvp ] =Driver( uvp ,nX,nY,initialX,initialY,NRe,finalX,finalY,tol, X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode,Qnnod,Lnnod  );

[ omega ] =vorticity(  uvp ,nX,nY,initialX,initialY,NRe,finalX,finalY,tol, X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode,Qnnod,Lnnod ,omega );

[ streamline ] =streamfunction(  uvp,nX,nY,initialX,initialY,NRe,finalX,finalY,tol, X,Y,QnodeX,QnodeY,LnodeX,LnodeY, tLnode,tQnode,Qnnod,Lnnod ,omega ,streamline);
[ Ufinal,Vfinal ] = allplots( uvp,QnodeX,QnodeY ,tQnode,tLnode,nX,nY,X,Y,omega,streamline);
toc,
end

