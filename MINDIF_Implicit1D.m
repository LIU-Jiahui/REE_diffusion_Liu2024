function [MtxXi] = MINDIF_Implicit1D(Xi,Yi,ti,D)
%

dt = ti(2)-ti(1);
dX = Xi(2)-Xi(1);
F = D*dt/dX^2;

MtxXi = zeros(numel(Xi),numel(ti));
MtxXi(:,1) = Yi;

tic
for it = 2:numel(ti)
    
    A = sparse([]);
    A(1,1) = 1;
    A(numel(Xi),numel(Xi)) = 1;
    
    for i = 2:numel(Xi)-1
        A(i,i-1) = -F;
        A(i,i+1) = -F;
        A(i,i) = 1+2*F;
    end
    
    b = MtxXi(:,it-1);
    
    Y = A\b;
    
    MtxXi(:,it) = Y;
end    

t = toc; 
end

%figure
%plot(Xi,Yi,'.-r'), hold on
%plot(Xi,Y,'-k'), hold off
%title(['Implicit (',num2str(t),' sec) | t = ',num2str(ti(it)),' / ',num2str(ti(end))])
%drawnow