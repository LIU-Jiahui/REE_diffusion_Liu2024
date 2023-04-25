function [Residual] = Multiprof_ObjectiveFct_opt_D(D,Xi,Yi,Yref,t)

dX = Xi(2)-Xi(1);
dt = t/1000;
% dt = dX^2*0.5/D;
ti = 0:dt:t;

for i = 1:size(Yi,2)
    [MtxXi] = MINDIF_Implicit1D(Xi,Yi(:,i),ti,D);
    delta = MtxXi(:,end)-Yref(:,i);
    delta(isnan(delta))=0;
    Residual_single(:,i) = sqrt(sum(delta.^2));
end
Residual = sum(Residual_single);

end
