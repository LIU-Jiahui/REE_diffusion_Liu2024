function [Residual] = Multiprof_ObjectiveFct_opt_t(t,Xi,Yi,Yref,D)

dX = Xi{1}(2)-Xi{1}(1);
dt = dX^2*0.5/max(D);
ti = 0:dt:t;

for i = 1:length(Xi)
    for j = 1:size(Yi{1},2)
        [MtxXi{i}] = MINDIF_Implicit1D(Xi{i},Yi{i}(:,j),ti,D(j));
        delta{i} = MtxXi{i}(:,end)-Yref{i}(:,j);
        delta{i}(isnan(delta{i}))=0;
        Residual_single{i}(:,j) = sqrt(sum(delta{i}.^2));
    end
    Residual_p{i} = sum(Residual_single{i});
end

Residual = 0;
for i = 1:length(Residual_p)
    Residual = Residual + Residual_p{i};
end

end