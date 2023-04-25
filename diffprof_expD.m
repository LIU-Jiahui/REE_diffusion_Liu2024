function [Y_model_epD,R2_epD] = diffprof_expD(HREE_Name,Timemeangeo,D_experi_simu_meangeo,Xi,Yi,Yref)


t = Timemeangeo;
dx = Xi(2)-Xi(1);
dt = dx^2*0.5./max(D_experi_simu_meangeo);
ti = 0:dt:t;
Yref_realnumber = Yref;
Yref_realnumber(isnan(Yref_realnumber))=0;
Yref_realnumber_location = find(Yref_realnumber(:,1));

for i = 1:length(HREE_Name)
    [MtxXi] = MINDIF_Implicit1D(Xi,Yi(:,i),ti,D_experi_simu_meangeo(i));
    Y_model_epD(:,i) = MtxXi(:,end);
    R2_epD(:,i) = 1 - sum((Yref_realnumber(Yref_realnumber_location,i) - MtxXi(Yref_realnumber_location,end)).^2)/sum((Yref_realnumber(Yref_realnumber_location,i) - mean(Yref_realnumber(Yref_realnumber_location,i))).^2);
end

end
