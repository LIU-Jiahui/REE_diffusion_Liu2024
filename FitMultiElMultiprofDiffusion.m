function [Time,Y] = FitMultiElMultiprofDiffusion(D,Xi,Yi,Yref)
% Diffusion modeling of multi-elements & multi-profiles together

% Calculate diffusion time for each simulated D group
fun = @Multiprof_ObjectiveFct_opt_t;
Options = optimset('Display','iter','TolX',1e-6,'MaxFunEvals',500);
X0 = 500*20^2*0.5/max(D);
X = fminsearch(fun,X0,Options,Xi,Yi,Yref,D);
Time = X;

% Calculate effective D for each diffusion time
tfixed = X;
fun = @Multiprof_ObjectiveFct_opt_D;
Options = optimset('Display','iter','TolX',1e-6,'MaxFunEvals',500);
Y0 = D;

% Y = cell(1,length(Xi));
for i = 1:length(Xi)
    for j = 1:length(D)
        Y(j,i) = fminsearch(fun,Y0(j),Options,Xi{i},Yi{i}(:,j),Yref{i}(:,j),tfixed);
    end
end

end