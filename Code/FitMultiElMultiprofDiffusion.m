% REE_diffusion_Liu2023 is a free code for rare earth element diffusion
% modeling, which developed for Liu et al. 2023 manuscript "Ultra-fast mineral growth
% during regional metamorphism"
%
% Copyright © 2022-2023 University of Bern, Institute of Geological
% Sciences, Pierre Lanari, Jiahui Liu
%
% REE_diffusion_Liu2023 is free code: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any 
% later version.
%
% REE_diffusion_Liu2023 is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with REE_diffusion_Liu2023. If not, see https://www.gnu.org/licenses.


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