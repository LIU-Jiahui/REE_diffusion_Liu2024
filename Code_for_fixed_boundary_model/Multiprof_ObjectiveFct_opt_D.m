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
