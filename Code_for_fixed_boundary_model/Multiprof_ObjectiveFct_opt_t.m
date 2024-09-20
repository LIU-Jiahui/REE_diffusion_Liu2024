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