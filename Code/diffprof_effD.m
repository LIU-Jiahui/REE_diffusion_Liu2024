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

function [Y_model_effD,R2_effD] = diffprof_effD(HREE_Name,Timemeangeo,D_eff_meangeo,Xi,Yi,Yref)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

t = Timemeangeo;
dx = Xi(2)-Xi(1);
dt = dx^2*0.5./max(D_eff_meangeo);
ti = 0:dt:t;
Yref_realnumber = Yref;
Yref_realnumber(isnan(Yref_realnumber))=0;
Yref_realnumber_location = find(Yref_realnumber(:,1));

for i = 1:length(HREE_Name)
    [MtxXi] = MINDIF_Implicit1D(Xi,Yi(:,i),ti,D_eff_meangeo(i));
    Y_model_effD(:,i) = MtxXi(:,end);
    R2_effD(:,i) = 1 - sum((Yref_realnumber(Yref_realnumber_location,i) - MtxXi(Yref_realnumber_location,end)).^2)/sum((Yref_realnumber(Yref_realnumber_location,i) - mean(Yref_realnumber(Yref_realnumber_location,i))).^2);
end

end
