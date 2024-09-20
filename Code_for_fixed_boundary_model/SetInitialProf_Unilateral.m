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


function [Yi] = SetInitialProf_Unilateral(Yref,Count_mins,HREE_Name)
% Set initial compositional profile for unilateral transect

Y1 = max(Yref);
Y2 = mean(Yref(end-Count_mins:end,:));
Y1_2_ones = ones(2,length(Y1));
Y3_end_ones = ones(length(Yref)-2,length(Y1));
Yi = zeros(length(Yref),length(HREE_Name));
Yi(1:2,:) = Y1.*Y1_2_ones;
Yi(3:end,:) = Y2.*Y3_end_ones;

end