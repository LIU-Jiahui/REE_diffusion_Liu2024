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

function [Xi,Yref] = Extractdata(Filename,HREE_Name,Pointstart,Pointend,Pointeliminate)
% Extract data


% Input data
File_Data = readtable(Filename);
File_Name = File_Data.Properties.VariableNames;
File_Comp = table2array(File_Data);

% Extract data
[IsM_Name,Ind_Name] = ismember(HREE_Name,File_Name);
REE_Comp = File_Comp(:,Ind_Name);
REE_Comp(Pointeliminate,:) = NaN;
Yref = REE_Comp(Pointstart:Pointend,:); 
Xi = File_Comp(Pointstart:Pointend,1);
end