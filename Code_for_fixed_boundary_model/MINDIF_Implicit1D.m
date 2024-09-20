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

function [MtxXi] = MINDIF_Implicit1D(Xi,Yi,ti,D)
%

dt = ti(2)-ti(1);
dX = Xi(2)-Xi(1);
F = D*dt/dX^2;

MtxXi = zeros(numel(Xi),numel(ti));
MtxXi(:,1) = Yi;

tic
for it = 2:numel(ti)
    
    A = sparse([]);
    A(1,1) = 1;
    A(numel(Xi),numel(Xi)) = 1;
    
    for i = 2:numel(Xi)-1
        A(i,i-1) = -F;
        A(i,i+1) = -F;
        A(i,i) = 1+2*F;
    end
    
    b = MtxXi(:,it-1);
    
    Y = A\b;
    
    MtxXi(:,it) = Y;
end    

t = toc; 
end

%figure
%plot(Xi,Yi,'.-r'), hold on
%plot(Xi,Y,'-k'), hold off
%title(['Implicit (',num2str(t),' sec) | t = ',num2str(ti(it)),' / ',num2str(ti(end))])
%drawnow