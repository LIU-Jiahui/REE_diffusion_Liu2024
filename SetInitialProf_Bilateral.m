function [Yi] = SetInitialProf_Bilateral(Yref,Count_max_left,Count_max_right,Count_mins,HREE_Name)
% Set initial compositional profile for bilateral transect

Y1 = max(Yref(1:Count_max_left,:));
Yend = max(Yref(end-Count_max_right:end,:));
Ysort = sort(Yref);
Ymid = mean(Ysort(1:Count_mins,:));
Y1_2_ones = ones(2,length(Y1));
Yend_2_end_ones = ones(2,length(Y1));
Y3_end3_ones = ones(length(Yref)-4,length(Y1));
Yi = zeros(length(Yref),length(HREE_Name));
Yi(1:2,:) = Y1.*Y1_2_ones;
Yi(end-1:end,:) = Yend.*Yend_2_end_ones;
Yi(3:end-2,:) = Ymid.*Y3_end3_ones;

end