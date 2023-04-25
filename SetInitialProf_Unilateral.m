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