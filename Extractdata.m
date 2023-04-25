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