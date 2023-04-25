
%% Diffusion modeling for HREE+Y in melt (inherited by garnet)
clear all
clc
close all

%% Set HREE elements and diffusion conditions
HREE_Name = {'Lu','Yb','Er','Ho','Dy','Y'};

T_C = 800;                                                % Temperature in °C
T_K = T_C + 273.15;                                       % Temperature in Kelvin
T_K_sigma = 25;                                           % Uncertainties of temperature
logD0 = [-3.74,-5.23,-3.39,-4.65,-4.60,-5.09];            % D0 at 800C (J/K), Table 3, Holycross & Watson, 2018
logD0_sigma = [1.40,0.84,1.49,1.07,1.04,1.01];            % Uncertainties of D0
Ea = [198.97,164.47,205.25,175.99,176.47,165.74];         % Ea at 800C (KJ), Table 3, Holycross & Watson, 2018
Ea_sigma = [31.24,18.78,33.36,23.78,23.27,22.51];         % Uncertainties of Ea

% Simulated experimental D using Monte Carlo
N = 20;                                                 % No. of Tries 
T_K_simu = T_K + T_K_sigma.*randn(N,1);                   % Simulated temperature
logD0_simu = logD0 + logD0_sigma.*randn(N,1);             % Simulated logD0
Ea_simu = Ea + Ea_sigma.*randn(N,1);                      % Simulated Ea
D_experi_simu = (10.^(logD0_simu).*exp(-Ea_simu.*1000./8.314./T_K_simu)).*10.^12./(1./(3600*24)); % Simulated experimental D using Monte Carlo, uncertainties from Temerature, logD0 and Ea.

% make sure each experimental D < 1.9e7, otherwise it exceeds the statistical range
for i =1:size(D_experi_simu)
    if max(D_experi_simu(i,:)) > 1.9e7
        while 1
            T_K_simu2 = T_K + T_K_sigma.*randn(1,1);
            logD0_simu2 = logD0 + logD0_sigma.*randn(1,1);
            Ea_simu2 = Ea + Ea_sigma.*randn(1,1);
            D_simu2 = (10.^(logD0_simu2).*exp(-Ea_simu2.*1000./8.314./T_K_simu2)).*10.^12./(1./(3600*24));
            if max(D_simu2) < 1.9e7
                D_experi_simu(i,:) = D_simu2;
                break
            end
        end
    end
end
D_experi_simu_meangeo = exp(mean(log(D_experi_simu)));    % μ*, median (geometric mean) of simulated experimental D
D_experi_simu_sigmageo = exp(std(log(D_experi_simu)));    % σ*, multiplicative standard deviation (geometric "deviation")
D_experi_simu_upper = D_experi_simu_meangeo.*D_experi_simu_sigmageo;     % μ*×σ*
D_experi_simu_lower = D_experi_simu_meangeo./D_experi_simu_sigmageo;     % μ*/σ*
uppererror = D_experi_simu_upper - D_experi_simu_meangeo;
lowererror = D_experi_simu_meangeo - D_experi_simu_lower;

%% Extract data & Set initial profiles
Profile_Name = {'Profile1','Profile2','Profile3','Profile4','Profile5','Profile6'};
[Filename,Pointstart,Pointend,Pointeliminate,Count_max_left,Count_max_right,Count_mins,Xi,Yi,Yref] = deal(cell(1,length(Profile_Name)));

% Profile 1
Filename{1} = 'Data1.txt';
Pointstart{1} = 13;
Pointend{1} = 76;
Pointeliminate{1} = [34,35,61,62,44,45,46,58,59,60,65,66];
[Xi{1},Yref{1}] = Extractdata(Filename{1},HREE_Name,Pointstart{1},Pointend{1},Pointeliminate{1});
Count_mins{1} = 9;
[Yi{1}] = SetInitialProf_Unilateral(Yref{1},Count_mins{1},HREE_Name);

% Profile 2
Filename{2} = 'Data2.txt';
Pointstart{2} = 13;
Pointend{2} = 67;
Pointeliminate{2} = [28,32,41,45,33,40,56,57,58,59,60];
[Xi{2},Yref{2}] = Extractdata(Filename{2},HREE_Name,Pointstart{2},Pointend{2},Pointeliminate{2});
Count_mins{2} = 6;
[Yi{2}] = SetInitialProf_Unilateral(Yref{2},Count_mins{2},HREE_Name);

% Profile 3
Filename{3} = 'Data3.txt';
Pointstart{3} = 4;
Pointend{3} = 87;
Pointeliminate{3} = [59,60,71];
[Xi{3},Yref{3}] = Extractdata(Filename{3},HREE_Name,Pointstart{3},Pointend{3},Pointeliminate{3});
Count_mins{3} = 7;
[Yi{3}] = SetInitialProf_Unilateral(Yref{3},Count_mins{3},HREE_Name);

% Profile 4
Filename{4} = 'Data4.txt';
Pointstart{4} = 3;
Pointend{4} = 47;
Pointeliminate{4} = [26,27,24,25,38,39];
[Xi{4},Yref{4}] = Extractdata(Filename{4},HREE_Name,Pointstart{4},Pointend{4},Pointeliminate{4});
Count_mins{4} = 7;
[Yi{4}] = SetInitialProf_Unilateral(Yref{4},Count_mins{4},HREE_Name);

% Profile 5
Filename{5} = 'Data5.txt';
Pointstart{5} = 3;
Pointend{5} = 49;
Pointeliminate{5} = [12,34,11,13,14,16,17,28,29,30,31,32,33];
[Xi{5},Yref{5}] = Extractdata(Filename{5},HREE_Name,Pointstart{5},Pointend{5},Pointeliminate{5});
Count_max_left{5} = 5;
Count_max_right{5} = 5;
Count_mins{5} = 5;
[Yi{5}] = SetInitialProf_Bilateral(Yref{5},Count_max_left{5},Count_max_right{5},Count_mins{5},HREE_Name);

% Profile 6
Filename{6} = 'Data6.txt';
Pointstart{6} = 8;
Pointend{6} = 74;
Pointeliminate{6} = [30,31,48,49,65];
[Xi{6},Yref{6}] = Extractdata(Filename{6},HREE_Name,Pointstart{6},Pointend{6},Pointeliminate{6});
Count_max_left{6} = 5;
Count_max_right{6} = 5;
Count_mins{6} = 5;
[Yi{6}] = SetInitialProf_Bilateral(Yref{6},Count_max_left{6},Count_max_right{6},Count_mins{6},HREE_Name);

%% Estimate diffusion time, effective diffusion coeffiecient, and geometry parameter
[D_eff,D_eff_meangeo,G] = deal(cell(1,length(Profile_Name))); % effective D, mean effective D and geometry parameter in the porosity medium
Y = zeros(N,length(HREE_Name),length(Xi));

for i = 1:N
    [Time(i),Y(i,:,:)] = FitMultiElMultiprofDiffusion(D_experi_simu(i,:),Xi,Yi,Yref);
end
for i = 1:length(Xi)
    D_eff{i} = Y(:,:,i);
end

% Estimate diffusion time with uncertainties. For log-normal distribution, the interval (68% probability) is [μ*/σ*,μ*×σ*]
Timemeangeo = exp(mean(log(Time)));               % μ*(geometric "mean") = exp(μ)
Timesigmageo = exp(std(log(Time)));               % σ*(multiplicative standard deviation, geometric "deviation") = exp(σ)
Timeupper = Timemeangeo*Timesigmageo;             % μ*×σ*
Timelower = Timemeangeo/Timesigmageo;             % μ*/σ*
Timeconfiinterval = find(Time >= Timelower & Time <= Timeupper);      % confidence interval
disp(['Method 1  ',num2str(Timemeangeo),'   + ',num2str(Timeupper-Timemeangeo), '   -',num2str(Timemeangeo-Timelower),'   n=',num2str(length(Timeconfiinterval))]);

for i = 1:length(D_eff_meangeo)
    D_eff_meangeo{i} = exp(mean(log(D_eff{i})));  % median (geometric mean) of effective D_simu (effective)
end

% Plot experimental D and effective D (meangeo)
xerrorbar = 1:1:5;
errorbar(xerrorbar,D_experi_simu_meangeo(xerrorbar),lowererror(xerrorbar),uppererror(xerrorbar));
hold on
Dff_mean = {'error','Deff_mean_tr1','Deff_mean_tr2','Deff_mean_tr3','Deff_mean_tr4','Deff_mean_tr5','Deff_mean_tr6'}
for i = 1:length(D_eff_meangeo)
    plot(xerrorbar,D_eff_meangeo{i}(1:length(xerrorbar)),'o')
    hold on
end
legend(Dff_mean)
ax = gca;
ax.YScale = 'log';
axis([0,7,1e1,1e6])

% Estimate geometry parameter
porosity = 0.1;                                   % ε (porosity), assumption as melt content = 10% in phase equilibrium modeling
for i = 1:length(G)
    G{i} = porosity./(D_eff{i}./D_experi_simu);   % G means geometry parameter, the effective diffusion coefficient in porous media. Deff/Dexp = ε/G
end

G_HREE_Name = {'Lu','Yb','Er','Ho','Dy'};
X_plot = 1:1:length(G_HREE_Name);                 % plot G for each profile
figure
for i = 1:length(G)
    plot(X_plot,G{i}(:,1:5),'.','MarkerSize',8)
    plot(X_plot,mean(G{i}(:,1:5)),'o-','LineWidth',1.5)
    hold on
end

title('Geometric parameter')
xlim([0 6])
ylim([0 2.5])
xticklabels(G_HREE_Name)
legend(Profile_Name)
newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#F90'};
colororder(newcolors)
set(gca,'FontSize',40);
set(gcf,'PaperSize',[25 25])

%% Plot diffusion profile using fixed time & D
[Y_model_expD,Y_model_effD] = deal(cell(1,length(Profile_Name)));
for i = 1:length(Xi)
    [Y_model_expD{i},R2_epD{i}] = diffprof_expD(HREE_Name,Timemeangeo,D_experi_simu_meangeo,Xi{i},Yi{i},Yref{i}); 
    [Y_model_effD{i},R2_eff{i}] = diffprof_effD(HREE_Name,Timemeangeo,D_eff_meangeo{i},Xi{i},Yi{i},Yref{i});
end

Xi_nor = Xi;   % Reset the initial position of X (distance) as start from 0 micrometer.
for i = 1:length(Xi)
    Xi_nor{i} = Xi{i}-Xi{i}(1);
end

% Plot diffusion profile using experimental D
for i = 1:length(HREE_Name)
    figure
    count = 0;
    for j = 1:length(Xi)
        count = count + 1;
        plot3(Xi_nor{j},count.*(ones(length(Xi_nor{j}),1)),Yref{j}(:,i),'k.','MarkerSize',10)
        hold on
        plot3(Xi_nor{j},count.*(ones(length(Xi_nor{j}),1)),Y_model_expD{j}(:,i),'b-','Linewidth',0.5)
        hold on
        plot3(Xi_nor{j},count.*(ones(length(Xi_nor{j}),1)),Yi{j}(:,i),'r-','Linewidth',0.5)
        hold on
    end
    title('Diffusion modeling profile (using experimental D)', HREE_Name(i))
    xlabel('Distance (μm)')
    ylabel('Profile No.')
    zlabel('Mass fraction (μg/g)')
    yticks([0:1:6])
    set(gca,'FontSize',40);
end

% Plot diffusion profile using effective D
for i = 1:length(HREE_Name)
    figure
    count = 0;
    for j = 1:length(Xi)
        count = count + 1;
        plot3(Xi_nor{j},count.*(ones(length(Xi_nor{j}),1)),Yref{j}(:,i),'k.','MarkerSize',10)
        hold on
        plot3(Xi_nor{j},count.*(ones(length(Xi_nor{j}),1)),Y_model_effD{j}(:,i),'b-','Linewidth',0.5)
        hold on
        plot3(Xi_nor{j},count.*(ones(length(Xi_nor{j}),1)),Yi{j}(:,i),'r-','Linewidth',0.5)
        hold on
    end
    title('Diffusion modeling profile (using effective D)', HREE_Name(i))
    xlabel('Distance (μm)')
    ylabel('Profile No.')
    zlabel('Mass fration (μg/g)')
    yticks([0:1:6])
    set(gca,'FontSize',40);
end
