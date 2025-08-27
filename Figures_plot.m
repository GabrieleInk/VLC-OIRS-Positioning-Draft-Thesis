%exportgraphics(gcf,'esempio.pdf','ContentType','vector')
%% Scenario LoS/NLoS
clear variables
close all
load("config.mat");

% Set entities position
LED     = [x_max/2      , y_max/2      , z_max    ];

RIS(1,:)    = [x_max/2      , 0            , z_max/2];
Norm(1,:)=[0 1 0];
alpha(1) = 0;
beta(1) = 0;
w(1)=1;
h(1)=1;

PD_test = [x_max/2, 0, 0];

plotCoverage(LED, RIS, x_max, y_max, z_max, alpha, beta, Norm, w, h);

%% Figure 2

%clear variables
clear; clc; % close all;

figure(2)
load("figure2.mat");
% f1=figure();
figure(2).Position= [0 0 1200 600];
plot(SNR_db,RMSE_RML([1 3 5],:),"LineWidth",2);
hold on
plot(SNR_db,RMSE_ML([1 3 5],:),"LineWidth",2);
hold on
plot(SNR_db,CRLB([1 3 5],:),"LineWidth",2);
hold on
legend("RML - 1 samples","RML - 3 samples","RML - 5 samples",...
    "ML - 1 samples","ML - 3 samples","ML - 5 samples",...
    "DEB - 1 samples","DEB - 3 samples","DEB - 5 samples",'NumColumns',3)
ylabel("RMSE [m]")
xlabel("SNR [dB]")
fontsize(16,"points") %Legend font size
set(gca,"FontSize",22,"FontName","Times New Roman") %Axis font size
grid on
ylim([1e-2 0.5]);
set(gca, 'YScale', 'log')
hold off
title('Fig. 2');

%exportgraphics(gcf,'Fig2.eps','ContentType','vector','BackgroundColor','none')


warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig2_tikz.tex','width','\plotWidth');

%(100./RMSE_RML(1,:)).*RMSE_RML(3,:)
%% Figure 3

%clear variables
clear; clc; % close all;

figure(3)
load("figure3.mat");
% f1=figure(1);

figure(3).Position= [0 0 1200 600];
plot(angles,RMSE_RML([1 3 5],:),"LineWidth",2);
hold on
plot(angles,RMSE_ML([1 3 5],:),"LineWidth",2);
hold on
plot(angles,CRLB([1 3 5],:),"LineWidth",2);
hold on
xlim([67.5 112.5]);
legend("RML - 5 samples","RML - 20 samples","RML - 80 samples",...
    "ML - 5 samples","ML - 20 samples","ML - 80 samples",...
    "DEB - 5 samples","DEB - 20 samples","DEB - 80 samples",'NumColumns',4,'Location','north')

ylabel("RMSE [m]")
xlabel("Angle [°]")
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
grid on
ylim([3e-3 1e-1]);
set(gca, 'YScale', 'log')
hold off
title('Fig. 3');

%exportgraphics(gcf,'Fig3.eps','ContentType','vector','BackgroundColor','none')

% warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig3_tikz.tex','width','\plotWidth');

%% Figure 4

%clear variables
clear; clc; % close all;

figure(4);
load("figure4.mat");

figure(4).Position= [0 0 1200 600];
plot(span,RMSE_RML([1 3 5],:),"LineWidth",2);
hold on
plot(span,RMSE_ML([1 3 5],:),"LineWidth",2);
hold on
plot(span,CRLB([1 3 5],:),"LineWidth",2);
hold on
xlim([1.26 5]);
legend("RML - 5 samples","RML - 20 samples","RML - 80 samples",...
    "ML - 5 samples","ML - 20 samples","ML - 80 samples",...
    "DEB - 5 samples","DEB - 20 samples","DEB - 80 samples",'NumColumns',3,'Location','southeast')

ylabel("RMSE [m]")
xlabel("Y-coordinate [m]")
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
grid on
set(gca, 'YScale', 'log')
hold off
title('Fig. 4');

%exportgraphics(gcf,'Fig4.eps','ContentType','vector','BackgroundColor','none')

warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig4_tikz.tex','width','\plotWidth');
%% Figure 5.1

%clear variables
clear; clc; % close all;

figure(5);
load("figure5_1.mat");

figure(5).Position= [0 0 1200 600];
plot(SNR_db,RMSE_RML(:,:),"LineWidth",2,"LineStyle","--");
hold on
plot(SNR_db,RMSE_ML(:,:),"LineWidth",2,"LineStyle",":");
hold on
plot(SNR_db,CRLB(:,:),"LineWidth",2);
hold on
legend("RML - 1 samples","RML - 3 samples","RML - 5 samples",...
    "ML - 1 samples","ML - 3 samples","ML - 5 samples",...
    "DEB - 1 samples","DEB - 3 samples","DEB - 5 samples",'NumColumns',3)
ylabel("RMSE [m]")
xlabel("SNR [dB]")
fontsize(16,"points")
xlim([30 60]);
set(gca,"FontSize",22,"FontName","Times New Roman")
grid on
set(gca, 'YScale', 'log')
hold off
title('Fig. 5(1)');

%exportgraphics(gcf,'Fig5.eps','ContentType','vector','BackgroundColor','none')

warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig5_tikz.tex','width','\plotWidth');
%% Figure 5.2

%clear variables
clear; clc; % close all;

figure(6);
load("figure5_2.mat");

figure(6).Position= [0 0 1200 600];
plot(SNR_db,RMSE_RML(:,:),"LineWidth",2,"LineStyle","--");
hold on
plot(SNR_db,RMSE_ML(:,:),"LineWidth",2,"LineStyle",":");
hold on
plot(SNR_db,CRLB(:,:),"LineWidth",2);
hold on
legend("RML - 1 samples","RML - 3 samples","RML - 5 samples",...
    "ML - 1 samples","ML - 3 samples","ML - 5 samples",...
    "DEB - 1 samples","DEB - 3 samples","DEB - 5 samples",'NumColumns',3)
ylabel("RMSE [m]")
xlabel("SNR [dB]")
fontsize(16,"points")
xlim([30 60]);
set(gca,"FontSize",22,"FontName","Times New Roman")
grid on
set(gca, 'YScale', 'log')
hold off
title('Fig. 5(2)');

%exportgraphics(gcf,'Fig6.eps','ContentType','vector','BackgroundColor','none');

warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig6_tikz.tex','width','\plotWidth');
%% Figure 6

%clear variables
clear; clc; % close all;

figure(7);
load("figure6.mat");

figure(7).Position= [0 0 1200 600];
plot(1:5,rmse_1(1:2:3,:),"LineWidth",2)
hold on
plot(1:5,crlb_1,"LineWidth",2)
hold on
plot(1:5,rmse_2(1:2:3,:),"LineWidth",2)
hold on
plot(1:5,crlb_2,"LineWidth",2)
legend("RML ILS 1000 Lumen","RML IWLS 1000 Lumen","PEB 1000 Lumen","RML ILS 3000 Lumen","RML IWLS 3000 Lumen","PEB 3000 Lumen",'NumColumns',2)
ylabel("RMSE [m]")
xlabel('Iterations [#]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
grid on
set(gca, 'YScale', 'log')
hold off
title('Fig. 6');

%exportgraphics(gcf,'Fig6.eps','ContentType','vector','BackgroundColor','none');

warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig6_tikz.tex','width','\plotWidth');
%% Figure 7.1
clear variables

figure(8);
load("figure7_1.mat");

figure(8).Position= [0 0 900 700];
s=pcolor(x,y,result(1).crlb');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('Fig. 7(1)');


f2=figure();
f2.Position= [0 0 900 700];
s=pcolor(x,y,result(2).crlb');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off

f3=figure();
f3.Position= [0 0 900 700];
s=pcolor(x,y,result(3).crlb');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off

%% Figure 7.2
clear variables

load("figure7_2_1.mat");

RMSE_1_full = RMSE_1;

n= length(RMSE_1);
%Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_1_full(i, j) )
            % Apply all reflections
            RMSE_1_full(n+1-i, j) = RMSE_1_full(i, j); % Vertical reflection
            RMSE_1_full(i, n+1-j) = RMSE_1_full(i, j); % Horizontal reflection
        end
    end
end

figure(9);
f1.Position= [0 0 900 700];
s=pcolor(x,y,RMSE_1_full');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('Fig. 7(2)');


load("figure7_2_2.mat");

RMSE_2_full = RMSE_2;

n= length(RMSE_2);
% Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_2_full(i, j) )
            % Apply all reflections
            RMSE_2_full(n+1-i, j) = RMSE_2_full(i, j); % Vertical reflection
            RMSE_2_full(i, n+1-j) = RMSE_2_full(i, j); % Horizontal reflection
        end
    end
end

f2=figure();
f2.Position= [0 0 900 700];
s=pcolor(x,y,RMSE_2_full');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off

load("figure7_2_3.mat");

RMSE_3_full = RMSE_3;

n= length(RMSE_3);
% Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_3_full(i, j) )
            % Apply all reflections
            RMSE_3_full(n+1-i, j) = RMSE_3_full(i, j); % Vertical reflection
            RMSE_3_full(i, n+1-j) = RMSE_3_full(i, j); % Horizontal reflection
        end
    end
end

f3=figure();
f3.Position= [0 0 900 700];
s=pcolor(x,y,RMSE_3_full');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('Fig. 7_2');

%% Figure 1-subs

%clear variables
clear; clc; % close all;

figure(10)
f.Position= [0 0 780 1200];
colormap(jet); % Apply jet colormap
tiledlayout(3, 2,'TileSpacing','tight');

caxis_vals = []; % Per rendere le scale coerenti

load("figure7_1.mat");

ax = nexttile;
s=pcolor(x,y,result(1).crlb');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])



set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")



load("figure7_2_1.mat");

RMSE_1_full = RMSE_1;

n= length(RMSE_1);
%Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_1_full(i, j) )
            % Apply all reflections
            RMSE_1_full(n+1-i, j) = RMSE_1_full(i, j); % Vertical reflection
            RMSE_1_full(i, n+1-j) = RMSE_1_full(i, j); % Horizontal reflection
        end
    end
end

ax = nexttile;
s=pcolor(x,y,RMSE_1_full');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_1.mat");

ax = nexttile;
s=pcolor(x,y,result(2).crlb');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_2_2.mat");

RMSE_2_full = RMSE_2;

n= length(RMSE_2);
% Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_2_full(i, j) )
            % Apply all reflections
            RMSE_2_full(n+1-i, j) = RMSE_2_full(i, j); % Vertical reflection
            RMSE_2_full(i, n+1-j) = RMSE_2_full(i, j); % Horizontal reflection
        end
    end
end

ax = nexttile;
s=pcolor(x,y,RMSE_2_full');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_1.mat");

ax = nexttile;
s=pcolor(x,y,result(3).crlb');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_2_3.mat");

RMSE_3_full = RMSE_3;

n= length(RMSE_3);
% Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_3_full(i, j) )
            % Apply all reflections
            RMSE_3_full(n+1-i, j) = RMSE_3_full(i, j); % Vertical reflection
            RMSE_3_full(i, n+1-j) = RMSE_3_full(i, j); % Horizontal reflection
        end
    end
end

ax = nexttile;
s=pcolor(x,y,RMSE_3_full');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")

% Normalizzazione della scala colori
cmin = min(caxis_vals(:, 1));
cmax = max(caxis_vals(:, 2));
for i = 1:6
    nexttile(i);
    caxis([cmin cmax-0.007]);
end

% Creazione di una colorbar unica in alto
%cb = colorbar('northoutside');
cb.Layout.Tile = 'north';

% Sistemazione delle etichette e miglioramento della visualizzazione
ax = findall(gcf, 'Type', 'Axes');
for i = 1:6
ax(i).XTick = []; % Remove x-axis ticks
ax(i).YTick = []; % Remove y-axis ticks
end



exportgraphics(gcf,'Fig8.png','BackgroundColor','none','Resolution', 600);%,'ContentType','vector')

%% Figure 2-subs

%clear variables
clear; clc; % close all;

figure(11);
f.Position= [0 0 780 1200];
colormap(jet); % Apply jet colormap
tiledlayout(3, 2,'TileSpacing','tight');

caxis_vals = []; % Per rendere le scale coerenti

load("figure7_1.mat");

ax = nexttile;
%s=pcolor(x,y,result(1).crlb');
%s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])



%set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")



load("figure7_2_1.mat");

RMSE_1_full = RMSE_1;

n= length(RMSE_1);
%Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_1_full(i, j) )
            % Apply all reflections
            RMSE_1_full(n+1-i, j) = RMSE_1_full(i, j); % Vertical reflection
            RMSE_1_full(i, n+1-j) = RMSE_1_full(i, j); % Horizontal reflection
        end
    end
end

ax = nexttile;
%s=pcolor(x,y,RMSE_1_full');
%s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

%set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_1.mat");

ax = nexttile;
%s=pcolor(x,y,result(2).crlb');
%s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

%set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_2_2.mat");

RMSE_2_full = RMSE_2;
n= length(RMSE_2);
% Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_2_full(i, j) )
            % Apply all reflections
            RMSE_2_full(n+1-i, j) = RMSE_2_full(i, j); % Vertical reflection
            RMSE_2_full(i, n+1-j) = RMSE_2_full(i, j); % Horizontal reflection
        end
    end
end

ax = nexttile;
%s=pcolor(x,y,RMSE_2_full');
%s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

%set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_1.mat");

ax = nexttile;
%s=pcolor(x,y,result(3).crlb');
%s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

%set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_2_3.mat");

RMSE_3_full = RMSE_3;

n= length(RMSE_3);
% Reflect the filled 1/8 part across all symmetry axes
for i = 1:n
    for j = 1:n
        if ~isnan(RMSE_3_full(i, j) )
            % Apply all reflections
            RMSE_3_full(n+1-i, j) = RMSE_3_full(i, j); % Vertical reflection
            RMSE_3_full(i, n+1-j) = RMSE_3_full(i, j); % Horizontal reflection
        end
    end
end

ax = nexttile;
%s=pcolor(x,y,RMSE_3_full');
%s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

%set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")

% Normalizzazione della scala colori
cmin = min(caxis_vals(:, 1));
cmax = max(caxis_vals(:, 2));
for i = 1:6
    nexttile(i);
    caxis([cmin cmax-0.007]);
end

% Creazione di una colorbar unica in alto
cb = colorbar('northoutside');
cb.Layout.Tile = 'north';

% Sistemazione delle etichette e miglioramento della visualizzazione
ax = findall(gcf, 'Type', 'Axes');

%exportgraphics(gcf,'Fig8.png','BackgroundColor','none');%,'ContentType','vector')

warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig8_tikz.tex','width','\plotWidth');

