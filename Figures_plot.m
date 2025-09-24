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

%% Figure 2 RMSE confrontato con DEB degli stimatori LoS ML e RML per valori diversi di campioni di luce non coerente K, in funzione del SNR
%clear variables
clear; clc; % close all;

figure(2)
load("figure2.mat");
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

% exportgraphics(gcf,'Fig2.eps','ContentType','vector','BackgroundColor','none')

% warning('Have Matlab2TikZ in your path!');

matlab2tikz('showInfo', false,'Fig2_tikz.tex','width','\plotWidth');

%(100./RMSE_RML(1,:)).*RMSE_RML(3,:)
%% Figure 3 RMSE confrontato con DEB degli stimatori NLoS ML e RML per valori diversi di campioni di luce non coerente K_1 (nel caso NLoS lo chiamiamo K_1), in funzione dell'angolo azimutale, che descrive lo spostamento del PD lungo un arco di circonferenza attorno all'OIRS

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

%% Figure 4 RMSE confrontato con DEB degli stimatori NLoS ML e RML per valori diversi di campioni K_1 (nel caso NLoS lo chiamiamo K_1), in funzione della distanza orizzontale tra PD e OIRS (è la proiezione sul piano orizzontale della distanza tra PD e OIRS)  


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

% exportgraphics(gcf,'Fig4.eps','ContentType','vector','BackgroundColor','none')

% warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig4_tikz.tex','width','\plotWidth');
%% Figure 5.1 RMSE confrontato con DEB degli stimatori NLoS ML e RML per valori diversi di campioni di luce non coerente K_1 (nel caso NLoS lo chiamiamo K_1), in funzione del SNR quando OIRS ha l'orientamento iniziale non ottimizzato

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

% exportgraphics(gcf,'Fig5.eps','ContentType','vector','BackgroundColor','none')

% warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig5_tikz.tex','width','\plotWidth');
%% Figure 5.2 RMSE confrontato con DEB degli stimatori NLoS ML e RML per valori diversi di campioni di luce non coerente K_1 (nel caso NLoS lo chiamiamo K_1), in funzione del SNR quando OIRS ha l'orientamento ottimizzato, perfettamente puntato verso il PD

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

% exportgraphics(gcf,'Fig6.eps','ContentType','vector','BackgroundColor','none');

% warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig6_tikz.tex','width','\plotWidth');
%% Figure 6 RMSE confrontato con PEB degli stimatori IWLS e ILS per due diversi valori di flusso luminoso p (1000 e 3000 lumens), in funzione delle iterazioni dell'algoritmo di adaptive beam steering

clear variables
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

% exportgraphics(gcf,'Fig6.eps','ContentType','vector','BackgroundColor','none');

% warning('Have Matlab2TikZ in your path!');

% matlab2tikz('Fig6_tikz.tex','width','\plotWidth');
%% Figure 7.1 Performance di localizzazione basata sulla legenda di colori (freddo (~1mm) - caldo (~4mm) quindi errore commettibile più alto) con PEB, flusso luminoso fisso a 1000lm  e 2 iterazioni dell'algoritmo di adaptive beam steering

clear variables

load("figure7_1.mat");
f8=figure();
f8.Position= [0 0 900 700];
s=pcolor(x,y,result(1).crlb');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('4 OIRS CRLB (PEB)');

f9=figure();
f9.Position= [0 0 900 700];
s=pcolor(x,y,result(2).crlb');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('8 OIRS CRLB (PEB)');

f10=figure();
f10.Position= [0 0 900 700];
s=pcolor(x,y,result(3).crlb');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('12 OIRS CRLB (PEB)');

%% Figure 7.2 Performance di localizzazione basata sulla legenda di colori (freddo (~1mm) - caldo (~4mm) quindi errore commettibile più alto) con RMSE, flusso luminoso fisso a 1000lm  e 2 iterazioni dell'algoritmo di adaptive beam steering

%clear variables

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

f11=figure();
f11.Position= [0 0 900 700];
s=pcolor(x,y,RMSE_1_full');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('4 OIRS RMSE');


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

f12=figure();
f12.Position= [0 0 900 700];
s=pcolor(x,y,RMSE_2_full');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('8 OIRS RMSE');

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

f13=figure();
f13.Position= [0 0 900 700];
s=pcolor(x,y,RMSE_3_full');
s.FaceColor = 'interp';
set(s,'LineStyle','none')
colorbar
ylabel("y [m]")
xlabel('x [m]');
fontsize(16,"points")
set(gca,"FontSize",22,"FontName","Times New Roman")
hold off
title('12 OIRS RMSE');

%% Figure 1-subs 

%clear variables
clear; clc; % close all;

f14=figure();
f14.Position= [0 0 780 1200];
colormap(jet); % Apply jet colormap
tiledlayout(3, 2,'TileSpacing','tight');

caxis_vals = []; % Per rendere le scale coerenti

load("figure7_1.mat");


% Primo subplot, 4 OIRS CRLB
ax = nexttile;
s=pcolor(x,y,result(1).crlb');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")
title('4 OIRS CRLB (PEB)');

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


% Secondo subplot, 4 OIRS RMSE
ax = nexttile;
s=pcolor(x,y,RMSE_1_full');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")
title('4 OIRS RMSE');

load("figure7_1.mat");


% Terzo subplot, 8 OIRS CRLB
ax = nexttile;
s=pcolor(x,y,result(2).crlb');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")
title('8 OIRS CRLB (PEB)');

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


% Quarto subplot, 8 OIRS RMSE
ax = nexttile;
s=pcolor(x,y,RMSE_2_full');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")
title('8 OIRS RMSE');

load("figure7_1.mat");


% Quinto subplot, 12 OIRS CRLB
ax = nexttile;
s=pcolor(x,y,result(3).crlb');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")
title('12 OIRS CRLB (PEB)');

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


% Sesto subplot, 12 OIRS RMSE
ax = nexttile;
s=pcolor(x,y,RMSE_3_full');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")
title('12 OIRS RMSE');

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
for i = 1:6
ax(i).XTick = []; % Remove x-axis ticks
ax(i).YTick = []; % Remove y-axis ticks
end

exportgraphics(gcf,'Fig8.png','BackgroundColor','none','Resolution', 600);%,'ContentType','vector')

%% Figure 2-subs 

% clear variables
clear; clc; % close all;

f15=figure();
f15.Position= [0 0 780 1200];
colormap(jet); % Apply jet colormap
tiledlayout(3, 2,'TileSpacing','tight');

caxis_vals = []; % Per rendere le scale coerenti

load("figure7_1.mat");


% Primo subplot, 4 OIRS CRLB
ax = nexttile;
s=pcolor(x,y,result(1).crlb');
s.FaceColor = 'interp';
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

% Secondo subplot, 4 OIRS RMSE
ax = nexttile;
s=pcolor(x,y,RMSE_1_full');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

%set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_1.mat");


% Terzo subplot, 8 OIRS PEB
ax = nexttile;
s=pcolor(x,y,result(2).crlb');
s.FaceColor = 'interp';
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


% Quarto subplot, 8 OIRS RMSE
ax = nexttile;
s=pcolor(x,y,RMSE_2_full');
s.FaceColor = 'interp';
caxis_vals = [caxis_vals; caxis];
pbaspect([1 1 1])

%set(s,'LineStyle','none')
ylabel("y [m]")
xlabel('x [m]');set(gca,"FontSize",14,"FontName","Times New Roman")


load("figure7_1.mat");


% Quinto subplot, 12 OIRS PEB
ax = nexttile;
s=pcolor(x,y,result(3).crlb');
s.FaceColor = 'interp';
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


% Sesto subplot, 12 OIRS RMSE
ax = nexttile;
s=pcolor(x,y,RMSE_3_full');
s.FaceColor = 'interp';
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

% exportgraphics(gcf,'Fig8.png','BackgroundColor','none');%,'ContentType','vector')

% warning('Have Matlab2TikZ in your path!');

matlab2tikz('Fig8_tikz.tex','width','\plotWidth');

