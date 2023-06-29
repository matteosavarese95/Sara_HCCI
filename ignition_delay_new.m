close all;
clear all;

% Questo serve perche entriamo nella cartella del batch e senno non trova
% le funzioni
working_dir = pwd;
addpath(working_dir);

% Ci serve solo il vettore degli indici dove andare a prendere la ropa
load ignition_delay.mat ind

set(0, 'defaulttextfontsize', 16);
set(0, 'defaultaxesfontsize', 16);

% Metti il path al tuo OpenSMOKE
addpath /Users/matteosavarese/Desktop/Dottorato/OpenSmoke/opensmoke++suite-0.15.0.serial.osx/lib/;
addpath /Users/matteosavarese/Desktop/Dottorato/OpenSmoke/opensmoke++suite-0.15.0.serial.osx/bin/;    

% Anche qua
command_run = 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/matteosavarese/Desktop/Dottorato/OpenSmoke/opensmoke++suite-0.15.0.serial.osx/lib ; /Users/matteosavarese/Desktop/Dottorato/OpenSmoke/opensmoke++suite-0.15.0.serial.osx/bin/OpenSMOKEpp_BatchReactor.sh --input input.dic';

%% Parametric study definition

% Select the parameters we want to study
T_span = [1190 1215 1240 1265 1290];
P = 61*100000;
O2_span  = [0.21 0.30 0.40 0.50 0.70 1.0];

% Change directory
cd 01b-nonisothermal-constantvolume/

% Initialize variable of interest
igt = zeros(length(O2_span), length(T_span));
max_rr = zeros('like', igt);
Pmax = zeros('like', igt);

% Initialize ROPA dictionaries where all the information are stored
% Questo e un dizionario dove l'elemento (i,j) e la ropa fatta con O2(i) e
% T(j). Le nostre ROPA erano fatte a 1190 K quindi solo per j = 1, ma cosi
% almento ti ci puoi estrarre i reaction rate che ti servono dopo. Vedi
% dopo su come estrarre la roba da li, c'e un esempio di ropa plottata gia
ROPA_dic = cell(length(O2_span), length(T_span));

for i = 1 : length(O2_span)

    for j = 1 : length(T_span)
        
        % Check if the input exist, in that case delete it and write a new
        % one
        if isfile('input.dic')
            system('rm input.dic');
        end

        % Calculate phi
        x = O2_span(i);
        y = (1-x)/x;
        phi = (0.3/9.52)*2*(1+y);
        
        % Write the new input
        output = write_batch_in(O2_span(i), T_span(j), P, phi, 1e-3);
        
        % Run the simulation
        unix(command_run);

        % Load the output results
        data = importdata('Output/Output.out');
        val = data.data;
        
        Tout = val(:,5);
        t = val(:,1);

        O2 = val(end,15);
        H2O = val(end,24);

        dT = zeros(length(Tout)-1,1);
        for l = 1 : length(dT)

            dT(l) = (Tout(l+1) - Tout(l))/(t(l+1) - t(l));

        end

        Pmax(i,j) = max(val(:,6))/100000;
        
        [~,id] = max(dT);
        igt(i,j) = t(id);
        
        [~,idm] = max(val(:,14));
        tropa = t(id-ind(i));

        % Re-Run the simulation 
        % Write the new input
        output = write_batch_in(O2_span(i), T_span(j), P, phi, tropa);
        
        % Run the simulation
        unix(command_run);

        % Extract the ROPA values and data
        species_ropa = {'CH4', 'OH', 'HO2', 'H2O2'};

        for k = 1 : length(species_ropa)

            % Extract ROPA data for the selected species
            [rr, rs] = plot_ropa(species_ropa{k}, 'Output/ROPA.out');

            % Store the information
            ROPA_dic{i,j}.species_ropa{k}.RR = rr;
            ROPA_dic{i,j}.species_ropa{k}.Rnames = rs;
            ROPA_dic{i,j}.T = T_span(j);
            ROPA_dic{i,j}.O2value = O2_span(i);

        end

    end
        
end

cd ../

%% Plot ignition delay time

% Qui puoi cambiare la variabile che vuoi sull'asse X
xax = 'T'; % or O2

if strcmp(xax, 'O2')

    figure;
    leg = cell(length(T_span),1);
    count = 0;

    map = brewermap(length(T_span)+2, '-Greys');
    
    for j = 1 : length(T_span)
        plot(100*O2_span, igt(:,j), 'LineWidth', 2, 'Color', map(j,:), 'Marker', '*', "MarkerSize", 5, 'LineStyle','--');
        hold on;
        leg{j} = append('T ', num2str(T_span(j)), ' K');
        count = count + 1;
    end

    % Legend
    legend(leg, 'Box', 'off', 'Location','northeast');

    % Axis labels
    xlabel('O_2 [mol %]'); ylabel('IDT [s]');

    % Export figure
    if isfolder('NewFigures') == false
        mkdir('NewFigures');
    end
        
    fig = gcf;
    exportgraphics(fig, 'NewFigures/IDT_vs_T.png', 'Resolution',600);

elseif strcmp(xax, 'T')

    figure;
    leg = cell(length(O2_span),1);
    count = 0;

    map = brewermap(length(O2_span)+2, '-Greys');
    
    for j = 1 : length(O2_span)
    
        plot(T_span, igt(j,:), 'LineWidth', 2, 'Color', map(j,:), 'Marker', '*', "MarkerSize", 5, 'LineStyle','--');
        hold on;
        leg{j} = append(num2str(100*O2_span(j)), ' %');
        count = count + 1;
    
    end

    % Adjust legend
    legend(leg, 'Location','northeast');
    legend('boxoff');

    % Axis Labels
    xlabel('T [K]'); ylabel('IDT [s]');

    % Export figure
    if isfolder('NewFigures') == false
        mkdir('NewFigures');
    end
        
    fig = gcf;
    exportgraphics(fig, 'NewFigures/IDT_vs_O2.eps', 'Resolution',600);

end

%% Plot ROPA of CH4 for the main 2 reactions
r2 = {'O+CH4=OH+CH3'};
r1 = {'OH+CH4=H2O+CH3'};

r1_vals = zeros(length(O2_span),1);
r2_vals = zeros(length(O2_span),1);

for i = 1 : length(O2_span)

    % Esempio di come si estraggono nomi e valori
    rnames = ROPA_dic{i,1}.species_ropa{1,1}.Rnames;
    rvals  = ROPA_dic{i,1}.species_ropa{1,1}.RR;
    id1 = 0;
    id2 = 0;

    for j = 1 : length(rnames)

        if strcmp(rnames{j}, r1)
            id1 = j;
        elseif strcmp(r2, rnames{j})
            id2 = j;
        end
    end

    if id1 == 0 || id2 == 0
        error('One of the reaction was not found in the names');
    end

    r1_vals(i) = rvals(id1);
    r2_vals(i) = rvals(id2);
end

% Bar Plot
figure;
map = brewermap(length(O2_span)+2, '-Greys');
b1 = barh([r1_vals'; r2_vals']);

% Modify colors and create legend
leg = cell(length(O2_span),1);
for i = 1 : length(O2_span)
    b1(i).FaceColor = map(i,:);
    leg{i} = append(num2str(100*O2_span(i)), ' %');
end

% Change YTicks
ax = gca; ax.YTick = [];

% Create legend (not flipped)
ll = legend(leg, 'FontSize', 14, 'Location','northwest');
legend('boxoff');
% Flip order
legend(flip(b1), flip(leg));

% Axis label
xlabel('Production rate [kmol/m^3 s]');

% x-axis limits
xlim([-1050 0]);

% Put reactions as text
text(-1000, 1, r1); text(-450, 2, r2);

