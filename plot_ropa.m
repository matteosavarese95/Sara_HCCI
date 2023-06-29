function [rr, rs] = plot_ropa(sp_name, ropa_file)
% sp_name = string of the species
% ropa_file = 'ROPA.out'

cf = readlines(ropa_file);

rr = {}; % Value of reaction rate
rs = {}; % Reaction name

for i = 1 : length(cf)

    cfi = split(cf(i));

    if cfi(1) == sp_name && i > 180

        disp('Extracting data...');
        count = 2;
        while strcmp(cf(i+count), '') == false

            rstring = split(cf(i+count));

            rr{count-1} = str2double(rstring(2));
            rs{count-1} = rstring(7);

            count = count + 1;

        end

    end

end

rr = cell2mat(rr);

% Choose to plot
plotting = false;

if plotting == true
    h = barh(rr, 0.4); hold on;
    
    for j = 1 : length(rr)
    
        if rr(j) > 0
            text(-rr(j), j, rs{j}, 'FontSize',12);
    
        else
            text(0.1*rr(j), j, rs{j}, 'FontSize', 12);
        end
    end
    
    if max(rr) > abs(min(rr))
        xlim([-max(rr)-0.1*max(rr) max(rr)+0.1*max(rr)]);
    else
        xlim([min(rr)-0.1*(-min(rr)) -min(rr) + 0.1*(-min(rr))]);
    end
    
    tit = append('ROPA ', sp_name);

end







        





end

