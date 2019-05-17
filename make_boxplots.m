function make_boxplots(C_error, colors, legend_entries, x_pos, bar_width, varargin)
% MAKE_BOXPLOTS A function for making the box plots in our paper. This is a
% work-in-progress, so the function is currently not very well-documented
% at the moment.
%
%   MAKE_BOXPLOTS(C_error, colors, legend_entries, x_pos, bar_width)
%   creates a box plot like those presented in our paper. The function is
%   currently a work-in-progress, so it is best to look at the code below,
%   and the code that utilize this function, to figure out how to use it.

params = inputParser;
addParameter(params, 'reference_line', nan);
parse(params, varargin{:});

reference_line = params.Results.reference_line;

fig_handle = figure;
[no_methods, no_rec, no_trials] = size(C_error);

% Construct data vector
X = C_error(2:end, :, :);
x = [C_error(1, :, 1) X(:)'];

% Construct positions vector
positions =  zeros(1, no_methods*no_rec);
for rec = 1:no_rec
    positions((rec-1)*no_methods+1 : rec*no_methods) = rec + x_pos;
end

group = zeros(size(x));
group(1:no_rec) = 1 + (0:no_rec-1)*no_methods;
repvec = zeros(1, (no_methods-1)*no_rec);
cnt = 1;
for rec = 1:no_rec
    for met = 2:no_methods
        repvec(cnt) = met + (rec-1)*no_methods;
        cnt = cnt + 1;
    end
end
group(no_rec+1:end) = repmat(repvec, 1, no_trials);

% Construct group vector 2
% Old, currently unused code
%{
group2 = zeros(2, size(x, 2));
group2(2, :) = repelem((1:no_rec), 1, (no_methods-1)*no_trials+1);
repvec = 1;
for met = 2:no_methods
    repvec = [repvec met*ones(1, no_trials)]; %#ok<AGROW>
end
group2(1, :) = repmat(repvec, 1, no_rec);
%}

%{
for rec = 1:no_rec
    x((rec-1)*((no_methods-1)*no_trials+1)+1) = C_error(1, rec, 1);
    for met = 2:no_methods
        x(no_rec+1m : )
    end
end
%}

%{
x(4 : no_trials+3) = C_error(2, 1, :);
x(no_trials+4 : 2*no_trials+3) = C_error(2, 2, :);
x(2*no_trials+4 : 3*no_trials+3) = C_error(2, 3, :);
group = [1 3 5 2*ones(1, no_trials) 4*ones(1, no_trials) 6*ones(1, no_trials)];

positions = [1 1.25 2 2.25 3 3.25];
%}

%CC = [0 0 1; 1 0 0; .5 .5 0; 0 0 1];

bplot = boxplot(x, group, 'positions', positions, 'PlotStyle','traditional', 'widths', bar_width, 'colors' , colors, 'medianstyle', 'target');
set(bplot,{'linew'},{.5})

set(gca,'xtick', mean(reshape(positions, no_methods, no_rec)));
x_tick_cell = cell(1, no_rec);
for rec = 1:no_rec
    x_tick_cell{rec} = num2str(rec);
end
set(gca,'xticklabel', x_tick_cell)

%color = repmat(['k'], 1, no_methods*no_rec);
%face_color = {'red', 'none', 'red', 'none', 'red', 'none'};

line_width = repmat([ones(1, no_methods-1) 1], 1, no_rec);
%face_color = repmat(fliplr(face_color), 1, no_rec);
%edge_color = repmat(fliplr(edge_color), 1, no_rec);
%h = findobj(gca,'Tag','Box');
%for j=1:length(h)
%   patch('XData', get(h(j),'XData'), 'YData', get(h(j),'YData'), 'FaceColor', face_color{j}, 'EdgeColor', edge_color{j}, 'FaceAlpha',.5, 'linewidth', line_width(j));
%end

c = get(gca, 'Children');
h = findobj(gca,'Tag','Box');

%hleg1 = legend(c(1:no_methods), legend_entries, 'location', 'northwest');

xlabel('Number of recursions')
ylabel('Error')

delta = 1-(x_pos(end)+bar_width);
for rec = 2:no_rec
    xline(rec-(delta+bar_width)/2, 'linestyle', ':', 'linewidth', 1);
end

if isnan(reference_line)
    legend(h(no_methods:-1:1), legend_entries, 'location', 'northwest')
else
    hline = refline([0 reference_line]);
    hline.Color = 'black';
    hline.LineStyle = '--';
    hline.LineWidth = 2;
    legend_entries{end+1} = 'Standard';
    legend([h(no_methods:-1:1); hline], legend_entries, 'location', 'best')
end

shg

end