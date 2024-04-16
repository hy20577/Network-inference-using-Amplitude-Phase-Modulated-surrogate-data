function compactheatmap_fromTFPR(truePositiveRate,falsePositiveRate)

% Calculate the size scaling factor for the circles
maxCircleSize = 1e3; % You can adjust this based on your preference
circleSize = maxCircleSize * falsePositiveRate;

% Create the heatmap of true positive rate
imagesc(truePositiveRate);
xlabel('$pc_{2}$', 'Interpreter', 'latex');
ylabel('$pc_{1}$', 'Interpreter', 'latex');
% Set the font size for x-axis and y-axis labels

c = colorbar('Location', 'eastoutside');
c.Label.String = 'TPR'; % Set the colorbar label
c.Label.FontSize = 30;
c.Label.Interpreter = 'latex';
clim([0, 1]);
[numRows, numCols] = size(truePositiveRate);

hold on;

% Loop through each cell in the matrix and plot circles
for row = 1:numRows
    for col = 1:numCols
        centerX = col;  % X-coordinate of the cell center
        centerY = row;  % Y-coordinate of the cell center
        radius = circleSize(row, col); % Circle radius
        
        if radius > 0
            scatter(centerX, centerY, radius, 'r', 'filled');
        end
    end
end

xticks(1:numCols);
xticklabels(num2str(linspace(0, 100, numCols)', '%d'));

yticks(1:numRows);
yticklabels(num2str(linspace(100, 0, numRows)', '%d'));


% Set square aspect ratio
axis equal;

for row = 1:numRows
    for col = 1:numCols
        if all([truePositiveRate(row, col) == 1, falsePositiveRate(row,col) == 0])
            centerX = col; % + 0.5; % X-coordinate of the cell center
            centerY = row; %+ 0.5; % Y-coordinate of the cell center
            radius = 1e2; % Circle radius
                if radius > 0
                    scatter(centerX, centerY, radius, 'g', 'filled');
                end
        end
    end
end

hold off;

ax = gca; 
ax.FontSize =30;

% Adjust axis limits and labels
axis([0.5 numCols+0.5 0.5 numRows+0.5]);

lg = legend('FPR', 'Location', 'northoutside');
lg.Interpreter = 'latex';
set(lg, 'Fontsize', 30);

legend boxoff 

end