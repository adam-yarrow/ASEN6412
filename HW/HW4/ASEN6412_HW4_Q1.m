%{
    Running this function generates the tensor grids associated with the
    generation of the smolyak sparse CC grid from HW4 Q1.

    NOTE: The library path contains the codes posted on CANVAS.
%}

%% Plots tensor product grids
addpath(genpath("..\..\Library\"))

%% Make 1D Grids
level = [0, 1, 2, 3, 4];
x = cell(5, 1);
for i = 1:5
    if level(i) == 0
        Nl = 1;
    else
        Nl = 2^level(i) + 1;
    end
    x{i} = fclencurt(Nl, -1, 1);
end

%% Make tensor products
tensorProdY1 = [3, 2, 1, 0, 4, 3, 2, 1, 0];
tensorProdY2 = [0, 1, 2, 3, 0, 1, 2, 3, 4];
smolynakWeights = [-1, -1, -1, -1, 1, 1, 1, 1, 1];

tensorProdStoragePos = [];
tensorProdStorageNeg = [];

for prodI = 1:numel(tensorProdY1)
    l1 = tensorProdY1(prodI);
    l2 = tensorProdY2(prodI);
    
    name = sprintf('I%i x I%i', l1, l2);
    
    [y1, y2] = meshgrid(x{l1+1}, x{l2+1});
    if smolynakWeights(prodI) < 0
        tensorProdStorageNeg = vertcat(tensorProdStorageNeg,...
                                        [y1(:), y2(:)]);
    else
        tensorProdStoragePos = vertcat(tensorProdStoragePos,...
                                        [y1(:), y2(:)]);
    end

    % Process two separate batches
    if prodI == 1
        figure('Name','|l| = 3 Tensor Products');
        iCntr = 1;
        maxRows = 2;
    elseif prodI == 5
        figure('Name','|l| = 4 Tensor Products');
        iCntr = 1;
        maxRows = 3;
    end
    subplot(maxRows, 2, iCntr);
    plot(y1, y2, 'r.', 'MarkerSize',10);
    grid on;
    xlabel('y_1');
    ylabel('y_2');
    title(name);
    xlim([-1,1]);
    ylim([-1,1]);

    iCntr = iCntr + 1;
end

%% Check full product correct
finalGrid = tensorProdStoragePos;
tol = 1E-8;
for iNeg = 1:size(tensorProdStorageNeg,1)
    foundIdx = tolerancedFind(finalGrid, tensorProdStorageNeg(iNeg,:), tol);

    if ~isempty(foundIdx)
        finalGrid(foundIdx, :) = [];
    else
        error('Cant subtract something that doesnt exist');
    end
end

figure('Name','Verify Solution')
plot(finalGrid(:,1), finalGrid(:,2), 'r.', 'MarkerSize',10);
grid on;
xlabel('y_1');
ylabel('y_2');
title('Combined Sparse Grid');
xlim([-1,1]);
ylim([-1,1]);

%% helper function
function foundIdx = tolerancedFind(src, target, tol)
    % Finds target idx within src vector 
    % Assumes src = N x 2
    % Assumes target = 1 x 2

    srcDiff = NaN(size(src));
    srcDiff(:,1) = abs(src(:,1) - target(1));
    srcDiff(:,2) = abs(src(:,2) - target(2));

    condition = srcDiff(:,1) < tol & srcDiff(:,2) < tol;
    foundIdx = find(condition, 1,"first");
end

