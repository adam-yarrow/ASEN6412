%% Plots tensor product grids
addpath("..\..\ProfCodes\")

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

%% TODO - verify if time
% totalY1 = [];
% totalY2 = [];
% for i = 9:-1:1
%     totalY = tensorProdStorage
% end

