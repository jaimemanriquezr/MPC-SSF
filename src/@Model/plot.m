function graphPlot = plot(obj, options)
    arguments
        obj  Model
        options.Layout = 'layered';
        options.LineWidth = 1.5;
        options.EdgeColor = 'k';
        options.EdgeFontSize = 10;
        options.ArrowSize = 12;
        options.NodeFontWeight = 'bold';
    end
    reactant = obj.Order;
    sigma = obj.StoichiometricCoefficients;
    
    numComponents = length(obj.Components);
    numReactions = length(obj.Reactions);
    
    edgeWeights = [];
    edgeNodes = [];
    edgeLabels = [];
    for i = 1:numReactions
        startNode = mod(find(reactant(:, i)) - 1, numComponents) + 1;
        [endNode, ~, value] = find(sigma(:, i));
        direction = sign(value) * abs(sign(sigma(startNode, i)));
        if norm(direction) == 0
            startNode = find(sigma(:, i) == -1);
            endNode = find(sigma(:, i) == 1);
            value = sigma(startNode, i);
            direction = 1;
        end
        
        % flip edge if direction=-1
        edge = [(1 + direction)/2.*startNode + (1 - direction)/2.*endNode, ...
                (1 - direction)/2.*startNode + (1 + direction)/2.*endNode];
        label = repmat(obj.Reactions(i).Name, [size(edge, 1), 1]);

        edgeWeights = cat(1, edgeWeights, value);
        edgeLabels = cat(1, edgeLabels, label);
        edgeNodes = cat(1, edgeNodes, edge);
    end
    nodeNames = [obj.Components.Name].';
    nodeColors = ones(length(nodeNames), 3);
    nodeColors(ismember(nodeNames, [obj.Particles.Name]), :) = repmat([1 0 0], [length(obj.Particles), 1]);
    nodeColors(ismember(nodeNames, [obj.Liquids.Name]), :) = repmat([0 0 1], [length(obj.Liquids), 1]);

    edgeTable = table(edgeNodes, edgeWeights, edgeLabels, 'VariableNames', {'EndNodes', 'Weight', 'Label'});
    nodeTable = table((1:numComponents).', [obj.Components.Name].', 'VariableNames', {'EndNodes', 'Name'});
    reactionGraph = digraph(edgeTable, nodeTable);

    plotOptions = namedargs2cell(options);
    graphPlot = plot(reactionGraph, 'EdgeLabel', reactionGraph.Edges.Label, 'NodeColor', nodeColors, plotOptions{:});
end