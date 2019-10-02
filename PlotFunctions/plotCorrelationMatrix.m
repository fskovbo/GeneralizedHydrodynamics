function plotCorrelationMatrix(x_grid, h_j, h_i, corr_mat, x_margin, y_margin)

    if nargin < 6
        x_margin = 0;
        y_margin = 0;
    end

    h1 = subplot(3,3,[1 4]); % x plot
    plot(x_grid, h_i, 'k', 'Linewidth', 1.5)
    xlim([x_grid(1), x_grid(end) ])
    yl1 = ylim;
    xlabel('x')
    ylabel('hi')
    view([90 -90])
    
    
    h2 = subplot(3,3,[8 9]); % y plot
    plot(x_grid, h_j, 'k', 'Linewidth', 1.5)
    xlim([x_grid(1), x_grid(end) ])
    yl2 = ylim;
    xlabel('y')
    ylabel('hj')
    
    h3 = subplot(3,3,[ 2 3 5 6]); % correlations
    imagesc(x_grid, x_grid, corr_mat)
    set(gca,'YDir','normal')
    set(gca,'XTickLabels',[],'YTickLabels',[])
    xlim([x_grid(1), x_grid(end) ])
    ylim([x_grid(1), x_grid(end) ])
    
 
    % set equal ylimits for quantity plots 
    ylimits = [ min([yl1, yl2]) , 1.1*max([yl1, yl2])];
    set(h1,'ylim',ylimits)
    set(h2,'ylim',ylimits)
    
    set(findall(gcf,'-property','FontSize'),'FontSize',10)
    
    % calculate new positions
    p1 = h1.Position;
    p2 = h2.Position;
    p3 = h3.Position;
    
    p1_new = p1;
    p2_new = p2;
    p3_new = p3;
    
    p1_new(2) = p2(2) + p2(4) + y_margin;
    p1_new(4) = p1(4) + p1(2) - p1_new(2);
    
    p2_new(1) = p1(1) + p1(3) + x_margin;
    p2_new(3) = p2(3) + p2(1) - p2_new(1);
    
    p3_new(1) = p2_new(1);
    p3_new(2) = p1_new(2);
    p3_new(3) = p2_new(3);
    p3_new(4) = p1_new(4);
    
    set(h1, 'Position', p1_new ) % [x0 y0 width height]
    set(h2, 'Position', p2_new ) % [x0 y0 width height]
    set(h3, 'Position', p3_new ) % [x0 y0 width height]
    
    %
        
end