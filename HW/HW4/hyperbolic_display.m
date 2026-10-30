function hyperbolic_display ( family, growth, level )

%*****************************************************************************80
%
%% hyperbolic_display() illustrates the exactness of a 2D Smolyak sparse grid.
%
%  Discussion:
%
%    This version is for the hyperbolic-cross growth rate.
%
%    The underlying quadrature family is indexed as 
%      Q(i) = quadrature rule with i+1 point.
%
%    Thus, 
%      Q(0) is the 1 point rule, 
%      Q(1) is the 2 point rule, and so on.
%
%    Assuming a symmetric interval (as for 'CC' and 'L' families), 
%    the exactness of rule Q(i) is 
%      p(i) = i for even i (and odd order), 
%      p(i) = i-1 for odd i (and even order).
%
%    Thus,
%      Q(0) is the 1 point rule of exactness 1,
%      Q(1) is the 2 point rule of exactness 1,
%      Q(2) is the 3 point rule of exactness 3,
%      Q(3) is the 4 point rule of exactness 3, and so on.
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    11 March 2014
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    string FAMILY, indicates the 1D quadrature family used
%    to provide factors, default is 'cc';
%    'CC': Clenshaw-Curtis;
%    'L':  Gauss-Legendre
%
%    string GROWTH, indicates the growth rate, default is 'e'.
%    'E': exponential, for CC, L;
%    'HC': Hyperbolic Cross;
%    'L': linear growth, for CC, HC, L;
%    'LO', linear odd growth, for CC, L;
%    'SE': slow exponential, for CC, L.
%
%    integer LEVEL, the level of the sparse grid to be constructed.
%    0 <= LEVEL <= 5; default is 3.
%
  if ( upper ( growth ) ~= 'HC' )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'HYPERBOLIC_DISPLAY - Fatal error!\n' );
    fprintf ( 1, '  Only "HC" growth is legal for this function.\n' );
    fprintf ( 1, '  Your growth rate was "%s"\n', growth );
    error ( 'HYPERBOLIC_DISPLAY - Fatal error!' );
  end
%
%  Determine the total-degree exactness of this rule.
%
  if ( upper ( family ) == 'CC' )

    p = level;
    for a = 1 : level
      b = floor ( level / a );
      p = min ( p, a + b );
    end

  elseif ( upper ( family ) == 'L' )

    p = 1;
    for a = 1 : level
      pa = exactness_1d ( family, growth, a );
      b = floor ( level / a );
      pb = exactness_1d ( family, growth, b );
      if ( a == 1 )
        p = pa + pb;
      else
        p = min ( p, pa + pb );
      end
    end

  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Sparse grid exactness will be %d\n', p );
  fprintf ( 1, '\n' );
%
%  Request number of points in 1D rule of given level and growth rate.
%  This is the maximum number of points we will see.
%
  n1d_max = order_1d ( family, growth, level );
%
%  What is the maximum precision of the 1D rules we will use?
%
  p_max = exactness_1d ( family, growth, level );
%
%  Set GRID_MAX.
%
  grid_max = max ( p, p_max ) + 1;
%
%  Ready to begin plotting.
%
  clf

  n1 = grid_max;
  n2 = grid_max;

  g = zeros ( n1 + 1, n2 + 1 );

  p1 = n1;
  p2 = n2;
%
%  Make the plot.
%
  x_data_min = 0;
  x_data_max = n1;
  y_data_min = 0;
  y_data_max = n2;

  x_plot_min = x_data_min - 0.5;
  x_plot_max = x_data_max + 0.5;

  y_plot_min = y_data_min - 0.5;
  y_plot_max = y_data_max + 0.5;

  x_plot_range = x_plot_max - x_plot_min;
  y_plot_range = y_plot_max - y_plot_min;

  margin = 0.025 * max ( x_plot_range, y_plot_range );

  x_axes_min = x_plot_min - margin;
  x_axes_max = x_plot_max + margin;
  y_axes_min = y_plot_min - margin;
  y_axes_max = y_plot_max + margin;
%
%  Every box starts out WHITE.
%
  for x = 0 : n1

    for y = 0 : n2

      a = x - 0.44;
      b = x + 0.44;
      c = y - 0.44;
      d = y + 0.44;

      patch ( [ a, b, b, a ], [ c, c, d, d ], 'w', 'EdgeColor', 'black' );

    end
    
  end

  diag = 0;
  s = + 1;

  lp = level;
 
  lp1 = 0;
  lp1_save = lp1;
  lp2 = level;
 
  while ( lp1 <= level )

    if ( lp ~= 0 )

      while ( lp2 * ( lp1 + 1 ) <= level )
        lp1 = lp1 + 1;
        lp2 = floor ( level / lp1 );
      end

    end

    if ( diag == 1 )
      lp1 = lp1_save;
    end

    level1 = lp1 - 1;
    level2 = lp2 - 1;
    p1 = exactness_1d ( family, growth, level1 );
    p2 = exactness_1d ( family, growth, level2 );

    g(0+1:p1+1,0+1:p2+1) = g(0+1:p1+1,0+1:p2+1) + s;

    for x = 0 : p1
      for y = 0 : p2

        a = x - 0.44;
        b = x + 0.44;
        c = y - 0.44;
        d = y + 0.44;

        if ( g(x+1,y+1) <= -3 )
          rgb = [ 1.0, 0.0, 0.0 ];
        elseif ( g(x+1,y+1) <= -2 )
          rgb = [ 1.0, 0.2, 0.2 ];
        elseif ( g(x+1,y+1) <= -1 )
          rgb = [ 1.0, 0.4, 0.4 ];
        elseif ( g(x+1,y+1) < 0 )

          rgb = [ 1.0, 0.6, 0.6 ];

        elseif ( g(x+1,y+1) == 0 )
          if ( x + y <= p )
            rgb = [ 0.9, 0.9, 0.9 ];
          else
            rgb = [ 1.0, 1.0, 1.0 ];
          end

        elseif ( g(x+1,y+1) == 1 )
  
          if ( x + y <= p )
            rgb = [ 0.0, 0.0, 1.0 ];
          else
            rgb = [ 0.5, 0.5, 1.0 ];
          end

        elseif ( g(x+1,y+1) <= 2 )
          rgb = [ 0.5, 1.0, 0.5 ];
        elseif ( g(x+1,y+1) <= 3 )
          rgb = [ 0.3, 1.0, 0.3 ];
        elseif ( g(x+1,y+1) <= 4 )
          rgb = [ 0.1, 1.0, 0.1 ];
        else
          rgb = [ 0.0, 1.0, 0.0 ];
        end

        patch ( [ a, b, b, a ], [ c, c, d, d ], rgb, 'EdgeColor', 'none');

      end
    end

    line ( [-0.5, p+0.5 ], [ p+0.5, -0.50 ], 'LineWidth', 3, 'Color', 'k' )

    if ( s == +1 )
      title ( sprintf ( 'Add Q^{%d} x Q^{%d}, orders 0:%dx0:%d', ...
        lp1, lp2, p1, p2 ), 'Fontsize', 16 );
      filename = sprintf ( '%s_%s_%d_+%dx%d.png', ...
        family, growth, level, lp1, lp2 );
    else
      title ( sprintf ( 'Subtract Q^{%d} x Q^{%d}, orders 0:%dx0:%d', ...
        lp1, lp2, p1, p2 ), 'Fontsize', 16 );
      filename = sprintf ( '%s_%s_%d_-%dx%d.png', ...
        family, growth, level, lp1, lp2 );
    end

    axis equal
    axis square
    axis tight

    print ( '-dpng', filename );
    fprintf ( '  Graphics saved as "%s"\n', filename );

    pause ( 5 )

    if ( diag == 0 )
      lp1_save = lp1;
      diag = 1;
      s = -1;
    else
      lp1 = lp1_save;
      diag = 0;
      s = +1;
    end

    lp1 = lp1 + 1;
    lp2 = floor ( level / lp1 );

  end

  return
end
