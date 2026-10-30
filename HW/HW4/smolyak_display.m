function smolyak_display ( family, growth, level )

%*****************************************************************************80
%
%% smolyak_display() illustrates the exactness of a 2D Smolyak sparse grid.
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    22 April 2019
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
%    int LEVEL, the level of the sparse grid to be constructed.
%    0 <= LEVEL <= 5; default is 3.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'SMOLYAK_DISPLAY():\n' );
  fprintf ( 1, '  Illustrate the "construction" of the exactness diagram\n' );
  fprintf ( 1, '  of a Smolyak grid.\n' );

  if ( nargin < 1 )
    family = 'CC';
  end

  if ( nargin < 2 )
    growth = 'E';
  end

  if ( nargin < 3 )
    level = 3;
  elseif ( isstring ( level ) )
    level = s2num ( level );
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  FAMILY = %s\n', family );
  fprintf ( 1, '  GROWTH = %s\n', growth );
  fprintf ( 1, '  LEVEL =  %d\n', level );

  if ( growth == 'HC' )
    hyperbolic_display ( family, growth, level );
  else
    nonhyperbolic_display ( family, growth, level );
  end
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'SMOLYAK_DISPLAY\n' );
  fprintf ( 1, '  Normal end of execution.\n' );

  return
end

