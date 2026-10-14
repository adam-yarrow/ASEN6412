# shadedErrorBar

[![View raacampbell/shadedErrorBar on File Exchange](https://camo.githubusercontent.com/bc0a55e41ffd89fcd2884fb1a8b8f7d3eeaf1671627b6c764f3232a28a96addf/68747470733a2f2f7777772e6d617468776f726b732e636f6d2f6d61746c616263656e7472616c2f696d616765732f6d61746c61622d66696c652d65786368616e67652e737667)](https://uk.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar) [![Open in MATLAB Online](https://camo.githubusercontent.com/5047ced9e66f50cb2335e23dcd3242f1ae88a35ab5a32a32ef51d0298f784617/68747470733a2f2f7777772e6d617468776f726b732e636f6d2f696d616765732f726573706f6e736976652f676c6f62616c2f6f70656e2d696e2d6d61746c61622d6f6e6c696e652e737667)](https://matlab.mathworks.com/open/github/v1?repo=raacampbell/shadedErrorBar&file=demo_shadedErrorBar.mlx)

`shadedErrorBar` is a MATLAB function that creates a continuous shaded error region around a line rather than discrete bars. The error region can either be specified explicitly or calculated on the fly based on function handles. The appearance of the plot can be modified via input arguments and/or the handles of the returned plot objects. You will need OpenGL enabled for transparency to work.

# Examples

### Basic example with overlaid points

```
% Prepare data
y=randn(30,80)*5;
x=(1:size(y,2))-40;
yP = sin( linspace(-2*pi,2*pi,length(x)) )*20;
y = bsxfun(@plus,y,yP)+60;

% Make the plot
clf
shadedErrorBar(x,y,{@mean,@std}); 

% Overlay the raw data
hold on
plot(x,y,'.','color',[0.5,0.5,0.95])

grid on
```

[![](https://raw.githubusercontent.com/raacampbell/shadedErrorBar/master/./exampleImages/basic_with_overlay.png)](https://github.com/raacampbell/shadedErrorBar/blob/master/./exampleImages/basic_with_overlay.png)

### Overlaying multiple lines

```
% Prepare data for first line
y=ones(30,1)*x; 
y=y+0.06*y.^2+randn(size(y))*10;

clf
shadedErrorBar(x,y,{@mean,@std},'lineprops','-b','patchSaturation',0.33)

% Overlay second line
hold on
shadedErrorBar(x,2*y+20,{@mean,@std},'lineprops',{'-go','MarkerFaceColor','g'});

%Overlay third line
y=randn(30,80)*5; 
x=(1:size(y,2))-40;
yP = sin( linspace(-2*pi,2*pi,length(x)) )*20;
y = bsxfun(@plus,y,yP)+60;

% Make this line non-transparent
shadedErrorBar(x, y, {@mean,@std}, 'lineprops', '-r','transparent',false,'patchSaturation',0.075)

grid on
```

[![](https://raw.githubusercontent.com/raacampbell/shadedErrorBar/master/./exampleImages/multiple_lines.png)](https://github.com/raacampbell/shadedErrorBar/blob/master/./exampleImages/multiple_lines.png)

### Modifying plot properties post-hoc

```
% Prepare data
y=randn(30,80)*5; 
x=(1:size(y,2));
yP = sin( linspace(-2*pi,2*pi,length(x)) )*20;
y = bsxfun(@plus,y,yP);

% Plot
clf
s = shadedErrorBar(x, y, {@mean,@std}, 'lineprops', '-r');

% Set face and edge properties
set(s.edge,'LineWidth',2,'LineStyle',':')
s.mainLine.LineWidth = 5;
s.patch.FaceColor = [0.5,0.25,0.25];

% Overlay data points post-hoc
hold on
plot(s.mainLine.XData, s.mainLine.YData,'or','MarkerFaceColor','w')

grid on

set(gca,'XTickLabel',[],'YTickLabel',[])
```

[![](https://raw.githubusercontent.com/raacampbell/shadedErrorBar/master/./exampleImages/mod-handles.png)](https://github.com/raacampbell/shadedErrorBar/blob/master/./exampleImages/mod-handles.png)

# Change Log

15tth November 2022

- Allow targing to a specific axis

8th Novemeber 2019

- Merge Octave support from [JarlPed](https://github.com/JarlPed/)

24th November 2017

- Add `patchSaturation` as an optional argument.

6th August 2017

- Parameter/value pairs for input arguments greater than 3.
- Transparency on by default.
- No longer changes the renderer.

# Contributing

If you wish to contribute please fork the repo, check out the `dev` branch and make your changes there, then file a pull request describing what you did.

# Licensing

This code is licensed under the lesser GPL v3. See `LICENSE` file in the repository.