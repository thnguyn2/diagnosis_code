function F=makeGaussianFilters
  SUP=49;                 % Support of the largest filter (must be odd)
  SCALEX=sqrt(2).^[1:3];  % Sigma_{x} for the oriented filters
  NORIENT=8;              % Number of orientations

  NF=NORIENT*length(SCALEX);%Number of filter
  F=zeros(SUP,SUP,NF);
  hsup=(SUP-1)/2;
  [x,y]=meshgrid([-hsup:hsup],[hsup:-1:-hsup]);
  orgpts=[x(:) y(:)]'; %Coodinates of the orignial points

  count=1;
  for scale=1:length(SCALEX),
    for orient=0:NORIENT-1,
      angle=pi*orient/NORIENT;  % Not 2pi as filters have symmetry
      c=cos(angle);s=sin(angle);    %Rotation matrix
      rotpts=[c -s;s c]*orgpts;     %Coordinates of the new points
      F(:,:,count)=makefilter(SCALEX(scale),rotpts,SUP);
      count=count+1;
    end;
  end;
return

function f=makefilter(scale,pts,sup)
  x=pts(1,:);
  gx = exp(-x.^2/(2*(3*scale)^2));
  y=pts(2,:);
  gy = exp(-y.^2/(2*(scale)^2));
  f=(reshape(gx.*gy,sup,sup));
  f=f/sum(abs(f(:)));
return
