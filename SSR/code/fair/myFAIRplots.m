% =======================================================================================
% function FAIRplots(task,varargin);
% (c) Jan Modersitzki 2008/07/29, see FAIR and FAIRcopyright.m.
%
% FAIRplots: generates FAIR plots;
% the objective is to create the following 2-by-3 plot:
% _________________________________________________________________
% |                    |                    |                     |
% | titleR0            | titleT0            |  title(Tk)          |
% | showImage(R0)      | showImage(T0)      |  showImage(Tk)      |
% |                    |                    |                     |
% |_______________________________________________________________|
% |                    |                    |                     |
% | titleGrid          | titleT0            |  title(Tk)          |
% | showImage(T0)      | showDifference(D0) |  showDifference(Tk) |
% | and showGrid       |                    |                     |
% |_______________________________________________________________|
%
% the function uses persitent variables:
%   fig, plots, figname, gridHandle, 
%   showImage, showDifference, showGrid,
%   Rname, Tname, Dname, Gname
% and procides four tasks:
%
% - FAIRplots('clear')
%   clears all persistent variables and returns
% - FAIRplots('set',varargin)
%   initializes fig, plots, figname (all to empty), 
%   showImage=@viewImage, showDifference=@viewImage(255-abs(T-R)),
%   vecNorm=@norm; and possibly overwrites by the varargin list
%   initializes also gridHandle=[], showGrid (dimension dependent), and the 
%   name of the calling function, create a status line to be used as figname 
%   and the titles for the reference, template, distance, and grid
%
%   Examples: 
%
%   para = {Tc,Rc,omega,m,yc,normdy,Jc} a collection of the current data
%
%   set-up plots:   FAIRplots('set','mode','NPIR','fig',level);
%                   note no plots are shown using set
%
%   init plots:     FAIRplots('init',para);
%                   shows T and R
%                   R(xc)            -> (2,3,1) title: Rname
%                   T(xc)            -> (2,3,2) title: 'T(xc)'
%                   T(xc)            -> (2,3,4) title: 'T(xc)'
%                   T(xc)-R(xc)      -> (2,3,5) title: '|T(xc)-R|'
%
%   show stopping:  FAIRplots('stop',para);
%                   where -inf is just a dummy number and 
%                   T(yStop)         -> (2,3,3) title: 'Tstop'
%                   T(yStop)-R(xc)   -> (2,3,6) title: '|Tstop-R|=100%'
%                   yStop            -> (2,3,4) (replace)
%
%   show starting:  FAIRplots('start',para);
%                   T(y0)            -> (2,3,2) title: Tname
%                   T(y0)-R(xc)      -> (2,3,5) title: Dname
%                   y0               -> (2,3,4) (replace)
%
%   show current:   FAIRplots(iter,para);
%                   T(yc)            -> (2,3,3) title: Tanme
%                   T(yc)-R(xc)      -> (2,3,6) title: Dname
%                   yc               -> (2,3,4) (replace)
%
% =======================================================================================

function myFAIRplots(task,varargin)

% handle options
persistent fig plots figname gridHandle omega m 
persistent showImage showDifference showGrid showMI
persistent Rname Tname Dname Gname Jstop

% subplot
columns = 4;
plotIndex = @(r,c) (r-1)*columns + c;
MIname = @(iter) sprintf('MI(%d)',iter);

if strcmp(task,'set') || strcmp(task,'reset'),
  if strcmp(task,'reset'),
     fig = []; plots = []; figname = [];  gridHandle = [];  omega = [];  m = []; 
     showImage = [];  showDifference = [];  showGrid = []; showMI = [];
     Rname = [];  Tname = [];  Dname = [];  Gname = [];  Jstop = []; 
  end;
  

  % initialize fig, plots, figname, omega, m to []
  if ~exist('fig','var')      || isempty(fig),      fig     = []; end;
  if ~exist('plots','var')    || isempty(plots),    plots   = []; end;
  if ~exist('figname','var')  || isempty(figname),  figname = []; end;

  % initialize showImage to viewImage
  if ~exist('showImage','var')|| isempty(showImage),
    showImage = @viewImage;
  end;
  % initialize showDifference to viewImage(255-abs(T-R))
  if ~exist('showDifference','var') || isempty(showDifference),
    showDifference = @(T,R,omega,m) viewImage(255-abs(T-R),omega,m);
  end;
  % initialize showMI to MIcc
  if ~exist('showMI','var') || isempty(showMI),
    showMI = @(T,R,omega,m) MIcc(T,R,omega,m);
  end;

  mode    = 'MATLAB';
  
  for k=1:2:length(varargin), % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
  end;

  % initial gridhandle to [], and showGrid to plotGrid (dim==2) or none (dim==3)
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  gridHandle      = [];
  showGrid = @(yc,omega,m)  plotGrid(yc,omega,m,'spacing',ceil(m/32));
  

  % initialze figname and Rname depending on PIR*/NPIR*
  if length(mode) > 2 && strcmp(mode(1:3),'PIR'),
    w0 = trafo('w0');
    figname = @(omega,m)sprintf('%s: %s/%s/%s, %dD, m=%s',...
      mode,inter,distance,trafo,length(omega),dimstr(m));
    Rname = @(m)sprintf('R, %s, length(w)=%d',dimstr(m),length(w0));
  elseif  length(mode)>3 && strcmp(mode(1:4),'NPIR')
    figname = @(omega,m) sprintf('%s: %s/%s/%s, %dD, m=%s',...
      mode,inter,distance,regularizer,length(omega),dimstr(m));
    Rname = @(m) sprintf('R, %s, \\alpha=%s',...
      dimstr(m),num2str(regularizer('get','alpha')));
  else
    figname = @(omega,m) sprintf('%s: %dD, m=%s',...
      mode,length(omega),dimstr(m));
    Rname = @(m) sprintf('R, %s',dimstr(m));
  end;

    
  Tname = @(iter) sprintf('T(%d)',iter);
  Dname = @(j,Jc,Jstop) sprintf('|J(%d)/Jstop|=%s%%',j,num2str(abs(100*Jc/Jstop)));
  Gname = @(normdY) sprintf('T(xc), |dY|= %s',num2str(normdY));

  return;
end;

if ~plots, return; end;
if nargin>1, 
  T      = getField(varargin{1},'Tc'); 
  R      = getField(varargin{1},'Rc'); 
  omega  = getField(varargin{1},'omega'); 
  m      = getField(varargin{1},'m'); 
  yc     = getField(varargin{1},'yc');
  normdY = getField(varargin{1},'normdY');
  Jc     = getField(varargin{1},'Jc');
end;

switch task,
  case 'clear', % clear persistent variables and return
    clear fig plots figname gridHandle 
    clear showImage showDifference showGrid
    clear Rname Tname Dname Gname
    return;
    
  case 'init',
    % activate figure, set colordef and figname
    if ~ishandle(fig), fig = figure; else set(0,'CurrentFigure',fig); end;
    
    % extract variables
    xc = getCenteredGrid(omega,m);
    T  = inter(T,omega,xc);
    R  = inter(R,omega,xc);

    FAIRfigure(fig,'figname',figname(omega,m),'position','default');

    if length(omega) == 6, % disable showGrid
      showGrid = @(yc,omega,m) [];
    end;

    % plot
    % ______________________________________
    % | R       | T     |         | MI(0)  |
    % |_________|_______|_________|________|
    % | T +grid | T -R  |         |        |
    % |_________|_______|_________|________|
    if ~ishandle(fig), fig = figure; else set(0,'CurrentFigure',fig); end; clf;
    subplot(2,columns,plotIndex(1,1)); showImage(R,omega,m);           title(Rname(m));
    subplot(2,columns,plotIndex(1,2)); showImage(T,omega,m);           title('T(xc)');
    subplot(2,columns,plotIndex(2,1)); showImage(T,omega,m);           title('T(xc)'); hold on;
    subplot(2,columns,plotIndex(2,2)); showDifference(T,R,omega,m);    title('|T(xc)-R|');
    % MI
    [Dc,rho,dD,drho,d2psi,nT,nR] = MIcc(T,R,omega,m);
    logrho=log2(rho);
    logrho(isinf(logrho))=0;
    logrho(logrho~=0)=logrho(logrho~=0)-min(logrho(:));
    subplot(2,columns,plotIndex(1,4)); viewImage(logrho,[0 nT 0 nR],[nT nR],'scale',1); title(MIname(0));
    
    pause(1/100); drawnow;

  case 'stop',
    % plot
    % ______________________________________
    % |         |       | Tstop   |        |
    % |_________|_______|_________|________|
    % |   +grid |       | Tstop-R |        |
    % |_________|_______|_________|________|
    if ~ishandle(fig), fig = figure; else set(0,'CurrentFigure',fig); end;
    subplot(2,columns,plotIndex(1,3)); showImage(T,omega,m);     title('T^{stop}');
    subplot(2,columns,plotIndex(2,3)); showDifference(T,R,omega,m);
    title('|T^{stop}-R|, J^{stop}=100%');
    subplot(2,columns,plotIndex(2,1)); set(gridHandle,'visible','off');
    gridHandle = showGrid(yc,omega,m);
    pause(1/100); drawnow;
    Jstop = Jc;
    
  case 'start',
    % plot
    % ______________________________________
    % |         | T0    |         |        |
    % |_________|_______|_________|________|
    % |   +grid | T0-R  |         |        |
    % |_________|_______|_________|________|
    if ~ishandle(fig), fig = figure; else set(0,'CurrentFigure',fig); end;
    subplot(2,columns,plotIndex(1,2)); %{showImage(T,omega,m);
                            title(Tname(0));
    subplot(2,columns,plotIndex(2,2)); showDifference(T,R,omega,m);   
    title(Dname(0,Jc,Jstop));
    subplot(2,columns,plotIndex(2,1)); set(gridHandle,'visible','off');
    gridHandle = showGrid(yc,omega,m);
    pause(1/100); drawnow;

  otherwise, 
    if ~isnumeric(task),
      warning(['don''t no how to deal task <',task,'>!']);
      return;
    end;
    
    % plot
    % ______________________________________
    % |         |       | Tc      |        |
    % |_________|_______|_________|________|
    % |   +grid |       | Tc-R    | MI     |
    % |_________|_______|_________|________|
    set(0,'CurrentFigure',fig); % don't steel focus
    subplot(2,columns,plotIndex(2,1)); set(gridHandle,'visible','off');
    gridHandle = showGrid(yc,omega,m); 
    title(Gname(normdY));
    subplot(2,columns,plotIndex(1,3)); showImage(T,omega,m);             title(Tname(task));
    subplot(2,columns,plotIndex(2,3)); showDifference(T,R,omega,m);      
    title(Dname(task,Jc,Jstop));
    % MI
    [Dc,rho,dD,drho,d2psi,nT,nR] = MIcc(T,R,omega,m);
    logrho=log2(rho);
    logrho(isinf(logrho))=0;
    logrho(logrho~=0)=logrho(logrho~=0)-min(logrho(:));
    subplot(2,columns,plotIndex(2,4)); viewImage(logrho,[0 nT 0 nR],[nT nR],'scale',1); title(MIname(task));
    
    drawnow;    
end;

function v = getField(s,field);
if isfield(s,field), v = s.(field); else v = []; end;

%$=======================================================================================
%$  FAIR: Flexible Algorithms for Image Registration
%$  Copyright (c): Jan Modersitzki
%$  1280 Main Street West, Hamilton, ON, Canada, L8S 4K1
%$  Email: modersitzki@mcmaster.ca
%$  URL:   http://www.cas.mcmaster.ca/~modersit/index.shtml
%$=======================================================================================
%$  No part of this code may be reproduced, stored in a retrieval system,
%$  translated, transcribed, transmitted, or distributed in any form
%$  or by any means, means, manual, electric, electronic, electro-magnetic,
%$  mechanical, chemical, optical, photocopying, recording, or otherwise,
%$  without the prior explicit written permission of the authors or their 
%$  designated proxies. In no event shall the above copyright notice be 
%$  removed or altered in any way.
%$ 
%$  This code is provided "as is", without any warranty of any kind, either
%$  expressed or implied, including but not limited to, any implied warranty
%$  of merchantibility or fitness for any purpose. In no event will any party
%$  who distributed the code be liable for damages or for any claim(s) by 
%$  any other party, including but not limited to, any lost profits, lost
%$  monies, lost data or data rendered inaccurate, losses sustained by
%$  third parties, or any other special, incidental or consequential damages
%$  arrising out of the use or inability to use the program, even if the 
%$  possibility of such damages has been advised against. The entire risk
%$  as to the quality, the performace, and the fitness of the program for any 
%$  particular purpose lies with the party using the code.
%$=======================================================================================
%$  Any use of this code constitutes acceptance of the terms of the above
%$                              statements
%$=======================================================================================
