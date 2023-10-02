%% set position, size of figure; and set fonts
%  function set_positionFontsAll(nrow,ncol)
ftsz =15; mksz = 10; % if save to -dpdf
% ftsz =12; mksz = 8; % if save to -depsc
set(findall(gcf,'-property','FontSize'),'FontSize',ftsz);
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',mksz);
% tightfig;

%% get the number of subplots in gcf; 
cells_axes = findobj(gcf,'type','axes');
N = length(cells_axes);
pos1 = zeros(N,1); pos2 = zeros(N,1);
for n = 1:N
    pos1(n) = cells_axes(n).Position(1);
    pos2(n) = cells_axes(n).Position(2);
end
Ncols = numel(unique(pos1));
Nrows = numel(unique(pos2));

%% set figure size
width = 500+300*(Ncols-1); height = 300+200*(Nrows-1); 
set(gcf, 'Position',  [100, 1000, width, height]);    % set figure position + size [x , x, width height]
tightfig;