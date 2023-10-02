function im_to_gif(filename,im,idx)

for idx=1:idx
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
          imwrite(A,map,filename,'gif', 'Loopcount',inf,'DelayTime',0.5);
          %     elseif idx == 2
%           imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',2);
%     elseif idx == 3
%           imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',2);
    else
          imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);
    end
end
end
