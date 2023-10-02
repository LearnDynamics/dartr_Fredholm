function image_space_fsoi(B,V,eigAB,r,L_operator,y,y_true,titl)
% plot the image of the eigen-vectors

figure; subplot(1,2,1);
y_coordinate_V = V'*L_operator';   % size n_u *n_t
plot(y_coordinate_V(1:r,:)','linewidth',1);
title([titl,' eigen-vectors'])

coordinate_y      = y_coordinate_V*y;
coordinate_y_true = y_coordinate_V*y_true;

subplot(122)
plot(coordinate_y(1:r*2)); hold on; 
plot(coordinate_y_true(1:r*2),'k:','linewidth',2)
title([titl,' coefs'])