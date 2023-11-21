% a way to get intrinsinc (non redundant) coordinates: I use a PCA method
% PCA gives the direction with maximal variance, it is the 1st eigenvector 
% than PCA seeks the 2d direction of max variance, orthogonal to the 1st
% then the 3rd direction is orthogonal to the direction 1 & 2, thus is
% fully constrained by those eigenvectors 1 & 2
% julien Lagarde, Euromov

close all
% example with linear relation x & z
x = randn(100,1);
y = randn(100,1) + 50;
z = 2*x + 2 + randn(100,1);
data = [x y z];% rows = observations (time samples), columns = original coordinates basis (lab frame)
data = zscore(data);
% run the pca
% [eigvect, scores, eigval] = princomp(data);% previous version Matlab
[eigvect, scores, eigval] = pca(data);
%[COEFF,SCORE,latent] is [eigvect, scores, eigval]
% scores = observations projected in the new eigenvect basis
% eigvect = 3 X 3, 1st eigvect = first column = eigvect(:,1)
%
% to visualize the PCA new coordinates
% I plot the 3D data and lines to display the 3 eigenvectors
% 1st eigvect black, 2d green, 3rd cyan
figure(1)
% subplot(2,2,1)
plot3(data(:,1),data(:,2),data(:,3),'o')
grid on
hold on
line([mean(data(:,1)) (eigvect(1,1)) + mean(data(:,1))],...
 [mean(data(:,2)) (eigvect(1,2)) + mean(data(:,1))],...
 [mean(data(:,3)) (eigvect(1,3)) + mean(data(:,1))],...
 'Marker','none','LineStyle','-','LineWidth',2,'Color','k');
hold on
line([mean(data(:,1)) (eigvect(2,1)) + mean(data(:,1))],...
 [mean(data(:,2)) (eigvect(2,2)) + mean(data(:,1))],...
 [mean(data(:,3)) (eigvect(2,3)) + mean(data(:,1))],...
 'Marker','none','LineStyle','-','LineWidth',2,'Color','g');
hold on
line([mean(data(:,1)) (eigvect(3,1)) + mean(data(:,1))],...
 [mean(data(:,2)) (eigvect(3,2)) + mean(data(:,1))],...
 [mean(data(:,3)) (eigvect(3,3)) + mean(data(:,1))],...
 'Marker','none','LineStyle','-','LineWidth',2,'Color','c');
axis equal

% Projections: Is it the "scores" of the output of the princomp? Yes
% (modulo some little gap)
 data_proj = eigvect*data';
 data_proj =data_proj';
 
 figure(2)
 plot3(data_proj(:,1),data_proj(:,2),data_proj(:,3),'*')
 hold on
 plot3(scores(:,1),scores(:,2),scores(:,3),'o')
 grid on