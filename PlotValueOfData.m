function PlotValueOfData( rho, actualRelative, predicted )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

colors = {'b','g','r','c','m','y'};
actLine = ':';
predLine = '-';

[n,l] = size(predicted);

figure(1)
clf
hold on
for ii=1:n
    posVals = predicted(ii,:) > 0;
    plot(rho(posVals), 100*actualRelative(ii,posVals), strcat(colors{ii},predLine), ...
        'LineWidth',3)
end
yl = ylim;
yl(1) = min(0,yl(1));

for ii=1:n
    plot(rho, 100*actualRelative(ii,:), strcat(colors{ii},actLine), ...
        'LineWidth',4, 'MarkerSize',12)
end
ylim(yl)

xlabel('\rho', 'FontSize',14)
ylabel('Relative Decrease From Observation (%)', 'FontSize',14)

hold off

end

