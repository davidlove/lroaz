function PlotValueOfData( rho, actualRelative, predicted )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

assert(nnz(predicted) == nnz(predicted==1))

colors = {'b','g','r','c','m','y'};
actLine = ':';
predLine = '-';

[n,l] = size(predicted);

figure(1)
clf
hold on
for ii=1:n
    posVals = predicted(ii,:) == 1;
    plot(rho(posVals), 100*actualRelative(ii,posVals), strcat(colors{ii},predLine), ...
        'LineWidth',3)
end
yl = ylim;
yl(1) = min(0,yl(1));

for ii=1:n
    plot(rho, 100*actualRelative(ii,:), strcat(colors{ii},actLine), ...
        'LineWidth',3 )
    % Go over predicted values again to cover actual in that space
    posVals = predicted(ii,:) == 1;
    plot(rho(posVals), 100*actualRelative(ii,posVals), strcat(colors{ii},predLine), ...
        'LineWidth',3)
end

ylim(yl)

xlabel('\rho', 'FontSize',14)
ylabel('Cost Decrease (%)', 'FontSize',14)

hold off

end

