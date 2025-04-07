eps_range = logspace(log10(1/200), log10(1/100), 10);
fit_x = 1./eps_range;
mean_times_5 = [2.1759e4 1.8192e4 1.5965e4 1.3233e4 1.1635e4 1.0121e4 8.8772e3 7.5472e3 6.5578e3 5.9606e3];
mean_times_1 = [5.8185e4 4.0357e4 2.8903e4 2.1261e4 1.5688e4 1.2242e4 9.6482e3 7.3557e3 5.7020e3 4.8464e3];
fit_y_5 = log(mean_times_5);
fit_y_1 = log(mean_times_1);

p5 = polyfit(fit_x(1:10), fit_y_5(1:10), 1);
p1 = polyfit(fit_x(1:10), fit_y_1(1:10), 1);
p5(1)


plot(fit_x, fit_y_1, 'kx', fit_x, p1(1)*fit_x+p1(2), 'k-', fit_x, fit_y_5, 'bx', fit_x, p5(1)*fit_x+p5(2), 'b-', "LineWidth", 2, "MarkerSize", 10)
xlabel('1/\epsilon'); ylabel('log(\tau)');
legend('\alpha = 1', 'Slope= 0.0249', '\alpha = 5', 'Slope= 0.0130', "Location", 'northwest')
