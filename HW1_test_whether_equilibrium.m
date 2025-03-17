% Define the range of years
x = 1985:2005;

% Calculate p_CO2, pH, and actual [CO3^2-]
p_CO2 = 1.887 * x - 3427;
pH = -0.0018 * x + 11.7;
CO3_actu = -0.4867 * x + 1206;

% Calculate theoretical [CO3^2-]
K1 = 10^(-6.35);
K2 = 10^(-10.33);
K_CO2 = 10^(-1.46);
CO3_theo = (10^(-18.14) * (1.887 * x - 3427)) ./ 10.^(-2 * (-0.0018 * x + 11.7));

% Display the results in the command window
fprintf('Year\tActual [CO3^2-]\tTheoretical [CO3^2-]\n');
for i = 1:length(x)
    fprintf('%d\t%.2f\t\t%.2f\n', x(i), CO3_actu(i), CO3_theo(i));
end

% Plot the results
figure;
plot(x, CO3_actu, 'Color', [53/255 61/255 38/255], 'LineWidth', 2); % Actual [CO3^2-]
hold on;
plot(x, CO3_theo, 'Color', [189/255 203/255 177/255], 'LineWidth', 2); % Theoretical [CO3^2-]
xlabel('Year');
ylabel('[CO_3^{2-}] (\mu mol/kg)');
title('Comparison of Actual and Theoretical [CO_3^{2-}]', 'FontSize',20);
legend('Actual [CO_3^{2-}]', 'Theoretical [CO_3^{2-}]');
grid on;
hold off;