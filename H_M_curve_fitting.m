% Define the objective function with two FM terms. Add more FM terms as
% required to get a good fit. (H +- Hc) => try all combination of +- and
% chose the combination with lowest resnorm and without any discontinuity

objectiveFunc = @(x,H) (2 * x(3) / pi) * atan((H + x(1)) / x(1) * tan(pi * x(5) / 2)) + x(7) * H ...
        + (2 * x(4) / pi) * atan((H + x(2)) / x(2) * tan(pi * x(6) / 2));

% Provide the initial guess for the parameters. Change as per base material
% to get a good convergence
initialGuess = [7.1, 6.1, 0.3, 0.4, 0.01, 0.02, 2.3e-5];

% Load the experimental data (x: H values, y: corresponding experimental M measurements)
x = pure_CuO.H;   % Example H values
y = pure_CuO.M;  % Example corresponding experimental M measurements

% Perform curve fitting
[fitParams,resnorm] = lsqcurvefit(objectiveFunc, initialGuess, x, y);

% Extract the optimized parameter values
Hc1_optimized = fitParams(1);
Hc2_optimized = fitParams(2);
Ms1_optimized = fitParams(3);
Ms2_optimized = fitParams(4);
S1_optimized = fitParams(5);
S2_optimized = fitParams(6);
chi_optimized = fitParams(7);

% Display the optimized parameter values
disp('Optimized Parameter Values:');
disp(['Hc1: ', num2str(Hc1_optimized)]);
disp(['Hc2: ', num2str(Hc2_optimized)]);
disp(['Ms1: ', num2str(Ms1_optimized)]);
disp(['Ms2: ', num2str(Ms2_optimized)]);
disp(['S1: ', num2str(S1_optimized)]);
disp(['S2: ', num2str(S2_optimized)]);
disp(['chi: ', num2str(chi_optimized)]);
disp(['resnorm: ', num2str(resnorm)]);

% Separation of the desired fitted model variables
FM1 = (2 * Ms1_optimized / pi) * atan((x + Hc1_optimized) / Hc1_optimized * tan(pi * S1_optimized / 2));
FM2 = (2 * Ms2_optimized / pi) * atan((x + Hc2_optimized) / Hc2_optimized * tan(pi * S2_optimized / 2));
AFM = chi_optimized*x;
Total_M = FM1+FM2+AFM;

% Ploting not required for calculation. Modify as required.
plot(x,Total_M);
hold on;
scatter(x,y,10,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[0 0 .7],...
              'LineWidth',0.5);
hold on;