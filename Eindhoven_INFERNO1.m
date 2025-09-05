clear; % Clear all variables
clc; % Clear window

%CHEB compute D=differentiation matrix, x=Chebyshev grid
function [D, x] = cheb(N)
  if N == 0
    D = 0; x = 1; return
  end
  x = cos(pi * (0:N)' / N);
  c = [2; ones(N-1,1); 2] .* (-1).^(0:N)';
  X = repmat(x, 1, N+1);
  dX = X - X';
  D = (c * (1./c)') ./ (dX + eye(N+1));
  D = D - diag(sum(D, 2));
end

N = 20;  % Number of Grid points, fewer points for higher acuracy
[D, Y] = cheb(N);  % Y from 1 (top hot) to -1 (bottom cold)
D2 = D^2;  % Second derivative matrix in Y
D2_y = 4 * D2;  % Scaled for y* = (Y + 1)/2, d^2/dy*^2 = 4 d^2/dY^2
I = eye(N+1);  % Identity Matrix

u_star = 1.5 * (1 - Y.^2);  % Non-dim Poiseuille velocity u*/U, avg=1, max=1.5 at center
inv_u = 1 ./ u_star;  % 1/u* for equation

y_sensor = 0.5; % Sensor location (normalized height)

Pe_values = [0.1, 1, 10, 100];  % Peclet numbers to vary
colors = {[1 0.5 0], [0 0 1], [0 1 0], [1 1 0]};  % Orange, blue, green, yellow
total_xstar = 0.05;  % Non-dim sensor position x* = x/H; tuned for shapes
num_steps = 2000;  % Marching steps: high for accuracy in implicit scheme
dx_star = total_xstar / num_steps;  % Step size in x*; implicit stable for any size
x_plot = (1 - Y) / 2;  % Normalized position along sensor: 0 at top hot to 1 at bottom cold

figure; % Initializing the graph
hold on;
grid on;

for ip = 1:length(Pe_values)
  Pe = Pe_values(ip);  % Current Pe
  T = zeros(N+1, 1);  % Initial T* at inlet: uniform cold=0 interior
  T(1) = 1;  % Fixed top hot T*=1 (Y=1)
  T(end) = 0;  % Fixed bottom cold T*=0 (Y=-1)

  for n = 1:num_steps
    % Implicit matrix A = I - (dx*/Pe) * diag(1/u*) * D2_y
    A = I - (dx_star / Pe) * diag(inv_u) * D2_y;
    % Enforce Dirichlet BCs
    A(1, :) = 0; A(1, 1) = 1;
    A(end, :) = 0; A(end, end) = 1;
    % RHS = current T, with BCs
    R = T;
    R(1) = 1;
    R(end) = 0;
    % Solve for next T
    T = A \ R;
  end

  % Plot T* vs x_plot
  plot(x_plot, T, 'Color', colors{ip}, 'LineWidth', 2);
end

% Add sensor line
xline(y_sensor, '--k', 'Sensor Line', 'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5);

ax = gca;
set(ax, 'Color', 'white'); 
ax.GridColor = 'b';
ax.XColor = 'w';
ax.YColor = 'w';

title('GRAPH OF TEMPERATURE AGAINST THE LENGTH');
xlabel('x (sensor line at y=0.5)');
ylabel('Temperature');
legend(arrayfun(@(p) sprintf('Pe=%.1f', p), Pe_values, 'uni', 0), 'Location', 'northeast');
hold off;