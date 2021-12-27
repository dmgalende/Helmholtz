function plotOnSameFig(fig1, fig2, leg)

f = openfig(fig1);
H = findobj(f, 'Type', 'line');
x_data = get(H, 'xdata');
y_data = get(H, 'ydata');

g = openfig(fig2);
G = findobj(g, 'Type', 'line');
x1_data = get(G, 'xdata');
y1_data = get(G, 'ydata');

figure
plot(x_data, y_data, 'k-o')
hold on
plot(x1_data, y1_data, 'b-*')
if nargin == 3
    legend(leg)
else
    legend({fig1, fig2})
end