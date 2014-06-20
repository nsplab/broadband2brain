% Reconstructs Dr. Andrew Richardson's figure
% Nemo channel 5
function [] = tuning_curve_plot()

close all;

method_name = 'gibbs'; %'ker';
paramIDs = 34; %2;
data_sets = [2, 5:17, 19];  % Sync with base_ind etc below

[dir_real dir] = est_delay_b(method_name, paramIDs, data_sets, 1, [5]);
[dir_real_ff dir_ff] = est_delay_b(method_name, paramIDs, data_sets, 2, [5]);
[dir_real_wash dir_wash] = est_delay_b(method_name, paramIDs, data_sets, 3, [5]);

dir_real
dir_real_ff
dir_real_ff - dir_real

x_real = fix(dir_real_ff - dir_real);
y_real = fix(dir_real_wash - dir_real_ff);
x = fix(dir_ff - dir);
y = fix(dir_wash - dir_ff);

base_ind = [2, 3, 8:10, 14, 15]; %[5, 6, 11:13, 17, 19];
cw_ind = 4:7; %7:10;
ccw_ind = [1, 11:13]; %[2, 14:16];

figure;
hold on;
plot([x_real ; x], [y_real ; y], 'k:');
p1 = scatter(x_real(base_ind), y_real(base_ind), 'bo');
p2 = scatter(x_real(cw_ind), y_real(cw_ind), 'ro');
p3 = scatter(x_real(ccw_ind), y_real(ccw_ind), 'go');
p4 = scatter(x(base_ind), y(base_ind), 'b*');
p5 = scatter(x(cw_ind), y(cw_ind), 'r*');
p6 = scatter(x(ccw_ind), y(ccw_ind), 'g*');
legend([p1 p2 p3 p4 p5 p6], 'real (base)', 'real (cw)', 'real (ccw)', 'FRI (base)', 'FRI (cw)', 'FRI (ccw)');
xlabel('force - base');
ylabel('wash - force');
title('preferred direction shifts');
axis equal;
xlim([-3 3]);
ylim([-3 3]);

end

function v = fix(v)
for i = 1 : length(v)
    if v(i) > pi
        v(i) = v(i) - 2*pi;
    elseif v(i) < -pi
        v(i) = v(i) + 2*pi
    end
end
end