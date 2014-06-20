function [true_positives, false_positives, false_negatives, error_rate, tp_rate, fp_rate, fn_rate] = compare_spikes(spk, real_spk, tolerance)
% Compares 2 spike trains on a spike by spike basis

% Inputs:
% spk - vector of experimental spike times
% real_spk - vector of real spike times to be compared against
% tolerance - a spike in spk is considered to match a spike in real_spk if it differs in time by no more than tolerance

% Outputs:
% true_positives - indices (in spk) of spikes identified correctly
% false_positives - indices (in spk) of extra spikes detected that are not in real_spk
% false_negatives - indices (in real_spk) of spikes missed
% error_rate - errors per real spike

spk = sort(spk);
real_spk = sort(real_spk);

spk_ind = 1;
true_positives = zeros(min(length(spk), length(real_spk)), 1);
tp_ind = 1;
false_positives = zeros(length(spk), 1);
fp_ind = 1;
false_negatives = zeros(length(real_spk), 1);
fn_ind = 1;

for real_ind = 1 : length(real_spk)
    while 1
        if spk_ind <= length(spk) && abs(spk(spk_ind) - real_spk(real_ind)) <= tolerance
            % True positive
            true_positives(tp_ind) = spk_ind;
            tp_ind = tp_ind + 1;
            spk_ind = spk_ind + 1;
            break;
        elseif spk_ind <= length(spk) && spk(spk_ind) < real_spk(real_ind)
            % False positive
            false_positives(fp_ind) = spk_ind;
            fp_ind = fp_ind + 1;
            spk_ind = spk_ind + 1;
        else
            % False negative (missed spike)
            false_negatives(fn_ind) = real_ind;
            fn_ind = fn_ind + 1;
            break;
        end
    end
end
while spk_ind <= length(spk)
    % False positive
    false_positives(fp_ind) = spk_ind;
    fp_ind = fp_ind + 1;
    spk_ind = spk_ind + 1;
end

true_positives = true_positives(1 : tp_ind-1);
false_positives = false_positives(1 : fp_ind-1);
false_negatives = false_negatives(1 : fn_ind-1);
error_rate = (length(false_positives) + length(false_negatives)) / length(real_spk);  % Errors per real spike
tp_rate = length(true_positives) / length(real_spk);
fp_rate = length(false_positives) / length(real_spk);
fn_rate = length(false_negatives) / length(real_spk);

end
