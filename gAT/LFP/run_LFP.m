% batch-run nlx_spikephase (LFP spike locking analysis)

close all;

elec_vec = [10 11 12 13 14 14 15 16 3 8 9 9];
clus_vec = [1 1 1 2 2 4 1 1 1 2 1 2];

for i = 1 : length(elec_vec)
    nlx_spikephase(elec_vec(i), clus_vec(i), i);
end