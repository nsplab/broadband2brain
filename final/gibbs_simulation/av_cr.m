% Run run_cramer_rao a bunch of times (for a bunch of different random
% signals)
% Average results
% numsigs = num random signals to try
% all other params are for run_cramer_rao

function [devc_cr devc_gibbs devc_GA devc_ml devt_cr devt_gibbs devt_GA devt_ml devc_cr_psinc devt_cr_psinc devc_cad devt_cad devc_ml_psinc devt_ml_psinc] = av_cramer_rao(numsigs, numtrials, T, K, N, sigmah, sigmae, numiter, mean_iter, type);

devc_cr_vec = zeros(1, numsigs);
devc_gibbs_vec = zeros(1, numsigs);
devc_GA_vec = zeros(1, numsigs);
devc_ml_vec = zeros(1, numsigs);
devt_cr_vec = zeros(1, numsigs);
devt_gibbs_vec = zeros(1, numsigs);
devt_GA_vec = zeros(1, numsigs);
devt_ml_vec = zeros(1, numsigs);
devc_cr_psinc_vec = zeros(1, numsigs);
devt_cr_psinc_vec = zeros(1, numsigs);
devc_cad_vec = zeros(1, numsigs);
devt_cad_vec = zeros(1, numsigs);
devc_ml_psinc_vec = zeros(1, numsigs);
devt_ml_psinc_vec = zeros(1, numsigs);

for i = 1 : numsigs   % loop over random signals
    [devc_cr devc_gibbs devc_GA devc_ml devt_cr devt_gibbs devt_GA devt_ml devc_cr_psinc devt_cr_psinc devc_cad devt_cad devc_ml_psinc devt_ml_psinc] = run_cramer_rao(i, numsigs, numtrials, T, K, N, sigmah, sigmae, numiter, mean_iter, type);
    devc_cr_vec(i) = devc_cr;
    devc_gibbs_vec(i) = devc_gibbs;
    devc_GA_vec(i) = devc_GA;
    devc_ml_vec(i) = devc_ml;
    devt_cr_vec(i) = devt_cr;
    devt_gibbs_vec(i) = devt_gibbs;
    devt_GA_vec(i) = devt_GA;
    devt_ml_vec(i) = devt_ml;
    devc_cr_psinc_vec(i) = devc_cr_psinc;
    devt_cr_psinc_vec(i) = devt_cr_psinc;
    devc_cad_vec(i) = devc_cad;
    devt_cad_vec(i) = devt_cad;
    devc_ml_psinc_vec(i) = devc_ml_psinc;
    devt_ml_psinc_vec(i) = devt_ml_psinc;
end

devc_cr = mean(devc_cr_vec);
devc_gibbs = mean(devc_gibbs_vec);
devc_GA = mean(devc_GA_vec);
devc_ml = mean(devc_ml_vec);
devt_cr = mean(devt_cr_vec);
devt_gibbs = mean(devt_gibbs_vec);
devt_GA = mean(devt_GA_vec);
devt_ml = mean(devt_ml_vec);
devc_cr_psinc = mean(devc_cr_psinc_vec);
devt_cr_psinc = mean(devt_cr_psinc_vec);
devc_cad = mean(devc_cad_vec);
devt_cad = mean(devt_cad_vec);
devc_ml_psinc = mean(devc_ml_psinc_vec);
devt_ml_psinc = mean(devt_ml_psinc_vec);

%devt_cr_psinc_vec
%devt_cad_vec
%devt_ml_psinc_vec

end