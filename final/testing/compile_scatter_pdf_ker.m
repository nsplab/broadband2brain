% Outputs terminal command to combine pdfs
% Command to be run from ~/Desktop

sr = [500, 1000, 1500];

for sampling_rate = sr

    str = 'pdftk ';

    for data_set = [1, 2]

        for elec = channels_to_use(data_set)

            str = [str 'opt_ker/opt_ker_' int2str(data_set) '_' int2str(elec) '_' int2str(sampling_rate) '.pdf '];

        end

    end

    for data_set = [1, 2]

        str = [str 'scatter_ker_' int2str(sampling_rate) '_' int2str(data_set) '.pdf '];

    end

    str = [str 'cat output ker_params_' int2str(sampling_rate) '.pdf']

end
