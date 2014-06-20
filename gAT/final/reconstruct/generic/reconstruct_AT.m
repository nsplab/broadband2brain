function [t, c, m_w, elapsed_time] = reconstruct_AT(y, T, args, prev_data)

ticID = tic;

m_w = 0;  % not used
c = y(1);
t = T/2;

%if c == 0
%    t = 0;
%end

% so that you can see spikes on plot (for test_generic_reconstruct)
%c = c*50;

elapsed_time = toc(ticID);

end
