function [t, c, m_w, elapsed_time] = reconstruct_analog(y, T, args, prev_data)

ticID = tic;

m_w = 0;  % not used
c = y(1);
t = T - y(2)/c;

if c == 0 || t < 0
    t = 0;
end

% TESTING: replace c with actual spike integral
%c = y(3);

% so that you can see spikes on plot (for test_generic_reconstruct)
%c = c*50;

elapsed_time = toc(ticID);

% WARNING: sometimes t is slightly less than 0 (now fixed)

end
