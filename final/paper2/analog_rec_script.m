% Script to run reconstruct_and_save for AT and aFRI

% 10 to 20
for i = 20:-1:10
    reconstruct_and_save('AT', i, [1]);
    reconstruct_and_save('analog', i, [1]);
end