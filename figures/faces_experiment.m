function faces_experiment(outdir)

for j = 1:5
    face = random_face;

    % Degree 4
    [S, C, f, z, admmData] = minsec(face, 64, .1, 4, .1);
    save(fullfile(outdir, sprintf("face%d_deg4.mat", j)), "face", "S", "C", "f", "z", "admmData");

    % Degree 2
    [S, C, f, z, admmData] = minsec(face, 64, .1, 2, .1);
    save(fullfile(outdir, sprintf("face%d_deg2.mat", j)), "face", "S", "C", "f", "z", "admmData");
end

end