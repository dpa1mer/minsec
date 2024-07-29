function num_fully_processed = dataset_experiment(datadir)

mesh_files = dir(fullfile(datadir, "*.obj"));

num_fully_processed = 0;
for j = 1:length(mesh_files)
    the_file = mesh_files(j);
    filename = fullfile(the_file.folder, the_file.name);

    [~, pure_name, ~] = fileparts(the_file.name);
    out_name = fullfile(the_file.folder, strcat(pure_name, ".mat"));

    if exist(out_name, 'file')
        continue;
    end

    try
        [verts, faces] = readOBJ(filename);
        verts = verts(:, 1:3); % sometimes the vertices are higher-dimensional...?
        meshData = process_surface(verts, faces, NormalizeArea=1);

        % Cut along a plane if we don't already have a boundary
        if meshData.nb <= 10
            fprintf("No boundary! Cutting...\n");
            [verts, faces] = half_space_intersect(meshData.verts, meshData.faces, [0 0 0], [0 0 1], 'Cap', false);
            meshData = process_surface(verts, faces, NormalizeArea=1);
        end

        conn_cpts = connected_components(meshData.faces).';
        if max(conn_cpts) > 1
            fprintf("More than one connected component! Taking largest...\n");
            cc_size = accumarray(conn_cpts, 1, [max(conn_cpts) 1]);
            [~, biggest_ix] = max(cc_size);
            biggest_cpt = (conn_cpts == biggest_ix);
            biggest_cpt_vert_ix = nan(meshData.nv, 1);
            biggest_cpt_vert_ix(biggest_cpt) = 1:cc_size(biggest_ix);
            faces = reshape(biggest_cpt_vert_ix(meshData.faces), [], 3);
            faces = faces(all(~isnan(faces), 2), :);
            verts = meshData.verts(biggest_cpt, :);
            meshData = process_surface(verts, faces, NormalizeArea=1);
        end

        if max(connected_components(meshData.faces)) > 1
            warning("Still more than one connected component! Skipping...");
            continue;
        end

        [S, C, f, z, admmData] = minsec(meshData, 64, .1, 4, .2, UseGPU=(gpuDeviceCount>0));
        save(out_name, "meshData", "S", "C", "f", "z", "admmData");
        z_mbo = mbo_intrinsic(meshData, 0.1, 4);
        z_knoppel = mbo_intrinsic(meshData, -1, 4);
        save(out_name, 'z_mbo', 'z_knoppel', '-append');
        num_fully_processed = num_fully_processed + 1;
    catch E
        warning(E.message);
    end
end

end
