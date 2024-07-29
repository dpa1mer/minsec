function [verts, faces] = quad_mesh_quadwild_bimdf(quadwild_dir, working_dir, mesh_data, z)

% get absolute paths because we're going to cd
assert(isfolder(working_dir));
[~, info] = fileattrib(working_dir);
working_dir = info.Name;

assert(isfolder(quadwild_dir));
[~, info] = fileattrib(quadwild_dir);
quadwild_dir = info.Name;

trimesh_file = fullfile(working_dir, "trimesh.obj");
quadwild_patch_file = fullfile(working_dir, "trimesh_rem_p0.obj");
rosy_file = fullfile(working_dir, "field.rosy");
quadmesh_file = fullfile(working_dir, "trimesh_rem_p0_1_quadrangulation.obj");

% write data to files
writeOBJ(trimesh_file, mesh_data.verts, mesh_data.faces);
export_rosy(rosy_file, mesh_data, z);

% run meshing code
prevdir = pwd;
cd(quadwild_dir);
system(sprintf("%s %s %d %s %s", "./quadwild", trimesh_file, 3, "config/prep_config/basic_setup_larger.txt", rosy_file));
system(sprintf("%s %s %d %s", "./quad_from_patches", quadwild_patch_file, 1, "config/main_config/flow_larger.txt"));
cd(prevdir);

% load resulting quad mesh
[verts, faces] = readOBJ(quadmesh_file, Quads=true);

end