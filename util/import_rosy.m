function field = import_rosy(filename)

f = fopen(filename, "r");
nf = fscanf(f, "%d\n", 1);
four = fscanf(f, "%d\n", 1);
assert(four == 4);
field = textscan(f, "%f %f %f\n", nf, CollectOutput=true);
field = field{1};
fclose(f);

end