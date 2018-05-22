BASE_DIR = '/work/data/shapenet/';
OFF_DIR = [BASE_DIR 'off_cars/'];
FILE = fopen([BASE_DIR 'simplified.txt'], 'w');

files = dir([OFF_DIR '*.off']);
for i = 1: size(files)
    fprintf(FILE, sprintf('%d %s\n', i, files(i).name));
    fprintf(sprintf('%d %s\n', i, files(i).name));
end