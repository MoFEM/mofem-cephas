// Gmsh project created on Tue Jan 17 23:22:23 2023
SetFactory("OpenCASCADE");
lx = 20;
ly = 5;
Rectangle(1) = {0, 0, 0, lx, ly, 0};
Physical Curve("ESSENTIAL", 1) = {4};
Physical Surface("INITIAL", 2) = {1};
