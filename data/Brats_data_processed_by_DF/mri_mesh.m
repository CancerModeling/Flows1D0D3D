fname = "Brats18_2013_13_1_t1ce";
brain = niftiread(fname + ".nii.gz");
% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Doc/Basics
% v2m args
% arg1: input data
% arg2: threshold to define outside boundary
% arg3: maximum size of surface triangle
% arg4: maximum volume of element

% [node,elem,face] = v2m(brain,0.5,100,500);
% mesh.node = node; mesh.elem = elem; mesh.face = face;
% disp(size(mesh.node))
% disp(size(mesh.elem))
% save("mesh1.mat", '-struct', 'mesh');

% [node,elem,face] = v2m(brain,0.5,50,200);
% mesh.node = node; mesh.elem = elem; mesh.face = face;
% disp(size(mesh.node))
% disp(size(mesh.elem))
% save("mesh2.mat", '-struct', 'mesh');

% [node,elem,face] = v2m(brain,0.5,20,100);
% mesh.node = node; mesh.elem = elem; mesh.face = face;
% disp(size(mesh.node))
% disp(size(mesh.elem))
% save("mesh3.mat", '-struct', 'mesh');

[node,elem,face] = v2m(brain,0.5,5,25);
mesh.node = node; mesh.elem = elem; mesh.face = face;
disp(size(mesh.node))
disp(size(mesh.elem))
save("mesh4.mat", '-struct', 'mesh');

% [node,elem,face] = v2m(brain,0.5,4,20);
% mesh.node = node; mesh.elem = elem; mesh.face = face;
% disp(size(mesh.node))
% disp(size(mesh.elem))
% save("mesh5.mat", '-struct', 'mesh');

% [node,elem,face] = v2m(brain,0.5,3,18);
% mesh.node = node; mesh.elem = elem; mesh.face = face;
% disp(size(mesh.node))
% disp(size(mesh.elem))
% save("mesh6.mat", '-struct', 'mesh');

% [node,elem,face] = v2m(brain,0.5,2.5,16);
% mesh.node = node; mesh.elem = elem; mesh.face = face;
% disp(size(mesh.node))
% disp(size(mesh.elem))
% save("mesh7.mat", '-struct', 'mesh');

% [node,elem,face] = v2m(brain,0.5,2,12);
% mesh.node = node; mesh.elem = elem; mesh.face = face;
% disp(size(mesh.node))
% disp(size(mesh.elem))
% save("mesh8.mat", '-struct', 'mesh');
