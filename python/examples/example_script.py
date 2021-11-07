import os
import flows1d0d3d as f

def run():
    data_folder = '../../data'
    output_folder = './tmp'

    os.mkdir(output_folder)

    degree = 2
    tau = 1e-6
    t = 0

    s = f.HeartToBreast1DSolver()
    s.set_output_folder('./tmp')
    #s.set_path_inflow_pressures(os.path.join(data_folder, ""));
    s.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"));
    s.set_path_linear_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"));
    s.set_path_coupling_conditions(os.path.join(data_folder, "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"));
    s.setup(degree, tau)

    for i in range(1000):
        s.solve_flow(tau, t)
        t += tau
        if i % 100 == 0:
            print ('iter = {}, t = {}'.format(i, t))
            s.write_output(t)

if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()
