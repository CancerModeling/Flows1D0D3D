import os
import flows1d0d3d as f


def run():
    data_folder = '../../data'
    output_folder = './tmp'

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    tau = 1. / 2**16
    tau_out = 1. / 2**6
    t_end = 10.
    t = 0

    s = f.HeartToBreast1DSolver()
    s.set_output_folder('./tmp')
    #s.set_path_inflow_pressures(os.path.join(data_folder, ""));
    s.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"));
    s.set_path_linear_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"));
    s.set_path_coupling_conditions(os.path.join(data_folder, "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"));
    s.setup(degree, tau)

    for i in range(int(t_end / tau_out)):
        print ('iter = {}, t = {}'.format(i, t))
        t = s.solve_flow(tau, t, int(tau_out / tau))
        s.write_output(t)


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()
