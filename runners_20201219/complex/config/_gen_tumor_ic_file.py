
def gen_tumor_ic_file(file_dir, dp):

    # write to file
    inpf = open(file_dir + dp.tumor_ic_file, 'w')
    inpf.write("type, cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx, hyp_ry, hyp_rz\n")

    # type = 1 -- spherical tumor core
    # type = 2 -- elliptical tumor core (sharp)
    # type = 3 -- spherical tumor core and then spherical hypoxic core
    # type = 5 -- spherical tumor core (sharp)
    ic_type = dp.tumor_ic_type
    center = dp.tumor_ic_center
    r = dp.tumor_ic_radius
    for i in range(len(center)):
        inpf.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(ic_type[i], center[i][0], center[i][1], center[i][2], r[i][0], r[i][1], r[i][2], r[i][0], r[i][1], r[i][2]))

    inpf.close()
