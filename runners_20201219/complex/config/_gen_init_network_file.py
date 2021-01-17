
def gen_init_network_file(file_dir, dp):

    L = dp.L
    R = dp.vessel_radius
    pressures = dp.vessel_pressures
    l1 = dp.vessel_line_1
    l2 = dp.vessel_line_2

    # write to file
    inpf = open(file_dir + dp.network_init_file, 'w')
    inpf.write("DGF\n")

    # vertex
    inpf.write("Vertex\n")
    inpf.write("parameters {}\n".format(1))
    inpf.write("{} {} {} {}\n".format(l1[0], l1[1], l1[2], pressures[0]))
    inpf.write("{} {} {} {}\n".format(l1[3], l1[4], l1[5], pressures[1]))
    inpf.write("{} {} {} {}\n".format(l2[0], l2[1], l2[2], pressures[2]))
    inpf.write("{} {} {} {}\n".format(l2[3], l2[4], l2[5], pressures[3]))

    # segments
    inpf.write("#\n")
    inpf.write("SIMPLEX\n")
    inpf.write("parameters {}\n".format(2))
    inpf.write("{} {} {} {}\n".format(0, 1, R[0], 0.0075))
    inpf.write("{} {} {} {}\n".format(2, 3, R[1], 0.0075))
    # 
    inpf.write("#\n")
    inpf.write("BOUNDARYDOMAIN\n")
    inpf.write("default {}\n".format(1))

    # 
    inpf.write("#")
    inpf.close()

