import os
import numpy as np
import pandas as pd
import csv
import sys
import matplotlib.pyplot as plt
import importlib

import logging

from read_write_dgf import read_dgf 
from read_write_vtk import write_vtk


def convert_dgf_to_vtk(in_file, out_file):

    # read dgf file
    nodes, nodes_data, segments, segments_data = read_dgf(in_file)

    print("num nodes = {}".format(len(nodes)))
    # print(nodes)
    # print(nodes_data)

    print("num segments = {}".format(len(segments)))
    # print(segments)
    # print(segments_data)

    # write
    write_vtk(out_file, nodes, nodes_data, segments, segments_data)



convert_dgf_to_vtk("test.dgf", "test.vtk")
