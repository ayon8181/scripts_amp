#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Drop in replacement for pypaw-process_asdf
"""

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import mpi
from pypaw.utils import read_yaml_file, read_json_file
from tqdm import tqdm

import pyasdf
import argparse
import os

from pytomo3d.signal.process import process_stream
from pypaw.process import update_param


def run_with_mpi_multi(func, objects, title=None):
    """MPI helper

    :param func: function to call
    :type func: function
    :param objects: objects to pass to function
    :type objects: list
    :returns: results
    :rtype: dict
    """
    if mpi.rank == 0:  # server
        jobs = mpi.split(objects, mpi.size)
    else:
        jobs = []
    jobs = mpi.comm.scatter(jobs, root=0)

    results = {}
    if mpi.rank == 0:
        pbar = tqdm(total=len(jobs), desc=title)
    for _i, job in enumerate(jobs, 1):
        results.update(func(job))
        if mpi.rank == 0:
            pbar.update()

    return results


def process_station(sta):
    try:
        st = input_ds.waveforms[sta][input_tag]
        inv = input_ds.waveforms[sta].StationXML
    except:
        return {sta: (None, None)}
    try:
        process_stream(st, inv, **params)
    except Exception as e:
           print("ERROR:", st[0].id, e)
           st = None
           inv = None
           
    return {sta: (st, inv)}


def write_output(output_ds, results):
    for st, inv in tqdm(results.values(), "Writing"):
        if st is not None:
            output_ds.add_waveforms(st, tag=output_tag,
                                    event_id=event)
            output_ds.add_stationxml(inv)

parser = argparse.ArgumentParser()
parser.add_argument('-p', action='store', dest='params_file',
                    required=True, help="parameter file")
parser.add_argument('-f', action='store', dest='path_file',
                    required=True,
                    help="path file")
parser.add_argument('-v', action='store_true', dest='verbose',
                    help="verbose flag")
args = parser.parse_args()

params = read_yaml_file(args.params_file)
path = read_json_file(args.path_file)

input_asdf = path["input_asdf"]
input_tag = path["input_tag"]
output_asdf = path["output_asdf"]
output_tag = path["output_tag"]

if mpi.is_main_rank:
    if os.path.exists(output_asdf):
        print("OUTPUT EXISTS deleting...")
        os.remove(output_asdf)

mpi.comm.barrier()

input_ds = pyasdf.ASDFDataSet(input_asdf, compression=None, mode="r")
event = input_ds.events[0]
update_param(event, params)

jobs = [sta for sta in input_ds.waveforms.list()]
results = run_with_mpi_multi(process_station, jobs, "Process")

if mpi.is_main_rank:
    output_ds = pyasdf.ASDFDataSet(output_asdf, mpi=False, mode="a")
    output_ds.events = input_ds.events
    write_output(output_ds, results)
    del output_ds
    #del input_ds
    # pbar = tqdm(total=len(results), desc="Writing")
    for i in range(1, mpi.size):
        mpi.comm.send(True, dest=i)
        mpi.comm.recv(source=i)

else:
    mpi.comm.recv(source=0)
    if len(results) > 0:
        output_ds = pyasdf.ASDFDataSet(output_asdf, mpi=False, mode="a")
        write_output(output_ds, results)
        del output_ds
        #del input_ds
    mpi.comm.send(True, dest=0)

#del input_ds
