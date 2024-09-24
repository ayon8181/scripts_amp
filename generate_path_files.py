#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generate path files for pypaw

:copyright:
   Ridvan Orsvuran (orsvuran@geoazur.unice.fr), 2017
:license:
    GNU Lesser General Public License, version 3 (LGPLv3)
    (http://www.gnu.org/licenses/lgpl-3.0.en.html)
"""

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import argparse
import os
import json
import yaml


class FileOperator(object):
    """Documentation for FileOperator

    """
    def __init__(self):
        super(FileOperator, self).__init__()

    def load_json(self, filename):
        with open(filename) as f:
            data = json.load(f)
        return data

    def dump_json(self, data, filename):
        with open(filename, "w") as file_:
            json.dump(data, file_, indent=2, sort_keys=True)
        print("{} is written.".format(filename))

    def load_yaml(self, filename):
        with open(filename) as f:
            data = yaml.safe_load(f)
        return data

    def dump_yaml(self, data, filename):
        with open(filename, "w") as f:
            data = yaml.safe_dump(data, f, indent=2)
        return data
        print("{} is written.".format(filename))

    def load_events(self, filename):
        with open(filename) as f:
            events = [line.replace("\n", "") for line in f]
        return events

    def parse_folder(self, folder, parent):
        name, data = list(folder.items())[0]
        paths = {}
        files = {}
        if isinstance(data, str):
            path = data
            subfolders = []
            child_files = []
        else:
            path = data.get("path", name)
            subfolders = data.get("subfolders", [])
            child_files = data.get("files", [])

        paths[name] = os.path.join(parent, path)

        for child_file in child_files:
            fname, filepath = list(child_file.items())[0]
            files[fname] = os.path.join(paths[name], filepath)

        for subfolder in subfolders:
            sub_paths, sub_files = self.parse_folder(subfolder,
                                                     paths[name])
            paths.update(sub_paths)
            files.update(sub_files)

        return paths, files

    def load_path_config(self, filename):
        data = self.load_yaml(filename)
        root = os.getcwd()
        folders = {}
        files = {}
        for subfolder in data:
            subfolders, subfiles = self.parse_folder(subfolder,
                                                     parent=root)
            folders.update(subfolders)
            files.update(subfiles)
        return folders, files

    def makedir(self, dirname):
        try:
            os.makedirs(dirname)
        except OSError:
            pass


class FileGenerator(FileOperator):
    """Documentation for FileGenerator

    """
    def __init__(self):
        super(FileGenerator, self).__init__()
        self.generators = []

    def f(self, name, eventname=None, period_band=None):
        return self.files[name].format(eventname=eventname,
                                       period_band=period_band)

    def d(self, name, eventname=None, period_band=None):
        return self.folders[name].format(eventname=eventname,
                                         period_band=period_band)

    def generate_list_events(self):
        for eventname in self.events:
            print(eventname)

    def generate_list_period_bands(self):
        for period_band in self.settings["period_bands"]:
            print(period_band)

    def generate_folders(self):
        for eventname in self.events:
            for period_band in self.settings["period_bands"]:
                for name, path in self.folders.items():
                    self.makedir(path.format(eventname=eventname,
                                             period_band=period_band))

    def generate_converter(self):
        for eventname in self.events:
            for seis_type in ["prem_q16"]:#["obsd_glad", "synt"]:
                data = {
                    "filetype": "sac",
                    "output_file": self.f("raw_"+seis_type, eventname),
                    "tag": self.settings["raw_{}_tag".format(seis_type)],
                    "quakeml_file": self.f("quakeml", eventname),
                    "staxml_dir": self.d("staxml", eventname),
                    "waveform_dir": self.d("{}_sac".format(seis_type),
                                           eventname)
                }
                self.dump_json(data,
                               self.f("{}_converter_path".format(seis_type),
                                      eventname))

    def generate_proc(self):
        #tags=["obsd_3D_crust","obsd_glad","obsd_1D_crust","prem_3D_crust","obsd_3D_crust_2","obsd_3D_crust_15","prem_3D_atten","1D_ref"]
        tags=["obsd_3D_crust_18","obsd_3D_crust_1"]
        for eventname in self.events:
            for period_band in self.settings["period_bands"]:
                #for j,seis_type in enumerate(["obsd_3D_crust", "synt","obsd_glad","real_data","obsd_1D_crust","prem_3D_crust"]):
                for j,seis_type in enumerate(["obsd_3D_crust_18","obsd_3D_crust_1"]):#"synt","obsd_3D_crust","obsd_glad","obsd_1D_crust","prem_3D_crust","obsd_3D_crust_2","obsd_3D_crust_15","prem_3D_atten","1D_ref","prem_q16"]):
                    data = {
                        "input_asdf": self.f("raw_"+seis_type, eventname), # "/work2/09038/ayon8181/frontera/new_events/synthetics/asdf_data/"+eventname+"."+seis_type+".h5",#self.f("raw_"+seis_type, eventname),
                        "input_tag": self.settings["raw_{}_tag".format(tags[j])],  # NOQA
                        "output_asdf": self.f("proc_"+seis_type,
                                              eventname, period_band),
                        "output_tag": self.settings["proc_{}_tag".format(tags[j])]  # NOQA
                    }
                    self.dump_json(data,
                                   self.f("{}_proc_path".format(seis_type),
                                          eventname, period_band))
    


    def generate_windows(self):
        for eventname in self.events:
            for period_band in self.settings["period_bands"]:
                data = {
                    "figure_mode": True,
                    "obsd_asdf": self.f("proc_obsd", eventname, period_band),
                    "obsd_tag": self.settings["proc_obsd_tag"],
                    "output_file": self.f("windows_file", eventname, period_band),  # NOQA
                    "synt_asdf": self.f("proc_synt", eventname, period_band),
                    "synt_tag": self.settings["proc_synt_tag"]
                }
                self.dump_json(data,
                               self.f("window_path", eventname, period_band))

    def generate_measure(self):
        for eventname in self.events:
            for period_band in self.settings["period_bands"]:
                data = {
                    "figure_dir": self.d("measure_figures"),
                    "figure_mode": True,
                    "obsd_asdf": self.f("proc_obsd", eventname, period_band),
                    "obsd_tag": self.settings["proc_obsd_tag"],
                    "output_file": self.f("measure_file", eventname, period_band),  # NOQA
                    "synt_asdf": self.f("proc_synt", eventname, period_band),
                    "synt_tag": self.settings["proc_synt_tag"],
                    "window_file": self.f("windows_file", eventname, period_band)  # NOQA
                }
                self.dump_json(data,
                               self.f("measure_path", eventname, period_band))

    def generate_stations(self):
        for eventname in self.events:
            data = {
                "input_asdf": self.f("raw_synt", eventname),
                "outputfile": self.f("stations_file", eventname),
            }
            self.dump_json(data,
                           self.f("stations_path", eventname))

    def generate_filter(self):
        for eventname in self.events:
            for period_band in self.settings["period_bands"]:
                data = {
                    "measurement_file": self.f("measure_file", eventname, period_band),  # NOQA
                    "output_file": self.f("windows_filter_file", eventname, period_band),  # NOQA
                    "station_file": self.f("stations_file", eventname),
                    "window_file": self.f("windows_file", eventname, period_band)  # NOQA
                }
                self.dump_json(data,
                               self.f("filter_path", eventname, period_band))

    def generate_adjoint(self):
        for eventname in self.events:
            for period_band in self.settings["period_bands"]:
                data = {
                    "figure_dir": self.d("adjoint_figures"),
                    "figure_mode": False,
                    "obsd_asdf": self.f("proc_obsd", eventname, period_band),
                    "obsd_tag": self.settings["proc_obsd_tag"],
                    "output_file": self.f("adjoint_file", eventname, period_band),  # NOQA
                    "synt_asdf": self.f("proc_synt", eventname, period_band),
                    "synt_tag": self.settings["proc_synt_tag"],
                    "window_file": self.f("windows_filter_file", eventname, period_band)  # NOQA
                }
                self.dump_json(data,
                               self.f("adjoint_path", eventname, period_band))

    def count_windows(self, eventname):
        measurements = {}
        for period_band in self.settings["period_bands"]:
            windows_file = self.f("windows_filter_file", eventname, period_band)  # NOQA
            windows_data = self.load_json(windows_file)
            measurements[period_band] = dict((x, 0) for x in ["BHR", "BHT", "BHZ"])  # NOQA
            for station, components in windows_data.items():
                for component, windows in components.items():
                    c = component.split(".")[-1]
                    measurements[period_band][c] += len(windows)
        return measurements

    def get_ratio(self, eventname):
        counts = self.count_windows(eventname)
        for p in counts:
            for c in counts[p]:
                counts[p][c] = 1 / counts[p][c]
        return counts

    def generate_weight_params(self):
        template = self.load_yaml(self.f("weight_template"))
        for eventname in self.events:
            ratio = self.get_ratio(eventname)
            data = template.copy()
            data["category_weighting"]["ratio"] = ratio
            self.dump_yaml(data, self.f("weight_parfile", eventname))

    def generate_weight_paths(self):
        for eventname in self.events:
            data = {"input": {},
                    "logfile": self.f("weight_log", eventname)}
            for period_band in self.settings["period_bands"]:
                data["input"][period_band] = {
                    "asdf_file": self.f("proc_synt", eventname, period_band),
                    "output_file": self.f("weight_file", eventname, period_band),  # NOQA
                    "station_file": self.f("stations_file", eventname),
                    "window_file": self.f("windows_filter_file", eventname, period_band)  # NOQA
                }
            self.dump_json(data, self.f("weight_path", eventname))

    def generate_sum(self):
        for eventname in self.events:
            data = {"input_file": {},
                    "output_file": self.f("sum_adjoint_file", eventname)}
            for period_band in self.settings["period_bands"]:
                data["input_file"][period_band] = {
                    "asdf_file": self.f("adjoint_file", eventname, period_band),  # NOQA
                    "weight_file": self.f("weight_file", eventname, period_band),  # NOQA
                }
            self.dump_json(data, self.f("sum_adjoint_path", eventname))

    def run(self):
        steps = dict([(x[9:], x) for x in dir(self)
                      if x.startswith('generate_')])
        parser = argparse.ArgumentParser(
            description='Generate path files')
        parser.add_argument('-p', '--paths-file', default="paths.yml")
        parser.add_argument('-s', '--settings-file', default="settings.yml")
        parser.add_argument('-e', '--events-file', default="event_list")
        parser.add_argument("step",  help="file to generate",
                            choices=steps.keys())
        args = parser.parse_args()
        self.folders, self.files = self.load_path_config(args.paths_file)
        self.settings = self.load_yaml(args.settings_file)
        self.events = self.load_events(args.events_file)

        getattr(self, steps[args.step])()


if __name__ == '__main__':
    FileGenerator().run()
