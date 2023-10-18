#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 jobs on the queue using threads.

This version uses ibrun -n $nproc -o $offset
Ridvan Orsvuran, 2021
"""

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import os
import subprocess

import threading
import queue
import sys
from itertools import chain

def message(*args, **kwargs):
    print(*args, **kwargs)
    sys.stdout.flush()

nproc = int(sys.argv[1])
nworkers = int(sys.argv[2])
jobfile = sys.argv[3]

with open(jobfile) as f:
    jobs = [line.rstrip() for line in f.readlines()]

message("Number of jobs", len(jobs))
message("Number of processors", nproc)
if nworkers > len(jobs):
    nworkers = len(jobs)
    message("Number of workers (adjusted):", nworkers)
else:
    message("Number of workers", nworkers)

logdir=jobfile+"_logs"
try:
    os.makedirs(logdir)
except OSError:
    pass

q = queue.Queue()
failed = queue.Queue()

QUIT = "__QUIT__"

def run_command(command, thread_id, job_id):
    out_file = open(logdir+"/job_"+str(job_id)+"_"+str(thread_id), "w")
    stdout = out_file
    stderr = out_file
    # stdout = sys.stdout
    # stderr = sys.stderr
    # stdout = subprocess.PIPE
    # stderr = subprocess.PIPE
    cwd = os.getcwd()
    command = list(chain(["ibrun", "-n", str(nproc), "-o", str(thread_id*nproc)], command))
    try:
        p = subprocess.Popen(command,
                             stdout=stdout,
                             stderr=stderr,
                             cwd=cwd,
                             env=os.environ)
    except OSError as e:
        message(e)

    p.wait()
    out_file.close()
    return p.returncode


class MyThread(threading.Thread):
    def __init__(self, thread_id, **kwargs):
        super(MyThread, self).__init__(**kwargs)
        self.thread_id = thread_id

    def run(self):
        while True:
            jobid, item = q.get()
            if item == QUIT:
                break
            message("Started Job {jobid} with Worker {thread_id}".format(
                jobid=jobid, thread_id=self.thread_id))
            command = item.split()
            ret = run_command(command, self.thread_id, jobid)
            if ret > 0:
                message("Failed Job {jobid} with Worker {thread_id}".format(
                    jobid=jobid, thread_id=self.thread_id))
                failed.put((jobid, item))
            else:
                message("Finished Job {jobid} with Worker {thread_id}".format(
                    jobid=jobid, thread_id=self.thread_id))
            q.task_done()

        message("Worker Done {}".format(self.thread_id))
        q.task_done()


threads = []
for i in range(nworkers):
    threads.append(MyThread(i))
    threads[-1].start()
# MyThread(0).start()

jobid = 0
for job in jobs:
    q.put((jobid, job))
    jobid += 1

for i in range(nworkers):
    q.put((jobid, QUIT))
    jobid += 1
message('All task requests sent')

# block until all tasks are done
q.join()

# gather failed jobs, if any
failed_jobs = []
while not failed.empty():
    failed_id, failed_job = failed.get()
    failed_jobs.append(failed_job)

# report failed jobs
if len(failed_jobs) > 0:
    failed_jobfile = jobfile + "_failed"
    with open(failed_jobfile, "w") as f:
        for job in failed_jobs:
            f.write("{}\n".format(job))
    print("{} jobs failed. Failed jobs written in '{}'.".format(
        len(failed_jobs), failed_jobfile))
else:
    print("No failed jobs")

message('All work completed')

