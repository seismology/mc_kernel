#!/usr/bin/env python
# This Python file uses the following encoding: utf-8
#-------------------------------------------------------------------
# Filename: monitor_memory.py
#  Purpose: Monitor memory usage of a Kerner slave process
#   Author: Martin van Driel & Simon Stähler
#    Email: vandriel@erdw.ethz.ch
#
# Copyright (C) 2014 Martin van Driel and Simon Stähler
#        GNU Lesser General Public License, Version 3
#            (http://www.gnu.org/copyleft/lesser.html)
#---------------------------------------------------------------------

from subprocess import Popen, PIPE, STDOUT
import time
import numpy as np
import sys
import argparse


parser = argparse.ArgumentParser(description='Monitor memory usage of Kerner.')

parser.add_argument('-q', action="store", dest="queue", default="local", 
                    choices=('local', 'lsf'), help="local (default) or queue")

results = parser.parse_args()
queue = results.queue

if queue == 'local':
    print 'checking local run'
    #cmd_mem = "ps -e -orss=,args= |grep kerner| grep -v grep|  awk '{print $1}'"
    #cmd_rtime = "ps -e |grep kerner| grep -v grep| awk '{print $3}'"
    # The tail means that the process with the highest PID will be examined. Should not
    # matter, which slave is examined, as long as it is not the master process.
    cmd_mem = "ps -eo pid,fname,rss |grep kerner| grep -v grep| awk '{print $3}' | tail -n 1"
    cmd_rtime = "ps -e |grep kerner| grep -v grep| awk '{print $3}' | tail -n 1"
elif queue == 'lsf':
    print 'checking lsf run'
    cmd_mem = "bbjobs | grep 'Total Memory' | awk '{print $4}'"
    cmd_rtime = "bbjobs | grep 'Wall-clock' | awk '{print $3 " " $4}'"

print 'recording memory consumption for kerner'
print 'will stop as soon as kerner finishes'
print 'output written to ./memory_stats.dat'

f = open('./memory_stats.dat', 'w')

meml = []

while True:
    p = Popen(cmd_mem, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    mem = p.stdout.read()
    p = Popen(cmd_rtime, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    if queue == 'local':
        rtime = p.stdout.read().strip()
    elif queue == 'lsf':
        rtime = p.stdout.read().rstrip()
    
    if mem == '':
        break
    else:
        if queue == 'local':
            mem = int(mem) / 1024.**2
            rtime = time.strptime('1970:' + rtime, "%Y:%H:%M:%S")
            print >> f, '% 8.1f % 8.3f MB' % (time.mktime(rtime) + 3600., mem)
        elif queue == 'lsf':
            mem = int(mem) / 1024.
            print >> f, '%s % 8.3f' % (rtime, mem)

        meml.append(mem)

    f.flush()

    if queue == 'local':
        time.sleep(1)
    elif queue == 'lsf':
        time.sleep(10)

f.close()

try:
    mem = np.array(meml)
    max = mem.max()
    print 'Maximum of Memory usage: ', max, ' MB'
except ValueError:
    print 'Problem: Kerner seems not to be running!'
    
