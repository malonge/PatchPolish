#!/usr/bin/env python

import subprocess
import time
import sys


def log(message):
    """ Log messages to standard error. """
    sys.stderr.write(time.ctime() + ' --- ' + message + "\n")
    sys.stderr.flush()


def run_oe(cmd, out, err):
    """ Run a command and redirect stdout/stderr. """

    if not isinstance(out, str) or not isinstance(err, str):
        raise TypeError("out/err should be file names (strings)")

    f_out = open(out, "w")
    f_err = open(err, "w")

    log('Running: %s > %s 2> %s' % (" ".join(cmd), out, err))
    if subprocess.call(cmd, stdout=f_out, stderr=f_err) != 0:
        raise RuntimeError('Failed : %s > %s 2> %s' % (" ".join(cmd), out, err))

    log('Finished running : %s > %s 2> %s' % (" ".join(cmd), out, err))

    f_out.close()
    f_err.close()
