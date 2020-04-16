#!/usr/bin/python

import os
from subprocess import check_call as call

call("python pytheas%spytheas.py" % os.sep, shell=True)