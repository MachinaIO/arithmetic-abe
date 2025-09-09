#!/usr/bin/env sage -python
# import sys
from estimator.estimator import *
from estimator.estimator.lwe_parameters import *
from estimator.estimator.nd import *
import math
import datetime
import os
from decimal import Decimal, getcontext
from norms import CircuitNorms
import os
import subprocess
import json

getcontext().prec = 300
script_dir = os.path.dirname(os.path.abspath(__file__))
