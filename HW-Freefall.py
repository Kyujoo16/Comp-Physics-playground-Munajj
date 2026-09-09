import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("g", help="user chooses gravity constant ", type=float)
parser.add_argument("h", help="user chooses intial height", type=float)
args = parser.parse_args()
t= np.sqrt(2*(args.h/args.g))
if args.g < 0:
    print("Error gravity must be greater than zero value")
if args.h < 0:
    print("Error h must be greater or equal to zero")
else:
    print("time that ball hits ground =", t)
