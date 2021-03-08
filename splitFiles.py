import os
import argparse
import subprocess

argparser = argparse.ArgumentParser();
argparser.add_argument('a', type = int)
argparser.add_argument('b', type = int)
argparser.add_argument

args = argparser.parse_args()

a=args.a
b=args.b

aTobPath = '{}to{}'.format(a,b-1)
if not os.path.isdir('subd/input/{}'.format(aTobPath)):
    os.makedirs('subd/input/{}'.format(aTobPath))

for i in range(a,b):
    subprocess.run(['mv', 'subd/input/{}.poly'.format(i), 'subd/input/{}/'.format(aTobPath)])
