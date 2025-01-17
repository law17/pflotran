import sys
import os
try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
  sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft


filenames = []
filenames.append('orig-005.tec')
filenames.append('dev-005.tec')

f = plt.figure(figsize=(6,6))
plt.subplot(1,1,1)
f.suptitle("1D Calcite at 25 Years (with Flow)",fontsize=16)
plt.xlabel('X [m]')
plt.ylabel('Concentration [M]')

linestyle = []
linestyle.append('-')
linestyle.append('-.')

#plt.xlim(0.,1.)
#plt.ylim(0.,1.)
#plt.grid(True)
#plt.yscale('log')

for ifile in range(len(filenames)):
  columns = [7,8,9]
#  columns = [12]
  for icol in range(len(columns)):
    data = pft.Dataset(filenames[ifile],3,columns[icol])
    plt.plot(data.get_array('x'),data.get_array('y'),
             label=data.get_name('yname'),ls=linestyle[ifile])

#'best'         : 0, (only implemented for axis legends)
#'upper right'  : 1,
#'upper left'   : 2,
#'lower left'   : 3,
#'lower right'  : 4,
#'right'        : 5,
#'center left'  : 6,
#'center right' : 7,
#'lower center' : 8,
#'upper center' : 9,
#'center'       : 10,
plt.legend(loc=2)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#      plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#        plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

plt.twinx()
plt.ylabel('Volume Fraction [-]')
plt.ylim(0.,1.1e-5)

for ifile in range(len(filenames)):
  data = pft.Dataset(filenames[ifile],1,10)
  plt.plot(data.get_array('x'),data.get_array('y'),ls=linestyle[ifile], \
           label=data.get_name('yname'))

major_formatter = plt.FormatStrFormatter('%1.0e')
plt.gca().yaxis.set_major_formatter(major_formatter)

plt.legend(loc=1)
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#      plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#        plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

f.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.14,right=.82)

plt.show()
