from __future__ import print_function
import os
import inspect

from pkg_resources import resource_filename

#directories
__ROOT__ = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
# default library directory
# libsdir = os.path.abspath(os.path.join(__ROOT__, '../libs/'))
libsdir = resource_filename('pyphot', 'libs')