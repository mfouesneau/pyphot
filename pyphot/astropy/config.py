from __future__ import print_function
import os
import inspect


#directories
__ROOT__ = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
# default library directory
# libsdir = os.path.abspath(os.path.join(__ROOT__, '../libs/'))
# from pkg_resources import resource_filename
# libsdir = resource_filename('pyphot', 'libs')
from importlib import resources
libdir = os.path.join(resources.files('pyphot'), 'libs')