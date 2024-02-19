# version         0.1.7
# date            04.02.2024
# author          Alexander-Maurice Illig
# affilation      Chair of Biotechnology, RWTH Aachen
# email           a.illig@biotec.rwth-aachen.de

__version__ = '0.1.7'
__author__ = 'Alexander-Maurice Illig'

from ._encoding import Encode, get_data, generate_dataframe
from ._predictors import CombinedPredictor
from ._explore import Explore
from ._utils import *