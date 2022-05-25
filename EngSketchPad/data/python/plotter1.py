###################################################################
#                                                                 #
# plotter1.py -- sample line plotting                             #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

from pyEGADS import egads
from pyOCSM  import ocsm
from pyOCSM  import esp

import math
import time

esp.TimLoad("plotter", esp.GetEsp("pyscript"), "")

# one line
esp.TimMesg("plotter", "new|Example of one line|x-axis label|y-axis label|")
esp.TimMesg("plotter", "add|0;1;2;3;4|0;1;4;9;16|b-s|")
esp.TimMesg("plotter", "show")
esp.TimHold("pyscript", "plotter")

# two lines
esp.TimMesg("plotter", "new|Example of two lines|x-axis label|y-axis label|")
esp.TimMesg("plotter", "add|0;1;2;3;4|0;1;4;9;16|r-+|")
esp.TimMesg("plotter", "add|0;1;2;3;4|4;5;6;7;8|g:x|")
esp.TimMesg("plotter", "show")
esp.TimHold("pyscript", "plotter")

# all line types
esp.TimMesg("plotter", "new|Example of all lines/symbols|x-axis label|y-axis label|")
esp.TimMesg("plotter", "add|1;11;21|2;2;2|c- |")
esp.TimMesg("plotter", "add|1;11;21|4;4;4|m:o|")
esp.TimMesg("plotter", "add|1;11;21|6;6;6|b;x|")
esp.TimMesg("plotter", "add|1;11;21|8;8;8|k_+|")
esp.TimMesg("plotter", "add|1;11;21|10;10;10|r s|")
esp.TimMesg("plotter", "show")
esp.TimHold("pyscript", "plotter")

esp.TimQuit("plotter")
