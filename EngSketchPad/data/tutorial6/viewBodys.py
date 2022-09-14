###################################################################
#                                                                 #
# viewBodys.py - view the various Bodys in the model              #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

# get access to pyCAPS
import pyCAPS
from   pyOCSM import esp

# if running in serveESP, the following load is ignored
myProblem = pyCAPS.Problem(problemName = "RunWithoutESP",
                           capsFile    = "../data/tutorial6/model4.csm",
                           outLevel    = 1)

# load the viewer
esp.TimLoad("viewer", esp.GetEsp("pyscript"), "")

# view all Bodys on the stack
print("==> Viewing all Bodys on the stack")
esp.TimMesg("viewer", "MODL")

# view just the wing Bodys for AVL
print("==> Viewing AVL Bodys with intent=wing")
if ("avlAIM" not in myProblem.analysis):
    myProblem.analysis.create(aim        = "avlAIM",
                              capsIntent = "wing",
                              name       = "avlWing")

#foo = myProblem.analysis["avlWing"].output["Alpha"]
esp.TimMesg("viewer", "AIM|avlWing|")

# view all Bodys going into AVL
print("==> View all AVL Bodys")
if ("avlAIM" not in myProblem.analysis):
    myProblem.analysis.create(aim        = "avlAIM",
                              capsIntent = "",
                              name       = "avlAll")

esp.TimMesg("viewer", "AIM|avlAll|")

esp.TimQuit("viewer")
