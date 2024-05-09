# mechanism4.py
# written by John Dannenhoffer

# to be use as follows:
#    serveESP ../data/basic/mechanism4a
#       Tools->Pyscript  ../data/basic/mechanism4.py

# import pyCAPS module
import pyCAPS

# load geometry [.csm] file (if run from pthon prompt)
filename = "../data/basic/mechanism4a.csm"
print ('\n==> Loading geometry from file "'+filename+'"...')
myProblem = pyCAPS.Problem(problemName = "workDir_mechanism",
                           capsFile = filename,
                           outLevel = 0)

# show for 20 crank angles (retracted=288 to deployed=174)
print("fully retracted")

for iang in range(21):
    frac  = (iang) / (21-1)
    angle = (1-frac) * 288 + frac * 174
    print("angle =", angle)
    myProblem.geometry.despmtr["angle"].value = angle
    myProblem.geometry.view()

print("fully deployed")

