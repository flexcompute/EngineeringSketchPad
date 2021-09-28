# Make print statement Python 2-3 agnostic
from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
from pyCAPS import capsConvert

# Convert mile to feet 
value = capsConvert(19, "mile", "ft")
print(value, "ft")

# Convert gallon to milliliter 
value = capsConvert(5.0, "gallon", "ml")
print(value, "ml")

# Convert feet to kilograms - expect failure as this doesn't make sense   
try:
    value = capsConvert(1, "ft", "kg")
except: 
    print("Error occurred")

# Convert pounds to kilograms  
value = capsConvert(10.0, "lb", "kg")
print(value, "kg")

# Convert slug to kilogram  
value = capsConvert([1.0, 2.0, 4.0, 10], "slug", "kg")
print(value, "kg")

# Convert BTU (British Thermal Unit) to Joules   
value = capsConvert([[1.0, 2.0], [4.0, 10]], "btu", "J")
print(value, "J")

# Convert foot per second to meter per hour 
value = capsConvert(1, "foot per second", "m/h")
print(value, "m/h")

# Convert foot per second to meter per hour - scale ft/s by 1/6
value = capsConvert(1, "(5 foot)/(30 second)", "m/h")
print(value, "m/h")