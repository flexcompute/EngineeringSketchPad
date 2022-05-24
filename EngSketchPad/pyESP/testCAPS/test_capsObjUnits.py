import unittest

import os
import math

from pyCAPS import caps

from sys import version_info as pyVersion

class TestUnit(unittest.TestCase):

#=============================================================================-
    def test_convert(self):
        
        m1k = 1000*caps.Unit("m")
        km  = m1k.convert("km")

        self.assertEqual(km._value, 1)
        self.assertEqual(km._units, caps.Unit("km"))

#=============================================================================-
    def test_compare(self):

        m1k = caps.Unit("1000m")
        km  = caps.Unit("km")
        kg  = caps.Unit("kg")

        self.assertEqual(m1k, km)
        self.assertNotEqual(kg, km)

#=============================================================================-
    def test_binaryOp(self):
        
        s   = caps.Unit("s")
        s2  = caps.Unit("s^2")
        m   = caps.Unit("m")
        km  = caps.Unit("km")
        kg  = caps.Unit("kg")
        kgm = caps.Unit("kg.m")
        N   = caps.Unit("N")

        self.assertEqual(kgm, kg*m)
        self.assertEqual(s2, s**2)

        self.assertEqual(N, kg*m/s**2)

        m2  = 2*caps.Unit("m")
        m4  = 4*caps.Unit("m")
        self.assertEqual(m2, 1*m+1*m)
        self.assertEqual(1*m, m2-1*m)
        self.assertEqual(m2, 2*m)
        self.assertEqual(1*m, m2/2)

        with self.assertRaises(caps.CAPSError) as e:
            bad = m+s
        self.assertEqual(e.exception.errorName, "CAPS_UNITERR")

        with self.assertRaises(caps.CAPSError) as e:
            bad = m-s
        self.assertEqual(e.exception.errorName, "CAPS_UNITERR")

        self.assertEqual(caps.Unit(1)   , m/m)
        self.assertEqual(caps.Unit(1000), km/m)
        self.assertEqual(caps.Unit(1)   , m*m**-1)
        self.assertEqual(caps.Unit(1000), km*m**-1)


#=============================================================================-
    def test_offset(self):

        C = caps.Unit("celsius")
        K = caps.Unit("kelvin")
        
        self.assertEqual((1*C).convert("kelvin"), 274.15*K)
        self.assertEqual(1*C, 1*K + 273.15*K)
        self.assertEqual(C, caps.Unit("K @ 273.15"))

        kPa = caps.Unit("kPa")
        kg = caps.Unit("kg")
        m = caps.Unit("m")
        s = caps.Unit("s")
        
        rhoUnit = caps.Unit(1.225)*kg/m**3
        V_Unit = caps.Unit(10)*m/s
        dynpUnit = rhoUnit*V_Unit**2

        pinf = 101.3*kPa
        p = 0.99*pinf
        
        # Add offset to the unit
        CpUnit = dynpUnit + pinf/dynpUnit
        
        # Convert the pressure to Cp, and back to pressure
        Cp = p.convert(CpUnit)
        self.assertAlmostEqual(p, Cp.convert("Pa"), 5)

#         print()
#         print("Cp ", Cp)
#         print("1/dynp ", 1/dynpUnit)
#         print("p ", p)
#         print("(p - pinf)/dynp ", (p - pinf)/dynpUnit)
#         print("p.convert(Cp) ", p.convert(CpUnit))
#         print("Cp.convert(kPa) ", (-8.269387755102116*CpUnit).convert("kPa"))

#=============================================================================-
    def test_unaryOp(self):
        m   = caps.Unit("m")
        nm  = caps.Unit("-1 m")

        self.assertEqual(nm, -m)
        self.assertEqual( m, +m)
        
        self.assertEqual(abs(-1*m), +1*m)
        
        if pyVersion.major > 2:
            length = 9.123456*m
            self.assertEqual(round(length, 6), 9.123456*m)
            self.assertEqual(round(length, 5), 9.12346 *m)
            self.assertEqual(round(length, 4), 9.1235  *m)
            self.assertEqual(round(length, 3), 9.123   *m)
            self.assertEqual(round(length, 2), 9.12    *m)
            self.assertEqual(round(length, 1), 9.1     *m)
            self.assertEqual(round(length, 0), 9.      *m)

#=============================================================================-
    def test_inplaceOp(self):
        m  = caps.Unit("m")
        s  = caps.Unit("s")

        unit = 1*m

        unit += 1*m
        self.assertEqual(2*m, unit)
        unit -= 1*m
        self.assertEqual(1*m, unit)
        unit /= 2*s
        self.assertEqual(0.5*m/s, unit)
        unit **= 2
        self.assertEqual((0.5*m/s)**2, unit)

#=============================================================================-
    def test_angle(self):
        rad = caps.Unit("rad")
        deg = caps.Unit("degree")
        
        self.assertAlmostEqual(1*deg/rad, math.pi/180., 5)
        self.assertAlmostEqual(math.sin(5*deg), math.sin(5*deg/rad), 5)

#=============================================================================-
    def test_value(self):
        m  = caps.Unit("m")
        ft  = caps.Unit("ft")
        
        q = 10*m # Make a Quantity
        
        self.assertEqual(10, q/m)
        self.assertEqual(10, q.value())
        self.assertAlmostEqual(q.convert(ft).value(), q/ft, 5)
 
#=============================================================================-
    def test_list(self):
        m  = caps.Unit("m")
        s  = caps.Unit("s")

        lengths = [3, 4]*m
        
        self.assertEqual(lengths[0], 3*m)
        self.assertEqual(lengths[1], 4*m)

        speeds = lengths/(2.*s)

        self.assertEqual(speeds[0], 1.5*m/s)
        self.assertEqual(speeds[1], 2.0*m/s)

        accels = lengths/(2.*s**2)

        self.assertEqual(accels[0], 1.5*m/s**2)
        self.assertEqual(accels[1], 2.0*m/s**2)

if __name__ == '__main__':
    unittest.main() 
