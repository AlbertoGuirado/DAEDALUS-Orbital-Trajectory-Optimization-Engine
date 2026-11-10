"""Physical constants of DAEDALUS: one definition each, with unit and source (RD-01).

Every other module imports these names; none redefines a value. Constants
that were defined in dynamics/two_body.py and orbital/mechanics.py before this
module existed are still importable from there (they are re-exports).

All values are SI. They are conventional or model constants, not
measurements with an uncertainty budget; changing one changes the physical
model and must be recorded with every result that depends on it.
"""

# Geocentric gravitational constant GM (TCG-compatible conventional value).
# IERS Conventions (2010), IERS Technical Note 36, Table 1.1; m^3/s^2.
# https://iers-conventions.obspm.fr/content/tn36.pdf
EARTH_MU: float = 3.986004418e14

# Standard acceleration of gravity that defines specific impulse; m/s^2.
# Exact by definition (3rd CGPM, 1901); it is not the local gravity.
STANDARD_GRAVITY: float = 9.80665

# Equatorial radius of the WGS 84 reference ellipsoid; m.
# NIMA TR8350.2, 3rd ed. (2000), Table 3.1. Used as the reference radius of
# altitudes; the body is treated as a sphere of this radius (no flattening).
EARTH_EQUATORIAL_RADIUS: float = 6_378_137.0
