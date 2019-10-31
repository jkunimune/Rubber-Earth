import shapefile
import math


SPACING = 30
RADIUS = 4

CENTERS = [(ɸ, λ) for λ in range(-180, 180, SPACING) for ɸ in range(-90+SPACING, 90, SPACING)]
CENTERS.extend([(-90, 0), (90, 0)])


def obliquify(lat1, lon1, lat0, lon0):
	""" go from relative to absolute coordinates """
	lat1, lon1, lat0, lon0 = math.radians(lat1), math.radians(lon1), math.radians(lat0), math.radians(lon0)
	latf = math.asin(math.sin(lat0)*math.sin(lat1) - math.cos(lat0)*math.cos(lon1)*math.cos(lat1))
	innerFunc = math.sin(lat1)/math.cos(lat0)/math.cos(latf) - math.tan(lat0)*math.tan(latf)
	if lat0 == math.pi/2: # accounts for special case when lat0 = pi/2
		lonf = lon1+lon0
	elif lat0 == -math.pi/2: # accounts for special case when lat0 = -pi/2
		lonf = -lon1+lon0 + math.pi
	elif abs(innerFunc) > 1: # accounts for special case when cos(lat1) -> 0
		if (lon1 == 0 and lat1 < -lat0) or (lon1 != 0 and lat1 < lat0):
			lonf = lon0 + math.pi
		else:
			lonf = lon0
	elif math.sin(lon1) > 0:
		lonf = lon0 + math.acos(innerFunc)
	else:
		lonf = lon0 - math.acos(innerFunc)

	while lonf > math.pi:
		lonf -= 2*math.pi
	while lonf < -math.pi:
		lonf += 2*math.pi
		
	latf, lonf = math.degrees(latf), math.degrees(lonf)
	return latf, lonf


if __name__ == '__main__':
	with shapefile.Writer('../../data/tissots_indicatrix_30', shapetype=3) as shp:
		shp.field('Ellipse ID', 'C', size=5)
		for ɸ0, λ0 in CENTERS:
			shp.line([[reversed(obliquify(90-RADIUS, λ, ɸ0, λ0)) for λ in range(0, 360, 5)]])
			shp.record('{:02d}-{:02d}'.format(ɸ0, λ0))
