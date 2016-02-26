#!/home/eamon/anaconda/bin/python
import math

ellipsoids = {
    #Airy-1830: UK, International-1924: Europe, Clarke-1880: Africa
    #GRS-67: South America, GRS-80: America and Australia
    #name        major(m)   minor(m)            flattening factor
    'WGS-84':   (6378137,   6356752.3142451793, 298.25722356300003),
    'Airy-1830': (6377563.396, 6356256.910, 299.3249646),
    'International-1924': (6378388, 6356911.946, 297),
    'Clarke-1880': (6378249.145, 6356514.86955, 293.465),
    'GRS-80':   (6378137,   6356752.3141403561, 298.25722210100002),
    'GRS-67':   (6378160,   6356774.5160907144, 298.24716742700002),
}

def distanceHaversine(lat1, lon1, lat2, lon2):
    '''Author: Chen Meng Date:16/03/2013
    Reference: http://www.movable-type.co.uk/scripts/latlong.html
    Computes the Haversine distance (in kilometers) between two points on the earth.
    '''
    R = 6371 # km
    dLat = math.radians(lat2 - lat1)
    dLon = math.radians(lon2 - lon1)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)

    a = math.sin(dLat/2.0) * math.sin(dLat/2.0) + math.sin(dLon/2.0) * math.sin(dLon/2.0) * math.cos(lat1) * math.cos(lat2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = R * c
    
    y = math.sin(dLon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dLon)
    brng = math.atan2(y,x)
    return d,(brng/(2*math.pi))*360

def distanceVincenty(lat1, long1, lat2, long2, ellipsoid='WGS-84'):
    """From: http://gis.stackexchange.com/questions/32981/calculate-distance-in-km-to-nearest-points-given-in-lat-long-using-arcgis
    slove the bug can not compute same longitudes by Chen Meng
    Add the compution of the forward azimuth by Chen Meng
    Reference: http://www.movable-type.co.uk/scripts/latlong-vincenty.html
    Computes the Vicenty distance (in kilometers) between two points
    on the earth. Coordinates need to be in decimal degrees.
    """
    # Check if we got numbers
    # Removed to save space
    # Check if we know about the ellipsoid
    # Removed to save space
    major, minor, ffactor = ellipsoids[ellipsoid]
    # Convert degrees to radians
    x1 = math.radians(lat1)
    y1 = math.radians(long1)
    x2 = math.radians(lat2)
    y2 = math.radians(long2)
    # Define our flattening f
    f = 1 / ffactor
    # Find delta X
    deltaX = y2 - y1
    # Calculate U1 and U2
    U1 = math.atan((1 - f) * math.tan(x1))
    U2 = math.atan((1 - f) * math.tan(x2))
    # Calculate the sin and cos of U1 and U2
    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)
    # Set initial value of L
    L = deltaX
    # Set Lambda equal to L
    lmbda = L
    lmbdaP = lmbda - 1.0 #make sure lmbda - lmbdaP is not zero, so can start first loop
    #print lmbda
    # Iteration limit - when to stop if no convergence
    iterLimit = 100
    while abs(lmbda-lmbdaP) > 10e-12 and iterLimit >= 0:
        # Calculate sine and cosine of lmbda
        sin_lmbda = math.sin(lmbda)
        cos_lmbda = math.cos(lmbda)
        # Calculate the sine of sigma
        sin_sigma = math.sqrt(
                (cosU2 * sin_lmbda) ** 2 + 
                (cosU1 * sinU2 - 
                 sinU1 * cosU2 * cos_lmbda) ** 2
        )
	#print sin_sigma
        if sin_sigma == 0.0:
            # Concident points - distance is 0
            return 0.0
        # Calculate the cosine of sigma
        cos_sigma = (
                    sinU1 * sinU2 + 
                    cosU1 * cosU2 * cos_lmbda
        )
        # Calculate sigma
        sigma = math.atan2(sin_sigma, cos_sigma)
        # Calculate the sine of alpha
        sin_alpha = (cosU1 * cosU2 * math.sin(lmbda)) / (sin_sigma)
	#print sin_alpha
        # Calculate the square cosine of alpha
        cos_alpha_sq = 1 - sin_alpha ** 2
	#print cos_alpha_sq
        # Calculate the cosine of 2 sigma
        cos_2sigma = cos_sigma - ((2 * sinU1 * sinU2) / cos_alpha_sq)
        # Identify C
        C = (f / 16.0) * cos_alpha_sq * (4.0 + f * (4.0 - 3 * cos_alpha_sq))
	lmbdaP = lmbda
        # Recalculate lmbda now
        lmbda = L + ((1.0 - C) * f * sin_alpha * (sigma + C * sin_sigma * (cos_2sigma + C * cos_sigma * (-1.0 + 2 * cos_2sigma ** 2)))) 
        # If lambda is greater than pi, there is no solution
    	if (abs(lmbda) > math.pi):
            raise ValueError("No solution can be found.")
        iterLimit -= 1

    if iterLimit == 0 and lmbda > 10e-12:
        raise ValueError("Solution could not converge.")
    # Since we converged, now we can calculate distance
    # Calculate u squared
    u_sq = cos_alpha_sq * ((major ** 2 - minor ** 2) / (minor ** 2))
    # Calculate A
    A = 1 + (u_sq / 16384.0) * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)))
    # Calculate B
    B = (u_sq / 1024.0) * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)))
    # Calculate delta sigma
    deltaSigma = B * sin_sigma * (cos_2sigma + 0.25 * B * (cos_sigma * (-1.0 + 2.0 * cos_2sigma ** 2) - 1.0/6.0 * B * cos_2sigma * (-3.0 + 4.0 * sin_sigma ** 2) * (-3.0 + 4.0 * cos_2sigma ** 2)))
    # Calculate s, the distance
    s = minor * A * (sigma - deltaSigma)
    alpha1 = math.atan2(cosU2 * sin_lmbda, cosU1 * sinU2 - sinU1 * cosU2 * cos_lmbda)
    #alpha2 = math.atan2(cosU1 * sin_lmbda, cosU1 * sinU2 * cos_lmbda - sinU1 * cosU2)
    # Return the distance
    return s/1000.0,(alpha1/(2*math.pi))*360
