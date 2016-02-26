# geodistance
A python library for comupting distance and azimuth

Useage:

Haversine distance:
distance, azimuth = geodistance.distanceHaversine(lat1, lon1, lat2, lon2)

Vicenty distance:
distance, azimuth = geodistance.distanceVincenty(lat1, long1, lat2, long2, ellipsoid='WGS-84')
