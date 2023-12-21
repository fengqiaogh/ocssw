import numpy as np
from .. PacketUtils import *
from .. timestamp import *

'''
Ephemeris:  APID 128
OM_Ephem_Tlm packet (section 7.1)
'''

Fields = np.dtype([
  ('timestamp', CCSDS_timestamp),
  ('Padding',   '>u4'),   # Padding to maintain packet alignment

  ('FSWTime',      '>f8'),        # seconds since TAI epoch
  ('sunPosJ2000',  '>f8', (3,)),  # sun  position wrt ECI (m)
  ('moonPosJ2000', '>f8', (3,)),  # moon position wrt ECI (m)
  ('sunVelJ2000',  '>f8', (3,)),  # sun  velocity wrt ECI (m/sec)
  ('moonVelJ2000', '>f8', (3,)),  # moon velocity wrt ECI (m/sec)
  ('scPosJ2000',   '>f8', (3,)),  # PACE position wrt ECI (m)
  ('scVelJ2000',   '>f8', (3,)),  # PACE velocity wrt ECI (m/sec)
  ('PV_Time_GPS',  '>f8'),        # timestamp relative to GPS epoch (sec)
  ('DCM_ecef2eci', '>f8', (3,3)), # Direction cosine matrix converting from ECEF to ECI J2000 frame
  ('magfieldEci',  '>f8', (3,)),  # Earth mag field ECI (T)
  ('magfieldEcef', '>f8', (3,)),  # Earth mag field ECEF (T)
  ('isEclipse',    '|u1'),        # boolean in eclipse flag
  ('ephemValid',   '|u1'),        # spacecraft ephemeris validity flag
  ('Padding2',     '|u1', (6,)),

]) # total length = 298 bytes

def APID128(data):
  apid = 128
  check_packet_size(apid, Fields.itemsize, len(data))

  tmp = np.frombuffer(data, dtype=Fields, count=1)
  myDict = getDict(tmp)
  myDict['timestamp'] = parse_CCSDS_timestamp(myDict['timestamp'])

  return myDict

def derive_orbitparams(orb_recordlist):

    # adapted from Fred Patt's IDL routine get_ephem_from_hkt.pro
    posr = [np.matmul(rec['DCM_ecef2eci'], rec['scPosJ2000']) for rec in orb_recordlist]
    velr = [np.matmul(rec['DCM_ecef2eci'], rec['scVelJ2000']) for rec in orb_recordlist]

    omegae = 7.29211585494e-5
    for i in np.arange(len(orb_recordlist)):
        velr[i][0] += posr[i][1]*omegae
        velr[i][1] -= posr[i][0]*omegae

    # adapted from Fred Patt's IDL routine orb2lla.pro
    re = 6378.137      # Earth equatorial radius (km)
    rem = 6371.0       # Earth mean radius (km)
    f = 1./298.257     # Earth flattening factor
    omf2 = (1.0-f) * (1.0-f)

    xyz = np.array([p/1000.0 for p in posr], dtype=np.float64)  # convert meters to km
    x, y, z = [xyz[:,i] for i in np.arange(3)]  # separate coords

    # Compute longitude
    lon = np.arctan2(y, x)

    # Compute geodetic latitude
    rad = np.linalg.norm(xyz,axis=1)  # Euclidean distance
    omf2p = (omf2*rem + rad - rem)/rad
    pxy = x*x + y*y
    temp = np.sqrt(z*z + omf2p*omf2p*pxy)
    lat = np.arcsin(z/temp)

    # Compute altitude
    clatg = np.cos(np.arctan(omf2*np.tan(lat)))
    rl = re*(1.0-f)/np.sqrt(1.0-(2.0-f)*f*clatg*clatg)
    alt = rad - rl

    # return as dictionary
    orbitparams = {}
    orbitparams['posr'] = posr
    orbitparams['velr'] = velr
    orbitparams['lon']  = np.rad2deg(lon)
    orbitparams['lat']  = np.rad2deg(lat)
    orbitparams['alt']  = alt * 1000.0  # convert from km to meters

    return orbitparams
