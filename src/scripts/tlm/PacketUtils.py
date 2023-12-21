import sys
import numpy as np
from tlm.timestamp import *

def getDict(structured_array):
    '''
    convert numpy structured array to dict of scalars and numpy arrays
    '''
    myDict = {}

    for key in structured_array.dtype.names:
        if ( key.upper().startswith('PADDING') or
             key.upper().startswith('SPARE') ):
            continue  # skip meaningless fields
        myDict[key] = structured_array[key][0]

    return myDict


def getbits(var,i=0):
    return np.unpackbits(var[i,np.newaxis].view(dtype='|u1'))

def as_int(bitlist):
    # convert an arbitrary-length list of bits to integer
    result = 0
    for bit in bitlist:
        result = (result << 1) | bit
    return result

def pad_packet(data, length):
    # pad or truncate bytearray to expected length
    packet = bytearray(length)
    if len(data) > length:
        packet[:] = data[0:length]
    else:
        packet[0:len(data)] = data
    return packet

def check_packet_size(apid,expected,actual):
    if actual < expected:
        message = 'Error! packet too small; exiting.'
        print('APID {}: expected={}, actual={}\t{}'.format(apid,expected,actual,message))
        sys.exit(1)
    elif actual > expected:
        message = 'Warning! not all data will be parsed.'
        print('APID {}: expected={}, actual={}\t{}'.format(apid,expected,actual,message))

#-------------------------------------------------------------------------------

cfeHeaderFields = np.dtype([
    ('filetype', 'S4'),              # "type of file = "cFE1"
    ('subtype', '>u4'),              # 101d = downloaded PACE DS files
    ('length', '>u4'),               # length of CFE file header = 64 bytes
    ('SCID', 'S4'),                  # Spacecraft ID
    ('APID', '>u4'),                 # Application ID
    ('processorid', '>u4'),          # 1=Spacecraft, 2=OCI  (30 in ETE files)
    ('fileopen_seconds', '>u4'),     # Time when file created
    ('fileopen_subseconds', '>u4'),  # Time when file created (fractional seconds)
    ('description', 'S32'),          # = "DS data storage file"
])                                   # total length = 64 bytes

dsfHeaderFields = np.dtype([         # CFE/S Data Storage header
    ('fileclose_seconds', '>u4'),    # Time when file closed
    ('fileclose_subseconds', '>u4'), # Time when file closed (fractional seconds)
    ('filetable_index', '>u2'),      # 0=events, 1=housekeeping
    ('filename_type', '>u2'),        # 1=count, 2=time (PACE uses count)
    ('filename', 'S64'),             # fully qualified file path
])                                   # total length = 76 bytes

def readFileHeader(filehandle):

    # is a CFE header present?
    pos = filehandle.tell()
    fourchars = filehandle.read(4)
    filehandle.seek(pos) # rewind to original position
    if fourchars != b'cFE1':
        return None

    # read CFE header
    data = filehandle.read(cfeHeaderFields.itemsize)
    tmp = np.frombuffer(data, dtype=cfeHeaderFields, count=1)
    myDict = getDict(tmp)
    myDict['fileopen_ts'] = decode_timestamp(myDict['fileopen_seconds'],myDict['fileopen_subseconds'])
    myDict['fileopen'] = datetime_repr(tai58_as_datetime(myDict['fileopen_ts']))

    # append DS header, if it's there.
    if myDict['description'].startswith(b'DS'):
        data = filehandle.read(dsfHeaderFields.itemsize)
        tmp = np.frombuffer(data, dtype=dsfHeaderFields, count=1)
        myDict.update(getDict(tmp))
        myDict['fileclose_ts'] = decode_timestamp(myDict['fileclose_seconds'],myDict['fileclose_subseconds'])
        myDict['fileclose'] = datetime_repr(tai58_as_datetime(myDict['fileclose_ts']))

    return myDict

#-------------------------------------------------------------------------------

primaryHeaderFields = np.dtype([
    ('words', '>u2', 3),  # Unpack bits later
])                        # total length = 6 bytes

def primaryHeader(data):

    # unpack each 16-bit word
    tmp = np.frombuffer(data, dtype=primaryHeaderFields, count=1)
    words = tmp['words'][0]
    myDict = {}

    # WORD_00
    bits = getbits(words, 0)
    myDict['version']   = as_int(bits[0:0+3])  # Packet Version Number
    myDict['ptype']     = bits[3]              # Packet Type (0=tlm;1=cmd)
    myDict['secondary'] = bits[4]              # Secondary Header Flag
    myDict['APID']      = as_int(bits[5:5+11]) # Application Process Identifier

    # WORD_01
    bits = getbits(words, 1)
    myDict['grouping']  = as_int(bits[0:0+2])  # Sequence (00=cont., 01=first, 10=last, 11=only)
    myDict['sequence']  = as_int(bits[2:2+11]) # Packet Sequence Count

    # WORD_02
    myDict['length']    = words[2]             # Packet Data Length (total bytes - 1)

    # no conversions needed
    return myDict

#-------------------------------------------------------------------------------

def readPacket(filehandle):

    # Read packetField 1: packet primary header
    try:
        rawhdr = filehandle.read(primaryHeaderFields.itemsize)
        header = primaryHeader(rawhdr)
        datalen = header['length']
        assert datalen > 0
    except:
        return None, None, None

    # Read packetField 2: data
    data = filehandle.read(datalen+1)

    return header, data, rawhdr
