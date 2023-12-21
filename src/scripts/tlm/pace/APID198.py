import numpy as np
from .. PacketUtils import *
from .. timestamp import *

'''
Tilt:  APID 198
'''

Fields = np.dtype([
  ('timestamp', CCSDS_timestamp),
  ('Padding',   '>u4'),   # Padding to maintain packet alignment
  ('words',     '>u2', 32),    # Many words are packed bits; unpack later
])      # total length = 74 bytes

def APID198(data):
    apid = 198
    check_packet_size(apid, Fields.itemsize, len(data))

    tmp = np.frombuffer(data, dtype=Fields, count=1)
    myDict = {}
    myDict['timestamp'] = parse_CCSDS_timestamp(tmp['timestamp'][0])

    # unpack each 16-bit word
    words = tmp['words'][0]

    # WORD_00
    bits = getbits(words, 0)
    #myDict['RSVD_00']    = as_int(bits[0:0+2])   # Reserved bits 0-1
    myDict['MCE_REV']     = as_int(bits[2:2+2])   # MCE Revision
    myDict['FPGA_VER']    = as_int(bits[5:5+12])  # FPGA Version

    # WORD_01: SADA state
    bits = getbits(words, 1)
    myDict['SADA_SPCL_LOC']       = as_int(bits[0:0+2])   # 00=Not Special Loc 01=Home 10=NEG EOT 11=POS EOT
    #myDict['SADA_STP_ERR_FLG']   = bits[2]               # NOT USED FOR PACE Step Error pace_tsade_sht_dscFlag, 00=No Error, 01=Error
    myDict['SADA_RLY_STATUS']     = bits[3]               # 0=Open Relay, 1=Closed Relay
    myDict['SADA_STEP_PROFILE']   = as_int(bits[4:4+2])   # 00=Cardinal Stepping, 01=Smoothed Cardinal Stepping, 10=Sinusoidal
    myDict['SADA_CUR_DIR']        = bits[6]               # 0=Negative, 1=Positive
    myDict['SADA_LAST_COMM_STATE']= as_int(bits[7:7+3])   # 001=State 1, 010=State 2, 011=State 3, 100=Sate 4, 101=State 5, 110 State 6, 000 & 111=Not Valid
    myDict['SADA_LAST_STP_HOLD_TIME_STAT'] = bits[10]     # 0=Disabled (mtr pwr on full time), 1=Enabled (mtr pwr off @ end of step)
    myDict['SADA_H_BRIDGE_BOOST'] = bits[11]              # 0=Disabled, 1=Enabled
    myDict['SADA_PWR_LVL']        = as_int(bits[12:12+2]) # 00=50% maxV, 01=60% maxV, 10=80% maxV, 11=100% maxV
    myDict['SADA_HOME_INIT']      = bits[14]              # 0=Encoder not initialized, 1=initialized
    myDict['SADA_CNT_ENA']        = bits[15]              # Count enable

    # WORD_02 - WORD_03:  Motor 0 (SADA) position
    myDict['SADA_ENCPOS'] = words[2]                 # Motor 0 (SADA) reported encoder pos
    myDict['SADA_CMDPOS'] = words[3]                 # Motor 0 (SADA) Commanded pos,
    # x'3C4A'--x'7FFF' (-130 to -0.0075 degrees from home), x'8000' = home, x'8001--x'C3B6 = 0.0075 to 130

    # WORD_04 - WORD_06:  Motor 0 (SADA) Pulse
    myDict['SADA_PULSE_RATE']            = words[4]  # Motor 0 (SADA) Pulse rate
    myDict['SADA_PULSE_WIDTH']           = words[5]  # Motor 0 (SADA) Pulse Width
    myDict['SADA_LAST_STP_PULSE_LENGTH'] = words[6]  # Motor 0 (SADA) Last Step Pulse Length
    # x'0000'--x'007C' not valid, x'007D'--x'FFFF' 5msec to 2.62114 sec, # of 40usec clock cycles

    # WORD_07:  Motor 0 (SADA) flags
    bits = getbits(words, 7)
    myDict['SADA_INDEX']  = bits[0] # Motor 0 (SADA) Index Track,
    # D'15=Index track, D'14=Sine track, D'13=Cosine track, D'12=Hemisphere track, D'11=Almost Home track, D'10=Tavle limits track
    myDict['SADA_SIN']            = bits[1]       # Motor 0 (SADA) Sine Track
    myDict['SADA_COS']            = bits[2]       # Motor 0 (SADA) Cosine Track
    myDict['SADA_HEMI']           = bits[3]       # Motor 0 (SADA) Hemisphere Track
    myDict['SADA_ALMOST_HOME']    = bits[4]       # Motor 0 (SADA) Almost Home
    myDict['SADA_TRVL_LIM']       = bits[5]       # Motor 0 (SADA) Travel Limits
    myDict['SADA_EOT_OVRD_STAT']  = as_int(bits[6:6+2])   # Motor 0 (SADA) EOT Override Status
    #myDict['SADA_W7_RSVD']       = bits[8]       # Motor 0 (SADA) Reserved
    myDict['SADA_RAW_HEMI']       = bits[9]       # Motor 0 (SADA) Raw Hemi
    myDict['SADA_PROCESSED_HEMI'] = bits[10]      # Motor 0 (SADA) Processed Hemi
    myDict['SADA_INVERTED_HEMI']  = bits[11]      # Motor 0 (SADA) Inverted Hemi
    #myDict['SADA_W3_RSVD']       = bits[12]      # Motor 0 (SADA) Reserved
    myDict['SADA_LED_CMD_PWR_LVL']= as_int(bits[13:13+3]) # Motor 0 (SADA) LED Commanded Power Level, 0-011=commanded current lvl, 100-111=encoder off

    # WORD_08: SADA current
    myDict['SADA_DRVR_CUR'] = words[8].item()     # Motor 0 LED Driver Current, 1cnt=0.015259 mA

    # WORD_09: TILT state
    bits = getbits(words, 9)
    myDict['TILT_SPCL_LOC']       = as_int(bits[0:0+2])   # 00=Not a special location, 01=at home, 10=CW end of travel, 11, CCW end of travel
    #myDict['TILT_STP_ERR_FLG']   = bits[2]               # NOT USED FOR PACE Step Error Flag, 00=no error, 01=error
    myDict['TILT_RLY_STATUS']     = bits[3]               # 0=open (disabled), 1=closed (enabled)
    myDict['TILT_STEP_PROFILE']   = as_int(bits[4:4+2])   # 00=cardinal step, 01=smoothed cardinal stepping, 10=Sinusoidal
    myDict['TILT_CUR_DIR']        = bits[6]               # 0=Negative, 1=Positive
    myDict['TILT_LAST_COMM_STATE']= as_int(bits[7:7+3])   # 001=State 1, 010=State 2, 011=State 3, 100=Sate 4, 101=State 5, 110 State 6, 000 & 111=Not Valid
    myDict['TILT_LAST_STP_HOLD_TIME_STAT'] = bits[10]     # 0=Disabled (mtr pwr on full time), 1=Enabled (mtr pwr off @ end of step)
    myDict['TILT_H_BRIDGE_BOOST'] = bits[11]              # 0=Disabled, 1=Enabled
    myDict['TILT_PWR_LVL']        = as_int(bits[12:12+2]) # 00=50% maxV, 01=60% maxV, 10=80% maxV, 11=100% maxV
    myDict['TILT_HOME_INIT']      = bits[14]              # 0=Encoder not initialized, 1=initialized
    myDict['TILT_CNT_ENA']        = bits[15]              # Count enable

    # WORD_10 - WORD_11:  Motor 1 (TILT) position
    (myDict['TILT_ENCPOS'],                           # Motor 1 (TILT) reported encoder pos
     myDict['TILT_CMDPOS'],                           # Motor 1 (TILT) Commanded pos,
    ) = 0.0075 * (words[10:12].astype(int) - 32768)
    #) = (0.0075 * words[10:12]) - 245.76
    # x'3C4A'--x'7FFF' (-130 to -0.0075 degrees from home), x'8000' = home, x'8001--x'C3B6 = 0.0075 to 130
    # PolynomialConversion pace_tsade_enc_tilt           {coefficients={-245.76          , 0.0075}}

    # WORD_12 - WORD_14:  Motor 1 (TILT) Pulse
    myDict['TILT_PULSE_RATE']            = words[12]  # Motor 1 (TILT) Pulse rate
    myDict['TILT_PULSE_WIDTH']           = words[13]  # Motor 1 (TILT) Pulse Width
    myDict['TILT_LAST_STP_PULSE_LENGTH'] = words[14]  # Motor 1 (TILT) Last Step Pulse Length
    # x'0000'--x'007C' not valid, x'007D'--x'FFFF' 5msec to 2.62114 sec, # of 40usec clock cycles

    # WORD_15:  Motor 1 (TILT) flags
    bits = getbits(words, 15)
    myDict['TILT_INDEX']  = bits[0]       # Motor 1 (TILT) Index Track,
    # D'15=Index track, D'14=Sine track, D'13=Cosine track, D'12=Hemisphere track, D'11=Almost Home track, D'10=Tavle limits track
    myDict['TILT_SIN']            = bits[1]       # Motor 1 (TILT) Sine Track
    myDict['TILT_COS']            = bits[2]       # Motor 1 (TILT) Cosine Track
    myDict['TILT_HEMI']           = bits[3]       # Motor 1 (TILT) Hemisphere Track
    myDict['TILT_ALMOST_HOME']    = bits[4]       # Motor 1 (TILT) Almost Home
    myDict['TILT_TRVL_LIM']       = bits[5]       # Motor 1 (TILT) Travel Limits
    myDict['TILT_EOT_OVRD_STAT']  = as_int(bits[6:6+2])   # Motor 1 (TILT) EOT Override Status
    #myDict['TILT_W15_RSVD']      = bits[8]       # Motor 1 (TILT) Reserved
    myDict['TILT_RAW_HEMI']       = bits[9]       # Motor 1 (TILT) Raw Hemi
    myDict['TILT_PROCESSED_HEMI'] = bits[10]      # Motor 1 (TILT) Processed Hemi
    myDict['TILT_INVERTED_HEMI']  = bits[11]      # Motor 1 (TILT) Inverted Hemi
    #myDict['TILT_W3_RSVD']       = bits[12]      # Motor 1 (TILT) Reserved
    myDict['TILT_LED_CMD_PWR_LVL']= as_int(bits[13:13+3]) # Motor 1 (TILT) LED Commanded Power Level

    # WORD_16: TILT current
    myDict['TILT_DRVR_CUR'] = words[16].item()    # Motor 1 LED Driver Current, 1cnt=0.015259 mA

    # WORD_17  (First 5 bits indicate type in each CMD_ECHO.  Not a mistake.)
    bits = getbits(words, 17)
    myDict['CMD_ECHO_TYPE']       = as_int(bits[0:0+5])

    # WORDs 17 - 25
    (myDict['CMD_ECHO_0'],
     # Command Echo 0, 9 CMD words from last 1553 received cmd to the MCE RT, valid and invalid cmds are echoed
     myDict['CMD_ECHO_1'],        # Command Echo 1, ...
     myDict['CMD_ECHO_2'],        # Command Echo 2, ...
     myDict['CMD_ECHO_3'],        # Command Echo 3, ...
     myDict['CMD_ECHO_4'],        # Command Echo 4, ...
     myDict['CMD_ECHO_5'],        # Command Echo 5, ...
     myDict['CMD_ECHO_6'],        # Command Echo 6, ...
     myDict['CMD_ECHO_7'],        # Command Echo 7, ...
     myDict['CMD_ECHO_8'],        # Command Echo 8, ...
    ) = words[17:26]

    # WORD_26
    bits = getbits(words, 26)
    myDict['CMDS_NUMBER_SINCE_START'] = as_int(bits[0:8])
    # Number of i cmds received by MCE, x'0000' - x'FFFF'], can be cmd type or invalid parameter, or cmd when mtr not initialized
    myDict['CMDS_VLD_SINCE_START']    = as_int(bits[8:16])  # Number of valid cmds received by MCE, x'0000' - x'FFFF'

    # WORD_27
    bits = getbits(words, 27)
    #myDict['W27_RSVD']           = as_int(bits[0:0+2])   # Reserved bits 15-14
    myDict['TILT_SADA_CURRENT_AF']= bits[2]       # TILT & SADA Current action flag, if true, disables H-Bridge
    myDict['TILT_SADA_CURRENT_AE']= bits[3]       # TILT & SADA Current action enable, responds to error
    myDict['WD_WRITE_AF']         = bits[4]       # Watch Dog Write Action Flag
    myDict['WD_READ_AF']          = bits[5]       # Watch Dog Read Action Flag
    myDict['WD_WRITE_EF']         = bits[6]       # Watch Dog Write Error Flag
    myDict['WD_READ_EF']          = bits[7]       # Watch Dog Read Error Flag
    myDict['SADA_EOT_PGM_POS_EF'] = bits[8]       # SADA EOT Programable Postive Error Flag
    myDict['SADA_EOT_PGM_NEG_EF'] = bits[9]       # SADA EOT Programable Negative Error Flag
    myDict['TILT_EOT_PGM_POS_EF'] = bits[10]      # TILT EOT Programable Postive Error Flag
    myDict['TILT_EOT_PGM_NEG_EF'] = bits[11]      # TILT EOT Programable Negative Error Flag
    myDict['TILT_ERR_FLG']        = as_int(bits[12:12+2]) # Motor 1 (TILT) error flag, 00=no error, 01=No home detected, 10=No right position
    #myDict['W27_RSVD2']          = bits[14]      # Motor 1 (TILT) time contorl, 1=time cntrl enabled, 0=time cntrl disabled
    myDict['TILT_DIR']            = bits[15]      # Motor 1 (TILT) direction, 1=CW, 0=CCW

    # WORD_28
    bits = getbits(words, 28)
    myDict['ERR_CTR_NO_HOME_WHEN_TILT']           = as_int(bits[0:0+4])  # Error Counter no home detected when tilt
    myDict['ERR_CTR_NIT_NOT_INIT_POS_WHTILT']     = as_int(bits[4:4+4])  # Error Counter not init position when tilt
    myDict['NUM_OF_TILT_FSM_CMDS_DECODED']        = as_int(bits[8:8+8])  # Number of TILT FSM Commands Decoded

    # WORD_29
    (myDict['ABS_TILT_DELTA_ERR_POS'],            # Motor 1 (TILT) absolute delta error position
     myDict['ABS_SADA_DELTA_ERR_POS'],            # Motor 0 (SADA) absolute delta error position
    ) = words[29,np.newaxis].view(dtype='|u1')

    # WORD_30
    (myDict['SADA_POS_ED_EF'],            # SADA ED Error
     myDict['SADA_POS_ED_ENA'],
     myDict['TILT_SADA_CUR_EF'],
     myDict['TILT_EOT_LED_NEG_EF'],
     myDict['TILT_EOT_LED_POS_EF'],
     myDict['SADA_EOT_LED_NEG_EF'],
     myDict['SADA_EOT_LED_POS_EF'],
     myDict['TILT_EOT_ENC_NEG_ENA'],
     myDict['TILT_EOT_ENC_POS_ENA'],
     myDict['SADA_EOT_ENC_NEG_ENA'],
     myDict['SADA_EOT_ENC_POS_ENA'],
     myDict['TILT_DELTA_ERF'],
     myDict['TILT_TOTAL_TRAVEL_EF'],
     myDict['TILT_CRUISE_EF'],
     myDict['TILT_RAMP_DWN_EF'],
     myDict['TILT_RAMP_UP_EF'],
    ) = getbits(words, 30)

    # WORD_31
    (myDict['SADA_POS_ED_AE'],
     myDict['TILT_SADA_CUR_ENA'],
     myDict['TILT_DELTA_AE'],             # Motor 1 (TILT) Delta Action Enable
     myDict['TILT_TOTAL_TIME_AE'],        # Motor 1 (TILT) Total Time Action Enable
     myDict['TILT_CRUISE_AE'],            # Motor 1 (TILT) Crusie Action Enable
     myDict['TILT_RAMPDOWN_AE'],          # Motor 1 (TILT) Ramp Down Action Enable
     myDict['TILT_RAMPUP_AE'],            # Motor 1 (TILT) Ramp Up Action Enable
     myDict['TILT_DELTA_ENA'],
     myDict['TILT_TOT_TIME_ENA'],
     myDict['TILT_CRUISE_ENA'],
     myDict['TILT_RAMPDOWN_ENA'],
     myDict['TILT_RAMPUP_ENA'],
     myDict['WD_W_ERR_ENA'],      # Enables TSADE Error Detection for write communication to 1553 Bus
     myDict['WD_R_ERR_ENA'],      # Enables TSADE Error Detection for read communication to 1553 Bus
     myDict['WD_W_ACTION_ENA'],   # Enables TSADE Action if the WD Write Error Detection is triggered
     myDict['WD_R_ACTION_ENA'],   # Enables TSADE Action if the WD Read Error Detection is triggered
    ) = getbits(words, 31)

    return myDict
