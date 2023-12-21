import numpy as np
from .. PacketUtils import *
from .. timestamp import *

'''
Attitude:  APID 108
AC_AttFilterOut_tlm_msg packet (section 7.2)

Revised 2022-03-03 for updated packet structure
'''

Fields = np.dtype([
  ('timestamp', CCSDS_timestamp),
  ('Padding',   '>u4'),   # Padding to maintain packet alignment

  ('Est_Time',          '>f8'),         # estimate TAI timestamp
  ('q_EciToBrf_Est',    '>f8', (4,)),   # attitude estimate
  ('b_gyro_Brf_Est',    '>f8', (3,)),   # gyro bias estimate
  ('w_EciToBrf_Brf_Est','>f8', (3,)),   # angular rate estimate
  ('Covar',             '>f8', (6,6)),
  ('Covar_eigs',        '>f8', (6,)),   # covariance matrix eigenvalues
  ('Covar_cond',        '>f8'),         # covariance matrix condition number

  ('st1_latest_accept_time',    '>f8'), # timestamp of most recent ST1 accepted
  ('st1_latest_accept_q_EciToBrf',    '>f8', (4,)), # Most recent ST1 measurement accepted by MEKF
  ('st1_latest_mahalanobis',    '>f8'), # (Mahalanobis)^2 of most recent ST1 accepted
  ('st1_latest_reject_time',    '>f8'), # timestamp of most recent ST1 rejected
  ('st1_latest_reject_q_EciToBrf',    '>f8', (4,)), # Most recent ST1 measurement rejected by MEKF

  ('st2_latest_accept_time',    '>f8'), # timestamp of most recent ST2 accepted
  ('st2_latest_accept_q_EciToBrf',    '>f8', (4,)), # Most recent ST2 measurement accepted by MEKF
  ('st2_latest_mahalanobis',    '>f8'), # (Mahalanobis)^2 of most recent ST2 accepted
  ('st2_latest_reject_time',    '>f8'), # timestamp of most recent ST2 rejected
  ('st2_latest_reject_q_EciToBrf',    '>f8', (4,)), # Most recent ST2 measurement rejected by MEKF

  ('st3_latest_accept_time',    '>f8'), # timestamp of most recent ST3 accepted
  ('st3_latest_accept_q_EciToBrf',    '>f8', (4,)), # Most recent ST3 measurement accepted by MEKF
  ('st3_latest_mahalanobis',    '>f8'), # (Mahalanobis)^2 of most recent ST3 accepted
  ('st3_latest_reject_time',    '>f8'), # timestamp of most recent ST3 rejected
  ('st3_latest_reject_q_EciToBrf',    '>f8', (4,)), # Most recent ST3 measurement rejected by MEKF

  ('st1_latest_innov',          '>f4', (6,)), # Innovation (delta Gibbs and delta bias) computed from most recent ST1 measurement
  ('st1_latest_prefit_res',     '>f4', (3,)), # pre-fit residual (Gibbs) computed from most recent ST1 measurement
  ('st1_latest_postfit_res',    '>f4', (3,)), # post-fit residual (Gibbs) computed from most recent ST1 measurement
  ('st2_latest_innov',          '>f4', (6,)), # Innovation (delta Gibbs and delta bias) computed from most recent ST2 measurement
  ('st2_latest_prefit_res',     '>f4', (3,)), # pre-fit residual (Gibbs) computed from most recent ST2 measurement
  ('st2_latest_postfit_res',    '>f4', (3,)), # post-fit residual (Gibbs) computed from most recent ST2 measurement
  ('st3_latest_innov',          '>f4', (6,)), # Innovation (delta Gibbs and delta bias) computed from most recent ST3 measurement
  ('st3_latest_prefit_res',     '>f4', (3,)), # pre-fit residual (Gibbs) computed from most recent ST3 measurement
  ('st3_latest_postfit_res',    '>f4', (3,)), # post-fit residual (Gibbs) computed from most recent ST3 measurement

  ('st_err_bnd_rjct_consec_ctr','|u2'), # number consecutive ST meas rejected due to ExceedsExpectErrorBnd
  ('mekf_operation',            '|u1'), # Enable filter operation
  ('enable_ST3',                '|u1'), # ST3 processing enabled flag
  ('gyro_bias_est_enable',      '|u1'), # gyro bias estimation enabled flag
  ('aif_st1',                   '|u1'), # ST1 Auto/Inhibit/Force flag
  ('aif_st2',                   '|u1'), # ST2 Auto/Inhibit/Force flag
  ('aif_st3',                   '|u1'), # ST3 Auto/Inhibit/Force flag

  ('packed_bits',               '|u1', (3,)),  # to unpack later
  ('filter_is_init',            '|u1'), # filter is initialized flag
  ('gyro_buff_is_init',         '|u1'), # GyroBuffer is initialized flag
  ('st1_latest_reject_code',    '|u1'), # most recent ST1 reject code
  ('st2_latest_reject_code',    '|u1'), # most recent ST2 reject code
  ('st3_latest_reject_code',    '|u1'), # most recent ST3 reject code

  ('st1_mekf_disagree',         '|u1'), # ST1 meas disagree with MEKF est
  ('st2_mekf_disagree',         '|u1'), # ST2 meas disagree with MEKF est
  ('att_corr_excessive',        '|u1'), # 1 if attitude correction excessive
  ('gyro_bias_corr_excessive',  '|u1'), # 1 if gyro bias correction excessive
  ('gyro_bias_mag_excessive',   '|u1'), # 1 if gyro bias magnitude excessive
  ('mekf_err_code',             '|u1'), # General purpose MEKF error code
  ('gyro_meas_missing_ctr',     '|u1'), # number of gyro meas missing from GyroBuffer
  ('st1_accept_ctr',            '|u1'), # number of ST1 measurements that have been accepted

  ('st2_accept_ctr',            '|u1'), # number of ST2 measurements that have been accepted
  ('st3_accept_ctr',            '|u1'), # number of ST3 measurements that have been accepted
  ('st1_reject_ctr',            '|u1'), # number of ST1 measurements that have been rejected
  ('st2_reject_ctr',            '|u1'), # number of ST2 measurements that have been rejected
  ('st3_reject_ctr',            '|u1'), # number of ST3 measurements that have been rejected
  ('smugness_reinit_ctr',       '|u1'), # number of reinits due to smugness check
  ('gyro_buff_idx_new',         '|u1'), # GyroBuffer most recent entry index
  ('gyro_buff_idx_old',         '|u1'), # GyroBuffer oldest entry index

]) # total length = 888 bytes with primary packet header; 882 bytes without.


def APID108(data):
  apid = 108
  check_packet_size(apid, Fields.itemsize, len(data))

  tmp = np.frombuffer(data, dtype=Fields, count=1)
  myDict = getDict(tmp)
  myDict['timestamp'] = parse_CCSDS_timestamp(myDict['timestamp'])

  # unpack bit fields
  (
    myDict['excessive_gyro_meas_missing'], # 1 if excessive gyro missing
    myDict['invalid_covar_eigenvalues'],   # 1 if covariance eigenvalues too small
    myDict['excessive_covar_condition'],   # 1 if covariance condition number too big
    spare,                                 # ignore
    myDict['filter_initialized'],          # 1 if filter not init
    myDict['gyro_buffer_initialized'],     # 1 if gyro buffer not init
    myDict['excessive_att_uncert'],        # 1 if attitude estimate uncertainty excessive
    myDict['excessive_bias_uncert'],       # 1 if gyro bias estimate uncertainty excessive
  ) = np.unpackbits(myDict['packed_bits'][0])

  (
    myDict['att_X_uncert'],      # 1 if att X estimate uncertainty excessive
    myDict['att_Y_uncert'],      # 1 if att Y estimate uncertainty excessive
    myDict['att_Z_uncert'],      # 1 if att Z estimate uncertainty excessive
  ) = np.unpackbits(myDict['packed_bits'][1])[:3] # ignore last 5 bits

  (
    myDict['bias_X_uncert'],     # 1 if gyro X bias estimate uncertainty excessive
    myDict['bias_Y_uncert'],     # 1 if gyro Y bias estimate uncertainty excessive
    myDict['bias_Z_uncert'],     # 1 if gyro Z bias estimate uncertainty excessive
  ) = np.unpackbits(myDict['packed_bits'][2])[:3] # ignore last 5 bits

  del myDict['packed_bits']

  return myDict
