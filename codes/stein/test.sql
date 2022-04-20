#
#  Test SQL (Structured Query Language) to retreive data
#  from 2 tables in a databases:
#
#	table A = MMS2up - SC positions, MP crossing times, flags
#	table B = MMS2up_misc -HT analysis, Vn_MFR etc. 
#
#    A.DateStart,
#    format(B.Vn_MFR,2) as Velocity_MFR,
#    dot(A.x1x_mvab,A.x1y_mvab,A.x1z_mvab,  B.VHTi_tot_x,B.VHTi_tot_y,B.VHTi_tot_z) as
#    Vht_dot_n

SOURCE functions.sql;
SELECT
    A.EventId as EventId,
    A.Flagstr as FlagStr,
    format(A.X,2) as Xgse,
    format(A.Y,2) as Ygse,
    format(A.Z,2) as Zgse

FROM
    MMS2up A, MMS2up_misc B 
WHERE
    A.EventId = B.EventId AND
    A.FlagStr LIKE "mp%" AND
    A.DateStart BETWEEN "2018-04-01%" AND "2018-04-05%"
;


