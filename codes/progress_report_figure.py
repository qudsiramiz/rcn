from cProfile import label
import pyspedas as spd
import pytplot as ptt
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pytz

read_data = ""

if read_data:
    trange = ["2015-09-11 15:20:00", "2015-09-11 15:24:00"]

    _ = spd.mms.fgm(trange=trange, probe=3, time_clip=True, level='l2', data_rate="brst",
                    get_fgm_ephemeris=True)

    _ = spd.mms.fpi(trange=trange, probe=3, time_clip=True, level='l2', data_rate="brst",
                    latest_version=True)

b_gsm = ptt.get_data("mms3_fgm_b_gsm_brst_l2_bvec")[1:][0]
b_gsm_time = ptt.get_data("mms3_fgm_b_gsm_brst_l2")[0]

# Convert from UNIX time to datetime (UTC)
b_gsm_time = [datetime.datetime.utcfromtimestamp(t) for t in b_gsm_time]
b_gsm_time = [t.replace(tzinfo=pytz.UTC) for t in b_gsm_time]
b_gsm_time = np.array(b_gsm_time)

n_p = ptt.get_data("mms3_dis_numberdensity_brst")[1]
n_p_time = ptt.get_data("mms3_dis_numberdensity_brst")[0]

# Convert time from UNIX time to date time
n_p_time = [datetime.datetime.fromtimestamp(t) for t in n_p_time]
n_p_time = [t.replace(tzinfo=pytz.UTC) for t in n_p_time]
n_p_time = np.array(n_p_time)
# Shift n_p_time by 5 hours
n_p_time = n_p_time + datetime.timedelta(hours=4)

vp_gse = ptt.get_data("mms3_dis_bulkv_gse_brst")[1:4][0]
vp_gse_time = ptt.get_data("mms3_dis_bulkv_gse_brst")[0]

# Convert time from UNIX time to date time
vp_gse_time = [datetime.datetime.fromtimestamp(t) for t in vp_gse_time]
vp_gse_time = [t.replace(tzinfo=pytz.UTC) for t in vp_gse_time]
vp_gse_time = np.array(vp_gse_time)
# Shift vp_gse_time by 5 hours
vp_gse_time = vp_gse_time + datetime.timedelta(hours=4)

# Covert gse to gsm
_ = spd.cotrans(name_in=f'mms3_dis_bulkv_gse_brst',
                name_out=f'mms3_dis_bulkv_gsm_brst', coord_in='gse',
                coord_out='gsm')

vp_gsm = ptt.get_data(f'mms3_dis_bulkv_gsm_brst')[1:4][0]

# Stndardize the vp_gsm z-component data
vp_gsm_z = vp_gsm.T[2, :]
vp_gsm_z_std = (vp_gsm_z - vp_gsm_z.mean())

# Set the fontstyle to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

fontsize = 14
alpha = 0.1
plt.close("all")

fig = plt.figure(num=None, figsize=(6, 4), dpi=200, facecolor='w', edgecolor='k')
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0, hspace=0)

# Magentic field plot
axs1 = fig.add_subplot(4, 1, 1)
axs1.plot(b_gsm_time, b_gsm.T[0, :], c="b", ls="-", lw=.5, label=r"$B_{\rm x}$")
axs1.plot(b_gsm_time, b_gsm.T[1, :], c="g", ls="-", lw=.5, label=r"$B_{\rm y}$")
axs1.plot(b_gsm_time, b_gsm.T[2, :], c="r", ls="-", lw=.5, label=r"$B_{\rm z}$")
axs1.axvspan(vp_gse_time[260], vp_gse_time[370], color="r", alpha=alpha)
axs1.axvspan(vp_gse_time[-180], vp_gse_time[-50], color="r", alpha=alpha)
#axs1.set_xlabel(r"Time (UTC)")
axs1.set_ylabel(r"(nT)", fontsize=fontsize)
axs1.legend(loc="upper left")
axs1.set_xlim(b_gsm_time[0], b_gsm_time[-1])
axs1.set_ylim(-30, 55)

# Proton number density plot
axs2 = fig.add_subplot(4, 1, 2, sharex=axs1)
axs2.plot(n_p_time, n_p, c="k", ls="-", lw=.5, label=r"$n_{\rm p}$")
axs2.axvspan(vp_gse_time[260], vp_gse_time[370], color="r", alpha=alpha)
axs2.axvspan(vp_gse_time[-180], vp_gse_time[-50], color="r", alpha=alpha)
axs2.legend(loc="lower left")
axs2.set_ylim(4e-1, 1e2)
axs2.set_yscale("log")
axs2.set_ylabel(r"(cm$^{-3}$)", fontsize=fontsize)

# Bulk velocity plot
axs3 = fig.add_subplot(4, 1, 3, sharex=axs1)
axs3.plot(vp_gse_time, vp_gsm.T[0, :], c="b", ls="-", lw=.5, label=r"$v_{\rm x}$")
axs3.plot(vp_gse_time, vp_gsm.T[1, :], c="g", ls="-", lw=.5, label=r"$v_{\rm y}$")
axs3.plot(vp_gse_time, vp_gsm.T[2, :], c="r", ls="-", lw=.5, label=r"$v_{\rm z}$")
axs3.axvspan(vp_gse_time[260], vp_gse_time[370], color="r", alpha=alpha)
axs3.axvspan(vp_gse_time[-180], vp_gse_time[-50], color="r", alpha=alpha)
axs3.legend(loc="upper left")
axs3.set_ylabel(r"(km/s)", fontsize=fontsize)
axs3.set_ylim(-380, 400)

# Bulk velocity plot (standardized)
axs4 = fig.add_subplot(4, 1, 4, sharex=axs1)
axs4.plot(vp_gse_time, vp_gsm_z_std, c="r", ls="-", lw=1,
          label=r"($V_{\rm p,z} - \left< V_{\rm p,z} \right>)$")
axs4.set_ylabel(r"(km/sec)", fontsize=fontsize)
axs4.axvspan(vp_gse_time[260], vp_gse_time[370], color="r", alpha=alpha, label=r"jet location")
axs4.axvspan(vp_gse_time[-180], vp_gse_time[-50], color="r", alpha=alpha)
axs4.legend(loc="upper left")
axs4.set_ylim(-250, 250)
axs4.set_xlabel(r"Time (UTC)", fontsize=fontsize)

plt.savefig("mms3_fgm_bulkv_gsm_brst.png", dpi=200, bbox_inches="tight", pad_inches=0.1)