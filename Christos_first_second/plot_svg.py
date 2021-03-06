from model_phantom_DB import *
from linares_plot import * 
import itertools

#Palettes
print(sns.color_palette("tab10").as_hex())
sns.palplot(sns.color_palette("tab10"))
plt.show(block=False)
plt.close()

c_on = 'darkorange' #'#ff7f0e'
c_off = 'dodgerblue' #'#1f77b4'

onoff_pal=[c_on, c_off]
offon_pal = [c_off, c_on]

pal_cyan = sns.color_palette("RdBu_r", n_colors=200)[40:] #RdBu_r
ltc= 'gold'  #'springgreen'
lw_t=3
N=512
stimon = 50
stimoff = 50 + floor(250/2) ;
nsteps=int(floor(750));
p_targ2 = int((N * 90)/360)
p_dist = int((N * (90+57) )/360)
p_dist2 = int((N * (90+170) )/360)
sns.set_context("poster", font_scale=1.1)
sns.set_style("ticks")


def hemap(an):
    RE_sorted=flipud(an[4])
    sns.heatmap(RE_sorted, cmap=pal_cyan, vmin=0, vmax=17,  cbar=True)
    plt.gca().set_ylabel('')
    plt.gca().set_xlabel('')
    plt.gca().set_title('')
    plt.gca().plot([stimon, stimon+400], [p_targ2, p_targ2], ls='--', color =ltc, linewidth=lw_t) ## flipped, so it is p_target 
    #plt.gca().set_xticks([])
    plt.gca().set_xticks([0,an[4].shape[1]/2, an[4].shape[1]])
    plt.gca().set_xticklabels(['0', str(an[4].shape[1]), str(an[4].shape[1]*2)], rotation=0)
    #axn.set_xticks([0,simul[4].shape[1]/2, simul[4].shape[1]])
    plt.gca().set_xticklabels(['0s', str(an[4].shape[1]/1000) + 's', str( int(an[4].shape[1]*2/1000) ) + 's'], rotation=0)
    plt.gca().set_yticks([0, N/4, N/2,  3*N/4, N ])
    #plt.gca().set_yticklabels(['0','90','180', '270', '360'])
    plt.gca().set_yticklabels(['0','','$^\pi$', '', '2$^\pi$'])
    plt.gca().set_xlabel('time (ms)');
    plt.gca().set_ylabel('neuron ($^\circ$)');
    
##
def hemap2(an):
    RE_sorted=flipud(an[4])
    sns.heatmap(RE_sorted, cmap=pal_cyan, vmin=0, vmax=17,  cbar=True)
    plt.gca().set_ylabel('')
    plt.gca().set_xlabel('')
    plt.gca().set_title('')
    plt.gca().plot([500, 900], [p_dist, p_dist], ls='--', color =ltc, linewidth=lw_t) ## flipped, so it is p_target 
    #plt.gca().set_xticks([])
    plt.gca().set_xticks([0,an[4].shape[1]/2, an[4].shape[1]])
    plt.gca().set_xticklabels(['0', str(an[4].shape[1]), str(an[4].shape[1]*2)], rotation=0)
    #axn.set_xticks([0,simul[4].shape[1]/2, simul[4].shape[1]])
    plt.gca().set_xticklabels(['0s', str(an[4].shape[1]/1000) + 's', str( int(an[4].shape[1]*2/1000) ) + 's'], rotation=0)
    plt.gca().set_yticks([0, N/4, N/2,  3*N/4, N ])
    #plt.gca().set_yticklabels(['0','90','180', '270', '360'])
    plt.gca().set_yticklabels(['0','','$^\pi$', '', '2$^\pi$'])
    plt.gca().set_xlabel('time (ms)');
    plt.gca().set_ylabel('neuron ($^\circ$)');
    
    
def hemap2f(an):
    RE_sorted=flipud(an[4])
    sns.heatmap(RE_sorted, cmap=pal_cyan, vmin=0, vmax=17,  cbar=True)
    plt.gca().set_ylabel('')
    plt.gca().set_xlabel('')
    plt.gca().set_title('')
    plt.gca().plot([500, 900], [p_dist2, p_dist2], ls='--', color =ltc, linewidth=lw_t) ## flipped, so it is p_target 
    #plt.gca().set_xticks([])
    plt.gca().set_xticks([0,an[4].shape[1]/2, an[4].shape[1]])
    plt.gca().set_xticklabels(['0', str(an[4].shape[1]), str(an[4].shape[1]*2)], rotation=0)
    #axn.set_xticks([0,simul[4].shape[1]/2, simul[4].shape[1]])
    plt.gca().set_xticklabels(['0s', str(an[4].shape[1]/1000) + 's', str( int(an[4].shape[1]*2/1000) ) + 's'], rotation=0)
    plt.gca().set_yticks([0, N/4, N/2,  3*N/4, N ])
    #plt.gca().set_yticklabels(['0','90','180', '270', '360'])
    plt.gca().set_yticklabels(['0','','$^\pi$', '', '2$^\pi$'])
    plt.gca().set_xlabel('time (ms)');
    plt.gca().set_ylabel('neuron ($^\circ$)');



def hemap_p(an):
    RE_sorted=flipud(an[4])
    sns.heatmap(RE_sorted, cmap=pal_cyan, vmin=0, vmax=17,  cbar=True)
    plt.gca().set_ylabel('')
    plt.gca().set_xlabel('')
    plt.gca().set_title('')
    plt.gca().plot([500, 900], [p_targ2, p_targ2], ls='--', color =ltc, linewidth=lw_t) ## flipped, so it is p_target 
    #plt.gca().set_xticks([])
    plt.gca().set_xticks([0,an[4].shape[1]/2, an[4].shape[1]])
    plt.gca().set_xticklabels(['0', str(an[4].shape[1]), str(an[4].shape[1]*2)], rotation=0)
    #axn.set_xticks([0,simul[4].shape[1]/2, simul[4].shape[1]])
    plt.gca().set_xticklabels(['0s', str(an[4].shape[1]/1000) + 's', str( int(an[4].shape[1]*2/1000) ) + 's'], rotation=0)
    plt.gca().set_yticks([0, N/4, N/2,  3*N/4, N ])
    #plt.gca().set_yticklabels(['0','90','180', '270', '360'])
    plt.gca().set_yticklabels(['0','','$^\pi$', '', '2$^\pi$'])
    plt.gca().set_xlabel('time (ms)');
    plt.gca().set_ylabel('neuron ($^\circ$)');

#


fee=1
fei=1
fie=1
fii=1


# on_close_1= model(totalTime=3000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=57, tauE=20, tauI=10,  n_stims=2, I0E=0.05, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fii,
#            GEI=0.13*fei,
#            GIE=0.042*fie, 
#            sigE=7., sigI=5., k_noise=0.6,           
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=1.2, phantom_onset=50000, phnatom_duration=100)

# hemap(on_close_1)
# plt.show(block=False)
# #plt.savefig("1_o2_ON.svg")

# off_close_1= model(totalTime=3000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=57, tauE=20, tauI=10,  n_stims=2, I0E=-3.5, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fii,
#            GEI=0.13*fei,
#            GIE=0.042*fie, 
#            sigE=7., sigI=5., k_noise=0.6,           
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=2., phantom_onset=50000, phnatom_duration=100)

# hemap(off_close_1)
# plt.show(block=False)
# plt.savefig("1_close_off.pdf")


# on_far_1= model(totalTime= 3000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=0.05, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fii,
#            GEI=0.13*fei,
#            GIE=0.042*fie, 
#            sigE=7., sigI=5., k_noise=0.6,            
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=1.2, phantom_onset=50000, phnatom_duration=100)

# hemap(on_far_1)
# plt.show(block=False)
# plt.savefig("1_far_on_small_drift.pdf")



# off_far_1= model(totalTime= 3000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=-3.5, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fii,
#            GEI=0.13*fei,
#            GIE=0.042*fie, 
#            sigE=7., sigI=5., k_noise=0.6,            
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=1.2, phantom_onset=50000, phnatom_duration=100)

# hemap(off_far_1)
# plt.show(block=False)
# plt.savefig("1_far_off.pdf")




fee=0.94
fei=0.92
fie=1.14
fii=1.08


# on_close_2= model(totalTime=3000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=57, tauE=20, tauI=10,  n_stims=2, I0E=0.05, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.13*fie,
#            GIE=0.042*fii, 
#            sigE=7., sigI=5., k_noise=0.6,             
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=1.2, phantom_onset=50000, phnatom_duration=100)

# hemap2(on_close_2)
# plt.show(block=False)
# plt.savefig("2_close_on.pdf")



# off_close_2= model(totalTime=3000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=57, tauE=20, tauI=10,  n_stims=2, I0E=-2., I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.13*fie,
#            GIE=0.042*fii, 
#            sigE=7., sigI=5., k_noise=0.6,              
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=1.2, phantom_onset=50000, phnatom_duration=100)

# hemap2(off_close_2)
# plt.show(block=False)
# plt.savefig("2_close_off.pdf")


# on_far_2= model(totalTime=3000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=0.05, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.13*fie,
#            GIE=0.042*fii, 
#            sigE=7., sigI=5., k_noise=0.6,             
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=1.2, phantom_onset=50000, phnatom_duration=100)

# hemap2f(on_far_2)
# plt.show(block=False)
# plt.savefig("2_far_on.pdf")



off_far_2= model(totalTime=3000, targ_onset_1=100, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
           angle_separation=170, tauE=20, tauI=10,  n_stims=2, I0E=-2., I0I=0.5, 
           GEE=0.068*fee,
           GII= 0.13*fei,
           GEI=0.13*fie,
           GIE=0.042*fii, 
           sigE=7., sigI=5., k_noise=0.6,             
           kappa_E=45, 
           kappa_I=0.3, 
           kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
           plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
           phantom_st=1.2, phantom_onset=50000, phnatom_duration=100)

hemap2f(off_far_2)
plt.show(block=False)
# plt.savefig("2_far_off_drift.pdf") ##


# fee=0.94
# fei=0.92
# fie=1.14
# fii=1.08

# phantom_on= model(totalTime=2000, targ_onset_1=10000, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=0, tauE=20, tauI=10,  n_stims=2, I0E=0.05, I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.13*fie,
#            GIE=0.042*fii, 
#            sigE=7., sigI=5., k_noise=0.6,             
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=1.2, phantom_onset=100, phantom_on = 'on', phnatom_duration=300)

# hemap_p(phantom_on)
# plt.show(block=False)
# plt.savefig("phantom_on.pdf")

# phantom_off= model(totalTime=2000, targ_onset_1=10000, targ_onset_2=1000, angle_target_i=90, presentation_period=100,
#            angle_separation=0, tauE=20, tauI=10,  n_stims=2, I0E=-2., I0I=0.5, 
#            GEE=0.068*fee,
#            GII= 0.13*fei,
#            GEI=0.13*fie,
#            GIE=0.042*fii, 
#            sigE=7., sigI=5., k_noise=0.6,             
#            kappa_E=45, 
#            kappa_I=0.3, 
#            kappa_stim=40., N=512, stim_strengthE=9.4, stim_strengthI=0.,
#            plot_connectivity=False, plot_rate=False, plot_hm=False , plot_fit=False, 
#            phantom_st=1.2, phantom_onset=100, phantom_on = 'on', phnatom_duration=500)

# hemap_p(phantom_off)
# plt.show(block=False)
# plt.savefig("phantom_off.pdf")