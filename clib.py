import numpy as np
import struct
import sys
import os
import matplotlib.pyplot as plt
from itertools import combinations

# This code was originally written by Vikram Ravi.
# It has been since modified by multiple people, including Ewan Barr, Chris
# Flynn, Stefan Oslowski, Jamie Tsai, and Morgan Oneill
class CorrPlot(object):
    def __init__(self, pairs, antfile, best_idxs, good_idxs, mods, args, batch_size=1024, ntop=4, nukeF1=False, nramps=20, wcut=0.2):
        self.pairs = pairs
        self.best_idxs = best_idxs
        with open(antfile, "r") as fp:
            ants = [str.split(x) for x in fp.readlines()]
            ant_names = [ant_info[0] for ant_info in ants]
            ant_dists = [-float(ant[1]) if ant[0][0] == 'W' else float(ant[1]) for ant in ants]
            self.antenna_sorter = np.array(ant_dists).argsort()
            self.ants = np.array(ant_names)[self.antenna_sorter]


        self.ants_unsort = np.loadtxt(args.antfile,dtype='str').transpose()[0]

        self.batch_size = batch_size
        self.ntop = ntop
        self.nukeF1 = nukeF1
        self.nramps = nramps
        self.wcut = wcut
        self.fig = plt.figure(figsize=[14,10])
        self.map_ax_E = self.fig.add_subplot(3, 1, 1)
        self.map_ax_W = self.fig.add_subplot(3, 1, 2)
        self.cp_spec_ax = self.fig.add_subplot(3, 4, 9)
        self.phase_ax = self.fig.add_subplot(3, 4, 10)
        self.corr_ax = self.fig.add_subplot(3, 4, 11)
        self.fraction_ax = self.fig.add_subplot(3, 4, 12)
        self.good_idxs = good_idxs
        self.mods=mods
        self.args=args
        self.plot()

    def plot(self):
        nbins = self.pairs["b"].max()+1
        map_ar_EW = [] # Plot East and West separately
        map_ar_EW.append(np.zeros([self.ntop, nbins/2])) # West
        map_ar_EW.append(np.zeros([self.ntop, nbins/2])) # East


        for a,b,delay,sn,ar in self.pairs:
            sn = abs(ar).max()
            for i in xrange(0,self.ntop):
                if b == self.best_idxs[i]:
                    for j in range(0, nbins):
                        if j<=nbins/2-1:
                            if self.ants_unsort[a] == self.ants[j]:
                                map_ar_EW[0][i][j] = sn
                        elif self.ants_unsort[a] == self.ants[j]:
                            map_ar_EW[1][i][j-nbins/2] = sn
                elif a == self.best_idxs[i]:
                    for j in range(0, nbins):
                        if j<=nbins/2-1:
                            if self.ants_unsort[b] == self.ants[j]:
                                map_ar_EW[0][i][j] = sn
                        elif self.ants_unsort[b] == self.ants[j]:
                            map_ar_EW[1][i][j-nbins/2] = sn

        logsnr_W = np.log10(map_ar_EW[0])
        logsnr_E = np.log10(map_ar_EW[1])

        min_imshow, max_imshow = -1., np.log10(map_ar_EW).max()
        print min_imshow, max_imshow
        im = self.map_ax_W.imshow(logsnr_W,aspect="auto",
                                picker=True,
                                interpolation="nearest",cmap="jet", origin='low',vmin = min_imshow, vmax = max_imshow)
        plt.colorbar(im, ax=self.map_ax_W)
        self.map_ax_W.yaxis.set_label_position("right")
        self.map_ax_W.set_xticks(range(0, nbins/2, 16)) # set ticks at every 16th point, i.e. one per 4 modules
        self.map_ax_W.set_xticklabels(self.ants[range(0, nbins/2+1, 16)])
        self.map_ax_W.set_yticks(range(self.ntop))
        self.map_ax_W.set_yticklabels([self.ants[np.where(self.antenna_sorter==i)[0][0]] for i in self.best_idxs])
        plt.ylabel("vs West")


        im = self.map_ax_E.imshow(logsnr_E,aspect="auto",
                                picker=True,
                                interpolation="nearest",cmap="jet", origin='low',vmin = min_imshow, vmax = max_imshow)
        plt.colorbar(im, ax=self.map_ax_E)
        self.map_ax_E.yaxis.set_label_position("right")
        self.map_ax_E.set_xticks(range(0, nbins/2, 16)) # set ticks at every 16th point, i.e. one per 4 modules
        self.map_ax_E.set_xticklabels(self.ants[range(nbins/2, nbins, 16)])
        self.map_ax_E.set_yticks(range(self.ntop))
        self.map_ax_E.set_yticklabels([self.ants[np.where(self.antenna_sorter==i)[0][0]] for i in self.best_idxs])
        plt.ylabel("vs East")


        self.fig.canvas.mpl_connect("pick_event",self.onpick)
        
    def onpick(self,event):

# CF 09/12/16, added 0.5 to the event locations to get the clicks correctly centered
        x = int(event.mouseevent.xdata+0.5)
        y = int(event.mouseevent.ydata+0.5)

        # Figure out if we're clicking on reference vs West or vs East:
        west_antennas = False 
        # checking for east would be: event.mouseevent.inaxes == self.map_ax_E
        if event.mouseevent.inaxes == self.map_ax_W:
            west_antennas = True

        if x < 0:
            x = 0
        if y < 0:
            y = 0

        #print "Picked pair: %d -> %d"%(x,y)
        
        if not west_antennas:
            x += ( self.pairs["b"].max() + 1. ) / 2

        y = self.best_idxs[int(y)];
        ref_id = int(y)
        other_id = int(x)
        swapped=False
        if x < y:
            x, y = y, x
            swapped=True
        #print "Picked pair: %d -> %d"%(x,y)

        calib_delays_foo = open('calib.delays','r')
        antenna_distance_foo = open('../obs.antenna','r')
        modsdel_baselined = np.load('modsdel_baselined.npy','r')
        modsph_baselined  = np.load('modsph_baselined.npy','r')
# fix to printed out antenna name (CF: 09/12/16)
        if swapped:
            for i in calib_delays_foo:
                if self.ants[int(y)] == i[:5]:
                    print "Antenna : %s -> %s"%(self.ants[np.where(self.antenna_sorter==int(x))[0][0]], # [0][0] as np.where returns a tuple of arrays
                    self.ants[int(y)]), map(float, i.split(" ")[1:])[0]*10**9,"(ns), ",map(float, i.split(" ")[1:])[1], "(rad), ", map(float, i.split(" ")[1:])[2], "(weight)"
                    if self.args.verbose:
                        for j in range(len(self.good_idxs)):
                            if np.loadtxt(args.antfile,dtype='str').transpose()[0][self.good_idxs[j]] == self.ants[int(y)]:
                                print self.ants[int(y)], modsdel_baselined[j]*self.args.tsamp*1e3, '(ns)', modsph_baselined[j], '(rad)'


            for i in antenna_distance_foo:
                if i[:5] == self.ants[int(y)]:
                    location_a = float(i.split(" ")[1])
                    if i[0] =="W": location_a *= -1.
                if i[:5] == self.ants[np.where(self.antenna_sorter==int(x))[0][0]]:
                    location_b = float(i.split(" ")[1])
                    if i[0] =="W": location_b *= -1.
            #print location_a, location_b
            if abs( location_a - location_b )<50: print "short baseline, no contribution"

            
            if self.ants[np.where(self.antenna_sorter==int(x))[0][0]] == self.ants[int(y)]:
                print "same antennas!, returned"
                return
        else:
            if self.ants[np.where(self.antenna_sorter==int(y))[0][0]] == self.ants[int(x)]:
                print "same antennas!, returned"
                return
            for i in calib_delays_foo:
                if self.ants[int(x)] == i[:5]:
                    print "Antenna : %s -> %s"%(self.ants[np.where(self.antenna_sorter==int(y))[0][0]],
                    self.ants[int(x)]), map(float, i.split(" ")[1:])[0]*10**9,"(ns), ",map(float, i.split(" ")[1:])[1], "(rad), ", map(float, i.split(" ")[1:])[2], "(weight)"
                    if self.args.verbose:
                        for j in range(len(self.good_idxs)):
                            if np.loadtxt(args.antfile,dtype='str').transpose()[0][self.good_idxs[j]] == self.ants[int(x)]:
                                print self.ants[int(x)], -1*modsdel_baselined[j]*self.args.tsamp*1e3, '(ns)', -1*modsph_baselined[j], '(rad)'

            for i in antenna_distance_foo:
                if i[:5] == self.ants[int(x)]:
                    location_a = float(i.split(" ")[1])
                    if i[0] =="W": location_a *= -1.
                if i[:5] == self.ants[np.where(self.antenna_sorter==int(y))[0][0]]:
                    location_b = float(i.split(" ")[1])
                    if i[0] =="W": location_b *= -1.
            #print location_a, location_b
            if abs( location_a - location_b )<50: print "short baseline, no contribution"




        self.cp_spec_ax.cla()
        self.phase_ax.cla()
        self.corr_ax.cla()
        self.fraction_ax.cla()
        
        if swapped:
            idxs, = np.where((self.pairs["a"] == self.antenna_sorter[int(y)]) & (self.pairs["b"] == int(x)) )
        else:
            idxs, = np.where((self.pairs["b"] == self.antenna_sorter[int(x)]) & (self.pairs["a"] == int(y)) )
        if idxs.size == 0:
            idxs, = np.where((self.pairs["a"] == y) & (self.pairs["b"] == x) )

        pair = self.pairs[idxs[0]]
        ar = pair["ar"]

        # brute force cleanup of Telstra F1 band
        cp = np.fft.fft(ar)
        cpclean = cp
        f1start = 640
        f1end = 810
        if self.nukeF1:
            cpclean = np.fft.fft(ar)
            for i in range(len(cpclean)):
                if i>=f1start and i<=f1end :
                    cpclean[i] /= 1.0e8

        arclean = np.fft.ifft(cpclean)        
        cpclean = np.fft.fft(arclean)

        # Plot lag spectrum
        self.corr_ax.plot(abs(np.fft.fftshift(ar)),'g-')
        self.corr_ax.plot(abs(np.fft.fftshift(arclean)),'b-')
        self.corr_ax.set_xlim(0, self.batch_size)
        self.corr_ax.set_xlabel('Lag Spectrum')

        delay = np.linspace(-0.5,0.5,self.nramps)
        self.fraction_ax.plot(delay, abs(np.fft.ifft(delay_signal2(cpclean, delay))).max( axis=1) ,'rD--',markerfacecolor='none')
        self.fraction_ax.plot(delay, abs(            delay_signal2(cpclean, delay  ).mean(axis=1)),'bo--',markerfacecolor='none')
        self.fraction_ax.grid()
        self.fraction_ax.set_xlabel('Noisy Lag Checker')
        #There is a complex frequency bandpass data(phase temporal corrected).
        #The max of abs of ifft of this data is rD-
        #The abs of mean of this data is bo-
        
        npt = len(cp)

        # set the DC channel to something reasonable
        start = npt / 8 
        end = npt - start
        centre_mon = cp[start:end].min()
        cp[npt/2] = centre_mon

        # Plot cross power spectrum
        abs_cp = abs(cp)
        abs_cpclean = abs(cpclean)
        self.cp_spec_ax.scatter(range(abs_cp.size),np.log10(abs_cp),s=1,edgecolor="g")
        self.cp_spec_ax.scatter(range(abs_cpclean.size),np.log10(abs_cpclean),s=1,edgecolor="b")
        ymax = np.log10(abs_cpclean.max())
        ymin = np.log10(abs_cpclean.min())
        yrange = ymax - ymin
        ypad = yrange * 0.05;
        self.cp_spec_ax.set_ylim(ymin - ypad, ymax + ypad)
        self.cp_spec_ax.set_ylim(bottom=-2.0)
        self.cp_spec_ax.set_xlim(0,self.batch_size)
        self.cp_spec_ax.set_xlabel('Cross Power Spectrum')


        # Plot cross power phases
        cp_phase = np.angle(cpclean)
        self.phase_ax.scatter(range(cp_phase.size),cp_phase,s=1,edgecolor="g")
        for i in range(len(cp_phase)):
            if i>=f1start and i<=f1end :
                cp_phase[i] = -99.0  # null value
        self.phase_ax.scatter(range(cp_phase.size),cp_phase,s=1,edgecolor="b")
        self.phase_ax.set_ylim(-np.pi,np.pi)
        self.phase_ax.set_xlim(0, self.batch_size)
        self.phase_ax.set_xlabel('Residual Cross Power Phases')

        plt.draw()

def stupid_snr(ar,batch_size):
    baseline = np.hstack((ar[:batch_size/3],ar[2*batch_size/3:]))
    ar -= baseline.mean()
    return ar.max() / baseline.std()

def new_snr(row,batch_size,args):

    nch = batch_size
    nramps = args.nramps # over 10 samples
    wcut = args.wcut
    x = np.arange(nch)
    n_excise = 10

    # get bad fine channels
    powers = (np.abs(row)**2.)
    spowers = np.sort(powers)
    slocs = np.argsort(powers)
    steps = np.diff(spowers) #spowers[1:nch]-spowers[0:nch-1] # <- np.diff
    step_fracs = np.zeros(nch)
    step_fracs[1:nch] = steps/np.max(np.sort(steps)[0:nch/2]) 
    pm1 = np.zeros(nch)+1.
    pm1[slocs[np.where(step_fracs > 100.0)]] = 0.0

    # test each coarse channel
    edge_mask = np.zeros((args.nchan,nch/args.nchan))+1.
    edge_mask[:,0:3]=0.
    edge_mask[:,nch/args.nchan-3:nch/args.nchan]=0.
    power_stats = np.reshape(powers*edge_mask.reshape(nch),(args.nchan,nch/args.nchan)).mean(axis=1)
    locs = np.argsort(power_stats)

    spower_stats = np.sort(power_stats)
    steps = spower_stats[1:args.nchan] - spower_stats[0:args.nchan-1]
    cutoff=-1
    i = args.nchan-1-n_excise
    while cutoff==-1 and i<args.nchan-1:
        if steps[i]>4.*np.mean(steps[0:i]):
            cutoff=i
        i += 1
    power_mask = np.zeros(args.nchan)+1.
    #if cutoff > 0:
    #     power_mask[cutoff+1:args.nchan]=0.0
    pm2 = np.zeros(args.nchan)+1.
    pm2[locs[power_mask==0.0]] = 0.0
    pm2 = pm1*np.tile(pm2,(nch/args.nchan,1)).transpose().reshape(nch)

    ar = np.fft.ifft(row*pm2)
    ar = np.fft.fftshift(ar)

    # coarse delay
    dely_coarse = (x[np.log10(np.abs(ar)) == np.log10(np.abs(ar)).max()]*1.-nch/2.)*2.*np.pi
    coarse_corr = row*(np.cos(x*(dely_coarse)/(1.*nch))+1j*np.sin(x*(dely_coarse)/(1.*nch)))
    
    #coarse delay re-write:
    #dely_coarse = np.abs(ar).argmax() - nch/2.
    #ramp = np.e**(np.pi*2*1j*np.linspace(0,1,size,endpoint=False) * dely_coarse)
    #coarse_corr = row * ramp

    # delay ramps and fit
    delys = (np.arange(nramps)+0.5)*4.*np.pi/(1.*nramps)-2.*np.pi
    stds = np.zeros(nramps)
    for i in range(nramps):
        # Ewan: This is incorrect as it does not deal with phases near pi or -pi where there is a wrap
        stds[i] = np.std((np.sin(np.angle(coarse_corr*(np.cos(x*(delys[i])/(1.*nch))+1j*np.sin(x*(delys[i])/(1.*nch))))))[pm2 != 0.0])**2.
    good_dely = delys[stds == np.min(stds)]+dely_coarse

    # dat_good = decoarse and define delayed signal
    dat_good = row*(np.cos(x*(good_dely)/(1.*nch))+1j*np.sin(x*(good_dely)/(1.*nch)))
    
    mn = np.mean(dat_good)
    dat_good *= np.conj(mn/np.abs(mn)) # Zero the phase of the visibility
    avg_dgood = np.sum(np.reshape(np.angle(dat_good)*pm2,(nch/4,4)),1)/4. #eh? this is just the mean of the phase on a window of 4.
    # avg_dgood = np.angle(dat_good).reshape(nch/4,4).mean(axis=1)
    x2 = np.sum(np.reshape(x,(nch/4,4)),1)/4.
    # x2 = x.reshape(nch/4,4).mean(axis=1)

    # Now we move to fitting a version of the data that has been binned by a factor of 4.
    vals = np.polyfit(x2[avg_dgood!=0.0],avg_dgood[avg_dgood!=0.0],1)
    
    # This line appears to be pointless
    # f = np.poly1d(vals)

    # this is a secondary fine delay correction based on a 1d polyfit
    fin_dely = good_dely/(2.*np.pi)-vals[0]*nch/(2.*np.pi) # in samples

    # this is the final phase correction
    # this is not true phase however as the phases have been normalized about 0 prior to the poly fit
    fin_phase = np.angle(np.mean(row[pm2!=0.0]*(np.cos((x[pm2!=0.0]-nch/2.)*(fin_dely*2.*np.pi)/(1.*nch))+1j*np.sin((x[pm2!=0.0]-nch/2.)*(fin_dely*2.*np.pi)/(1.*nch)))))
    
    # Here we do a final derotation of phase
    fin_row = pm2*row*(np.cos(-fin_phase+(x-nch/2.)*(fin_dely*2.*np.pi)/(1.*nch))+1j*np.sin(-fin_phase+(x-nch/2.)*(fin_dely*2.*np.pi)/(1.*nch)))
    #fin_row = pm2*row*(np.cos(-fin_phase+x*(fin_dely*2.*np.pi)/(1.*nch))+1j*np.sin(-fin_phase+x*(fin_dely*2.*np.pi)/(1.*nch)))

    # Given that nobody seems to be using the args.nchan argument I am bemused by this
    # nch/args.nchan will default to 1280/20 so this will be 20 rows by 64 columns
    # Assuming that we actually used the correct args.nchans, this would be the MEAN of the real for each coarse channel/std of the imaginary
    fsnr = (np.sum(np.reshape(np.real(fin_row),(args.nchan,nch/args.nchan)),1)/(nch/args.nchan))/np.std(np.reshape(np.imag(fin_row),(args.nchan,nch/args.nchan)),axis=1)
    #Only select fsnrs of coarse channels that have a non zero mean in their real component (its a mask)
    #Signal is now phased to zero, so component in phase is astro and out of phase (i.e. imag component) is noise
    fsnr = fsnr[(np.sum(np.reshape(np.real(fin_row),(args.nchan,nch/args.nchan)),1)/(nch/args.nchan))!= 0.0]
    
    #Eh? What? Sort all the S/N to find the best channel and then select the second best channel?????
    #Is this just a mistake or is there some kind of reasoning here?
    fsnr = np.sort(fsnr)[len(fsnr)-2]
    
    if fsnr < 1.0: # note (10/50000)*sqrt(1.5e7*300)*sqrt(1/320) = 0.75
        fin_snr = 1e8
    else:
        fin_snr = 1./fsnr
    # What is fin_snr? Is it a reciprocal weight?
    return (fin_row,fin_dely,-fin_phase,fin_snr)
    
def solve_snr(x12,x13,x23):

    return(x12*x13/x23,x12*x23/x13,x13*x23/x12)

def makePFB(ants):

    new_ants = ants.copy()

    for i in range(len(ants)):
        if ants[i][0]=='E':
            grp='EG'
        if ants[i][0]=='W':
            grp='WG'
        
        num = int(ants[i][1:3])
        col = ants[i][4]
        if col=='B':
            col=0
        if col=='G':
            col=1
        if col=='Y':
            col=2
        if col=='R':
            col=3
        pfbnum = (num-1)*4+col
        in_pfb = str(pfbnum % 16)
        pfb = str(1+int(np.floor(pfbnum/16))).zfill(2)
        
        new_ants[i] = grp+pfb+'_'+in_pfb

    return new_ants
        

def proc_ants(mods,bests,mod_status,good_idxs,best_idxs,args,sca_sn):

    nbest = best_idxs.size
    ngood = good_idxs.size
    nant = mod_status.size
    tsamp = args.tsamp*1e-6 # default in seconds

    # read in obs.antenna file
    ants = np.loadtxt(args.antfile,dtype='str')
    dists = (ants.transpose()[1]).astype('float')
    dists[nant/2:nant] *= -1.
    ants = ants.transpose()[0]
    pfb = makePFB(ants)

    # calculate baseline weights - zero for < 50m
    baseline = np.zeros((mods.size/nbest,nbest))
    size_by_nbest_range = range(mods.size/nbest)
    #for i in range(mods.size/nbest):
    for i in size_by_nbest_range:
        for j in np.arange(nbest):
            baseline[i,j] = np.abs(dists[good_idxs[i]]-dists[best_idxs[j]])
    #np.save('baseline0',baseline)
    baseline[baseline < 50.]=0.
    baseline[baseline >= 50.]=1.

    #np.save('baseline1',baseline)

    # delays

    del_corrs = bests["del"][0,:]
    del_corrs[0]=0.

    for i in size_by_nbest_range:
        for j in np.arange(nbest-1)+1:

            if np.abs(np.abs(mods["del"][i,0])-np.abs(mods["del"][i,j]+del_corrs[j])) < np.abs(np.abs(mods["del"][i,0])-np.abs(mods["del"][i,j]-del_corrs[j])):
                mods["del"][i,j] += del_corrs[j]
            else:
                mods["del"][i,j] -= del_corrs[j]

            mods["del"][i,j] = np.abs(mods["del"][i,j])*np.sign(mods["del"][i,0])

    #print "delays"
    #print 'mods["del"]', mods["del"].shape, mods["del"]
    good_dels = np.sum(mods["del"]*baseline,axis=1)/np.sum(baseline,axis=1)

    #for i in range(good_idxs.size):
    #    print ants[good_idxs[i]], good_dels[i]*-1*tsamp, mods["del"][i]*baseline*tsamp
    #    print 

    np.save('modsdel_baselined', mods["del"]*baseline) # mods["del"]*baseline)
    #np.save('del_corrs',del_corrs)
    #print 'del_corrs', del_corrs.shape, del_corrs
    #print 'good_dels',good_dels.shape, good_dels

    # phases

    ph_corrs = bests["ph"][0,:]
    ph_corrs[0]=0.

    for i in size_by_nbest_range:
        for j in np.arange(nbest-1)+1:
            if np.abs(np.abs(mods["ph"][i,0])-np.abs(mods["ph"][i,j]+ph_corrs[j])) < np.abs(np.abs(mods["ph"][i,0])-np.abs(mods["ph"][i,j]-ph_corrs[j])):
                mods["ph"][i,j] += ph_corrs[j]
            else:
                mods["ph"][i,j] -= ph_corrs[j]

            mods["ph"][i,j] = np.abs(mods["ph"][i,j])*np.sign(mods["ph"][i,0])

    #print "phases"
    #print 'mods["ph"]', mods["ph"].shape, mods["ph"]
    good_phs =  np.sum(mods["ph"]*baseline,axis=1)/np.sum(baseline,axis=1)

    #for i in range(good_idxs.size):
    #    print ants[good_idxs[i]], mods['ph'][i], good_phs[i]
    #print ants[good_idxs[1]], good_phs[1], (mods["ph"]*baseline).shape

    np.save('modsph_baselined', mods['ph']*baseline) #mods["ph"]*baseline)
    #print 'ph_corrs', ph_corrs.shape,ph_corrs
    #print 'good_phs', good_phs.shape,good_phs
    # SNRs

    good_snrs = np.zeros((mods.size/nbest,nbest*(nbest-1)/2))
    best_snrs = np.zeros((nbest,(nbest-1)*mods.size/nbest))
    best_snct = np.zeros(nbest)
    combi = [i for i in combinations(range(nbest),2)]

    for i in size_by_nbest_range:
        for j in range(nbest*(nbest-1)/2):

            b1 = combi[j][0]
            b2 = combi[j][1]
            sns = solve_snr(mods["sn"][i,b1],mods["sn"][i,b2],bests["sn"][b1,b2])
            if baseline[i,b1]*baseline[i,b2]*baseline[b1,b2]==0 or mods["sn"][i,b1]==1e8 or mods["sn"][i,b2]==1e8 or bests["sn"][b1,b2]==1e8:
                good_snrs[i,j]=-1.
                best_snrs[b1,best_snct[b1]]=-1.
                best_snrs[b2,best_snct[b2]]=-1.
            else:
                good_snrs[i,j]=sns[0]
                best_snrs[b1,best_snct[b1]]=sns[1]
                best_snrs[b2,best_snct[b2]]=sns[2]
            best_snct[b1]+=1
            best_snct[b2]+=1
                   

    # should use np.min? median looks more like the standard curve

    alls = good_snrs.copy()
    good_snrs = np.zeros((mods.size/nbest))
    good_snr_errs_ma = np.zeros(mods.size/nbest)
    good_snr_errs_mi = np.zeros(mods.size/nbest)
    for i in size_by_nbest_range:
        if len((alls[i])[alls[i]>0.0])==0:
            good_snrs[i] = 0.0
            good_snr_errs_ma[i] = 0.0
            good_snr_errs_mi[i] = 0.0
        if len((alls[i])[alls[i]>0.0])==1:
            good_snrs[i] = ((alls[i])[alls[i]>0.0])[0]
            good_snr_errs_ma[i] = good_snrs[i]
            good_snr_errs_mi[i] = good_snrs[i]
        if len((alls[i])[alls[i]>0.0])>1:
            good_snrs[i] = np.median((alls[i])[alls[i]>0.0])
            good_snr_errs_ma[i] = np.max((alls[i])[alls[i]>0.0])
            good_snr_errs_mi[i] = np.min((alls[i])[alls[i]>0.0])

    alls = best_snrs.copy()
    best_snrs = np.zeros(nbest)
    best_snr_errs_ma = np.zeros(nbest)
    best_snr_errs_mi = np.zeros(nbest)
    for i in range(nbest):
        if len((alls[i])[alls[i]>0.0])==0:
            best_snrs[i] = 0.0
            best_snr_errs_ma[i] = 0.0
            best_snr_errs_mi[i] = 0.0
        if len((alls[i])[alls[i]>0.0])==1:
            best_snrs[i] = ((alls[i])[alls[i]>0.0])[0]
            best_snr_errs_ma[i] = best_snrs[i]
            best_snr_errs_mi[i] = best_snrs[i]
        if len((alls[i])[alls[i]>0.0])>1:
            best_snrs[i] = np.median((alls[i])[alls[i]>0.0])
            best_snr_errs_ma[i] = np.max((alls[i])[alls[i]>0.0])
            best_snr_errs_mi[i] = np.min((alls[i])[alls[i]>0.0])

    #good_snrs = np.median(good_snrs,axis=1)
    #best_snrs = np.median(best_snrs,axis=1)

    # final delays, phases and weights
    master_dels = np.zeros(nant)
    master_phs = np.zeros(nant)
    master_weights = np.zeros(nant)
    master_weight_errs_1 = np.zeros(nant)
    master_weight_errs_2 = np.zeros(nant)

    for i in range(nbest):
        master_dels[best_idxs[i]]=del_corrs[i]*tsamp
        master_phs[best_idxs[i]]=ph_corrs[i]
        master_weights[best_idxs[i]] = best_snrs[i]
        master_weight_errs_1[best_idxs[i]] = best_snr_errs_ma[i]
        master_weight_errs_2[best_idxs[i]] = best_snr_errs_mi[i]

    for i in range(ngood):
        master_dels[good_idxs[i]]=good_dels[i]*tsamp
        master_phs[good_idxs[i]]=good_phs[i]
        master_weights[good_idxs[i]] = good_snrs[i]
        master_weight_errs_1[good_idxs[i]] = good_snr_errs_ma[i] 
        master_weight_errs_2[good_idxs[i]] = good_snr_errs_mi[i] 


    # fudge to match current signs
    master_dels[best_idxs[0]:nant] *= -1.
    master_phs[best_idxs[0]:nant] *= -1.

    real_sefds = np.sqrt(2.*(1.5e7/640.)*args.tobs)*args.flux*master_weights.copy()

    # this is all a hack to get around zero division errors that are crocking up the data
    mask = master_weights > 0.0001
    master_weights[mask] = 1./master_weights[mask]
    master_weight_errs_1[mask] = 1./master_weight_errs_1[mask]
    master_weight_errs_2[mask] = 1./master_weight_errs_2[mask]
    _tmp = master_weight_errs_1[mask]-master_weight_errs_2[mask]
    _tmp[_tmp<0.0001] = 0.0001
    master_weight_errs_1[mask] = np.abs(master_weights[mask]/_tmp)

    #original code
    #master_weights[master_weights != 0.0] = 1./master_weights[master_weights != 0.0]
    #master_weight_errs_1[master_weights != 0.0] = 1./master_weight_errs_1[master_weights != 0.0]
    #master_weight_errs_2[master_weights != 0.0] = 1./master_weight_errs_2[master_weights != 0.0]
    #master_weight_errs_1[master_weights != 0.0] = np.abs(master_weights[master_weights != 0.0]/(master_weight_errs_1[master_weights != 0.0]-master_weight_errs_2[master_weights != 0.0]))
      


    master_weights /= master_weights.max()
    master_weights[master_weights <= 0.0] = 0.0
    master_weights[master_weights <= 0.001] = 0.0
    
    master_weights = np.sqrt(master_weights)
    if args.binary=='yes':
        master_weights[master_weights > 0.0] = 1.0

    # zero delays and phases for bad modules
    master_dels[master_weights==0.0] = 0.0
    master_phs[master_weights==0.0] = 0.0

    # print out delays file
    print 'ref '+ants[best_idxs[0]]
    for i in range(nant):
        print ants[i],master_dels[i],master_phs[i],master_weights[i]

    # actually save delays file
    f=open('calib.delays','w')
    f.write('ref '+ants[best_idxs[0]]+'\n')
    for i in range(nant):
        f.write(ants[i]+' '+str(master_dels[i])+' '+str(master_phs[i])+' '+str(master_weights[i])+'\n')
    f.close()

    # write analytics to file
    f=open('calib.out','w')
    f.write('#ref '+ants[best_idxs[0]]+'\n')
    best_ant_str = '#best ants: '
    for i in range(args.ntop):
        best_ant_str += ' '+ants[best_idxs[i]]
    f.write(best_ant_str+'\n')
    for i in range(nant):
        f.write(ants[i]+' '+pfb[i]+' '+str(master_dels[i])+' '+str(master_phs[i])+' '+str(master_weights[i])+' '+str(master_weight_errs_1[i])+' '+str(real_sefds[i])+' '+str(sca_sn[i])+'\n')
    f.close()

    if not args.noplot:
        plt.ion()
        plt.figure(1)

# form mask of bad ants
        badant = (master_weights**1)<args.wcut
        for i in range(len(badant)):
            if badant[i] is False and np.abs(master_dels[i]/tsamp)>0.5: 
                badant[i] = True

# delay correction versus antenna position
        plt.subplot(231)
        plt.xlabel('Antenna position (m)')
        plt.ylabel('Delay correction (samples)')
        #plt.title('Delays by position')
        plt.plot(dists,master_dels/tsamp,'go')
        plt.plot(dists[badant],(master_dels/tsamp)[badant],'ro')

# delay correction versus antenna weight
        plt.subplot(232)
        plt.xlabel('Antenna weight')
        plt.ylabel('Delay correction (samples)')
        #plt.title('Delays by weight')
        #plt.plot((master_weights**2),master_dels/tsamp,'go')
        #plt.plot((master_weights**2)[badant],(master_dels/tsamp)[badant],'ro')
        plt.plot((master_weights**1),master_dels/tsamp,'go')
        plt.plot((master_weights**1)[badant],(master_dels/tsamp)[badant],'ro')
        plt.ylim([-2.0,2.0]) # samples

# ranked module performance curve
        plt.subplot(233)
        plt.xlabel('Module rank')
        plt.ylabel('Module SNR')
        sefdmask = (np.sort((master_weights**1)))<args.wcut 
        plt.xlim(1,len(real_sefds))
        #plt.title('SNRs')
        plt.plot(np.arange(len(real_sefds)),np.sort(master_weights**2),'go',lw=2)
        plt.plot(np.arange(len(real_sefds))[sefdmask],np.sort(master_weights**2)[sefdmask],'ro',lw=2)
        plt.xlim(1,len(real_sefds))

# phase correction versus antenna position along array
        plt.subplot(234)
        plt.xlabel('Antenna position (m)')
        plt.ylabel('Phase correction (radians)')
        plt.ylim([-np.pi,np.pi])
        #plt.title('Phases by position')
        plt.plot(dists, master_phs, 'go')
        plt.plot(dists[badant], (master_phs)[badant], 'ro')

# phase correction versus antenna weight
        plt.subplot(235)
        plt.xlabel('Antenna weight')
        plt.ylabel('Phase correction (radians)')
        plt.ylim([-np.pi,np.pi])
        #plt.title('Phases by weight')
        plt.plot((master_weights**1),master_phs,'go')
        plt.plot((master_weights**1)[badant],master_phs[badant],'ro')

# antenna SNR versus position along array
        plt.subplot(236)
        plt.xlabel('Antenna position (m)')
        plt.ylabel('Module SNR')
        #plt.title('SNRs by position')
        plt.plot(dists,(master_weights**2.),'go')
        plt.plot(dists[badant],(master_weights**2.)[badant],'ro')

        ngood = len(real_sefds.nonzero()[0])
        print "ngood = ", ngood
        print "Adopted flux = ",args.flux
        print "Adopted tobs = ",args.tobs
        SEFD_average = ngood / (sum(1.0/real_sefds[np.nonzero(real_sefds)]))
        Tsys = SEFD_average*0.01
        print "Derived average SEFD ", SEFD_average
        print "Derived average Tsys ",Tsys


def delay_signal2(ar,delay):
    nbins = ar.size
    ramp = np.arange(0,np.pi*2,np.pi*2/nbins)
    ramp[nbins/2:]-=np.pi*2
    shift = np.array(np.exp(-1j*np.outer(delay,ramp)))
    return ar*shift

def delay_signal(ar,delay):
    nbins = ar.size
    ramp = np.arange(0,np.pi*2,np.pi*2/nbins)
    ramp[nbins/2:]-=np.pi*2
    shift = np.array(np.exp(-1j*ramp*delay))
    return ar*shift
    
def main(args):
    count = 0
    combi = []
    if args.acfile:
        count += args.nant
    if args.ccfile:
        combi = [i for i in combinations(range(args.nant),2)]
        count += len(combi)
        
    dtype = [("a","int32"),("b","int32"),("delay","float32"),
             ("sn","float32"),("ar","complex64",args.batch_size)]
    grp = [("sn","float32"),("del","float32"),("ph","float32")]
    pairs = np.recarray(count,dtype=dtype)
    
    nant = 0
    if args.acfile:
        acdata = np.fromfile(args.acfile,dtype="float32")
        nant = acdata.size/args.batch_size
        x = acdata.copy()
        acdata = acdata.reshape(nant,args.batch_size).astype("complex64")
        for ii,row in enumerate(acdata):
            ar = np.fft.ifft(row)
            sn = stupid_snr(np.fft.fftshift(abs(ar)),args.batch_size)
            pairs[ii] = (ii,ii,0.0,sn,ar)
            
    if args.ccfile:
        ccdata = np.fromfile(args.ccfile,dtype="complex64")
        ncombi = ccdata.size/args.batch_size
        ccdata = ccdata.reshape(ncombi,args.batch_size)
        
        # rank modules
        stupid_snrs = np.zeros(args.nant)
        stupid_number = np.zeros(args.nant)
        for ii,row in enumerate(ccdata):
            if not np.any(row==0j): # HIRES: originally: if np.sum(row)!=0j: 
                ar = np.fft.ifft(row)
                stupid_snr_value = stupid_snr(np.fft.fftshift(abs(ar)),args.batch_size)
                stupid_snrs[combi[ii][0]] += stupid_snr_value
                stupid_snrs[combi[ii][1]] += stupid_snr_value
                stupid_number[combi[ii][0]] += 1
                stupid_number[combi[ii][1]] += 1

        sca_sn = stupid_snrs/(1.*stupid_number)

        sca_sn = np.nan_to_num(sca_sn) # HIRES: remove nan's to identify best antennas
        for i in range(args.nant):
            if sca_sn[i] >= args.thresh and stupid_number[i]<args.nused:
                sca_sn[i] = args.thresh
                #print i, sca_sn[i],stupid_number[i],stupid_snrs[i]

        # module statii - -1 is dead, 0 is average, 1 is best-args.ntop
        mod_status = np.zeros(args.nant)
        mod_status[sca_sn < args.thresh] = -1
        mod_status[sca_sn > (np.sort(sca_sn))[args.nant-args.ntop-1]] = 1
        if (np.sort(sca_sn))[args.nant-args.ntop-1] <= args.thresh:
            print 'Dont have enough modules with a good SNR'
            print args.thresh, (np.sort(sca_sn))[args.nant-args.ntop-1]
            sys.exit(0)
        good_idxs = np.arange(mod_status.size)[mod_status==0]
        best_idxs = np.arange(mod_status.size)[mod_status==1]
        mods = np.recarray((good_idxs.size,args.ntop),dtype=grp)
        bests = np.recarray((args.ntop,args.ntop),dtype=grp)

        for ii,row in enumerate(ccdata):
            
            mydel = 0.0
            myph = 0.0
            sn = 0.0

            # if (mod_status[combi[ii][0]]==1 or mod_status[combi[ii][1]]==1):
            if (not np.any(row==0j) and (mod_status[combi[ii][0]]==1 or mod_status[combi[ii][1]]==1)) : # HIRES: original line above

                if np.sum(row)==0j:
                    print 'bad sum'
                    sys.exit(1)
                myfit = new_snr(row,args.batch_size,args)
                sn = myfit[3]
                ar = np.fft.ifft(myfit[0])
                mydel = myfit[1]
                myph = myfit[2]
                #print sn,mydel,myph
                            
                if mod_status[combi[ii][0]]==0.:
                    mods[good_idxs == combi[ii][0],best_idxs == combi[ii][1]]=(sn,mydel,myph)
                elif mod_status[combi[ii][1]]==0.:
                    mods[good_idxs == combi[ii][1],best_idxs == combi[ii][0]]=(sn,mydel,myph)
                else:
                    mult1 = np.arange(args.ntop)[best_idxs == combi[ii][0]]
                    mult2 = np.arange(args.ntop)[best_idxs == combi[ii][1]]
                    bests[mult1,mult2]=(sn,mydel,myph)
                    bests[mult2,mult1]=(sn,mydel,myph)
                
            pairs[ii+nant] =(combi[ii][0],combi[ii][1],mydel,sn,ar)


    proc_ants(mods,bests,mod_status,good_idxs,best_idxs,args,sca_sn)
    if not args.noplot:
        x = CorrPlot(pairs, batch_size=args.batch_size, antfile=args.antfile,
                ntop=args.ntop, best_idxs=best_idxs, nukeF1 = args.nukeF1,
                nramps=args.nramps, wcut=args.wcut, good_idxs=good_idxs, mods=mods, args=args)
        plt.ion()
        plt.show()
        raw_input(">>> Press to do something unexpected")
    return 1

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--batch_size", help="number of channels per spectrum",
                        type=int, default=1024)
    parser.add_argument("-n","--nant", help="number of antennas",
                        type=int, default=None)
    parser.add_argument("-a","--acfile", help="name of auto-correlation file (optional)",
                        type=str, default=None)
    parser.add_argument("-c","--ccfile", help="name of cross-correalation file (optional)",
                        type=str, default=None)
    parser.add_argument("-flux","--flux",help="source flux density in Jy (def 20.0)",type=float,default=20.0)
    parser.add_argument("-wcut","--wcut",help="lower limit for acceptable weights",type=float,default=0.2)
    parser.add_argument("-tobs","--tobs",help="observation length in seconds (def 300.0)",type=float,default=300.0)
    parser.add_argument("-top","--ntop",help="number of top antennas (default 4)",type=int,default=4)
    parser.add_argument("-nramps","--nramps",help="number ramps across single channel (default 20)",type=int,default=20)
    parser.add_argument("-antf","--antfile",help="obs.antenna file (no default)",type=str,default=None, required=True)
    parser.add_argument("-binw","--binary",help="binary weights (def no)",type=str,default='no')
    parser.add_argument("-noplot","--noplot",help="plot results of calibration (delay, phase, snr) default yes"
           , default=False, action="store_true")
    parser.add_argument("-tsamp","--tsamp",help="sample time in microseconds (def 0.064 us)",type=float,default=0.064)
    parser.add_argument("-thresh","--thresh",help="threshold for stupid_snr rejection (def=10)",type=float,default=10.)
    parser.add_argument("-nchan","--nchan",help="number of coarse channels (def 20)",type=int,default=20)
    parser.add_argument("-nused","--nused",help="stupid_number threshold (def 88)",type=int,default=88)
    parser.add_argument("-nukeF1","--nukeF1",help="Apply brute force Telstra F1 removal",
             default=False, action="store_true")
    parser.add_argument("-verbose","--verbose",help="Print out the values at pick up",
             default=None)
    args = parser.parse_args()

    if args.nant == None:
        raise Exception("Number of antennas must be specified")
    if not os.path.isfile(args.antfile):
        raise Exception("obs.antenna file does not exist.")

    main(args)
