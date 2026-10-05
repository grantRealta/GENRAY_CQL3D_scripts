import numpy as np
import matplotlib.pyplot as plt
from numpy import heaviside
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.rc('axes', labelsize = 16)
plt.rc('axes', titlesize = 16)
plt.rc('figure', titlesize = 14)
plt.rc('legend', fontsize = 12)

###
# Based off work by Andrew Seltzman
# see https://iopscience.iop.org/article/10.1088/1741-4326/ab22c8/meta
#
# Takes in either a target N|| and outputs what phase is required to produce that value
# or takes in a set of module phases and outputs the spectrum
# The Fourier transform can either be done analytically (much faster but requries that all module powers are equal)
# or numerically (slower but more general)
#
# Assumes that the produced N|| is negative. If the field direction means it's positive, just change the sign of the results in whatever script uses the outputs
#
# Returns  
#   peakNparas - the peaks of the N|| lobes
#   peakEdges - ~the minima on either side of the N|| lobes (only exactly true for symmetrical lobes)
#   directivities - the directivity of each peak

###

eps_0 = 8.854187817e-12
c = 299792458 #m/s
f_rf = 4.6e9 #Hz, RF frequency
w_rf = 2*np.pi*f_rf #rad/s, angular frequency

def generateSpectrum(target_npara=None, #if not supplying the phase shift, what is your target N||
                    modulePhaseShift = None, #if not supplying the target N||, what is the phase shift between modules
                    analytic = True, #if you want to calculate things analytically. Only valid for equal powers in all WGs
                    modulePowerRatio = None, #power ratios between in the modules
                    doPlot = False, #False -> no plot, 'spectrum' -> plot the spectrum, 'both' -> plot the spectrum and waveguide power+phase

                    delta = 0.001/2, #[m], septum width #septa on left and right should only be half a mm
                    w_wg = 0.005, #[m], waveguide width (width here being the horizontal dimension, which is the smaller dimension for LHCD grills)
                    w_spacer = 0.00647, #[m], spacer width between each module
                    num_wg = 6, #number of output waveguides module
                    num_module = 8, #number of launcher modules)
                    delta_phi = np.pi/2, #[rad/s], phase shift between waveguides ##This is due to how the antenna was designed. Cannot be tuned
                    numLobes = 4, #number of lobes to return data on

                    #GENRAY does not take in an arbitrary spectrum. Instead it uses basis functions
                    #if true, returned values produce the best match to the predicted spectrum when fed into GENRAY
                    returnGenrayBestFit = True 
                    ):

    #unless told otherwise, assume equal power across the modules
    if modulePowerRatio is None:
        modulePowerRatio = np.ones(num_module)

    modulePowerRatio = np.array(modulePowerRatio)#in case it was given as a list

    #width of a module
    # each aperture of the module has a septum on either side
    #the outer walls of the module are thus half the thickness of the inner walls
    w_module = num_wg*(delta + w_wg + delta) 
    #width of the full grill
    w_grill = num_module*w_module + (num_module-1)*w_spacer 

    #if the target N|| is specified
    if target_npara is not None and modulePhaseShift == None:
        optActiveWGPhaseShift = -target_npara*w_rf*(w_wg+delta*2)/c #optimal phase shift between active elements in a module
        optModuleWGPhaseShift =  num_wg*optActiveWGPhaseShift  - (target_npara*w_rf*(w_spacer)/c)#optimal phase shift between modules

        modulePhaseShift =optModuleWGPhaseShift

        print(f'actual, optimal phasing between active elements: {np.degrees(delta_phi)%360}, {np.degrees(optActiveWGPhaseShift)%360} degrees')
        print(f'Optimal phasing between modules: {np.degrees(modulePhaseShift)%360} degrees')

    if analytic:
        assert isinstance(modulePhaseShift, (int, float))

    #rectangular function of height 1 and width 1 centered at the origin
    def rect(x):
        return heaviside(x+1/2,.5)-heaviside(x-1/2,.5)

    #x axis for the output of the Fourier transform
    N_paras = np.linspace(-20, 20, int(5e4+1))#needs to be an odd number of points

    if doPlot == 'both':
        fig,axes = plt.subplots(nrows=2,gridspec_kw={'height_ratios': [3, 1]}, figsize = (6.5, 7))
        axes[0].set_ylabel(r'$P(N_{||})$')
        axes[0].set_xlabel(r'$N_{||}$')
        axes[1].set_xlabel(f'z (m)')
        axes[1].set_ylabel(f'E (normalized)')

    if doPlot == 'spectrum':
        fig,ax = plt.subplots()
        ax.set_ylabel(r'P$(N_{||})$ (Normalized)')
        ax.set_xlabel(r'$N_{||}$')

    assert len(modulePowerRatio) == num_module

    # Create the array of z points
    z_points = 2001
    z = np.arange(-w_wg, w_grill + w_wg, (w_grill + 2*w_wg) / z_points)

    # Differential element (width of array for integration)
    dz = z[1] - z[0]

    # Array of electric field and phase
    E_z = np.zeros_like(z)
    phi_z = np.zeros_like(z)

    z_wg = (delta + w_wg + delta)
    z_module = (num_wg * z_wg + w_spacer)
    # Loop through modules, groups, and waveguides

    #make the phase and electric field arrays across the grill face
    for module_index in range(num_module):
        delta_phi_module = 0
        if isinstance(modulePhaseShift,(int,float)):
            delta_phi_module =  module_index * modulePhaseShift
        else:
            delta_phi_module = modulePhaseShift[module_index]
        for wg_index in range(num_wg):

            wg_z_window = rect((z - (module_index * z_module + wg_index * z_wg + w_wg / 2)) / w_wg)

            E_j = modulePowerRatio[module_index] * wg_z_window  # Electric field at wg
            E_z += E_j

            phi_j = (wg_index * delta_phi +delta_phi_module) * wg_z_window + np.pi
            phi_z += phi_j

    if doPlot == 'both':
        axes[1].plot(z, E_z)
        phi_ax = axes[1].twinx()
        phi_ax.plot(z,(phi_z*180/np.pi)%360, color = 'r', linestyle = 'dashed')
        phi_ax.set_ylabel('phase (deg)')
        phi_ax.set_ylim([0,360])
        phi_ax.set_yticks([0,90,180,270,360])

    # require requested module powers are all equal to do analytic form
    if analytic == True and np.equal(modulePowerRatio,[1]*num_module).all():

        alpha_active = delta_phi + N_paras*w_rf*(delta + w_wg + delta)/c

        with np.errstate(divide='ignore'):
            P_active = (np.sin(N_paras*w_rf*w_wg/(2*c))**2 *
                    np.sin(num_wg*alpha_active/2)**2 /
                    (N_paras**2*np.sin(alpha_active/2)**2))

        alpha_module = modulePhaseShift+N_paras*w_rf*(1*(w_spacer + num_wg*(delta + w_wg + delta)))/c
        P_module = (np.sin(num_module*alpha_module/2)**2 /
                    (np.sin(alpha_module/2)**2))

        P_total = P_active*P_module

    else:
        FT_E = np.zeros_like(N_paras, dtype=complex)
        E_complex = E_z * np.exp(-1j * phi_z)

        # Compute the Fourier transform using numerical integration
        for n_index in range(len(N_paras)):
            FT_E[n_index] = np.sum(E_complex * np.exp(-1j * N_paras[n_index] * w_rf / c * z) * dz)

        # Compute P(n) = E(n) * conj(E(n))
        P_total = np.abs(FT_E * np.conj(FT_E))

        
    P_total = P_total/np.nanmax(P_total)


    for nanIndex in np.argwhere(np.isnan(P_total)):
        P_total[nanIndex] = (P_total[nanIndex-1] + P_total[nanIndex+1])/2

    totalArea = np.trapz(P_total, x=N_paras)

    peaksIndices, peakProps = find_peaks(P_total, height = 0.02)

    peakIndicesOfInterest = peaksIndices[np.flip(peakProps["peak_heights"].argsort())][:numLobes]
    minimaPeaks, _ = find_peaks(-P_total)

    directivities = np.zeros(numLobes)
    peakNparas = np.zeros(numLobes)
    peakEdges = np.zeros((numLobes,2))
    for i in range(numLobes):
        if i >= len(peakIndicesOfInterest):
            continue
        peak = peakIndicesOfInterest[i]

        min1 = minimaPeaks[minimaPeaks < peak][-1] #min to the left of our peak
        min2 = minimaPeaks[minimaPeaks > peak][0] #min to the right of our peak
        
        directivities[i] = np.trapz(P_total[min1:min2], x = N_paras[min1:min2])/totalArea
        peakNparas[i] = N_paras[peak]
        width = N_paras[min2] - N_paras[min1]

        peakEdges[i,0] = peakNparas[i]-width/2#N_paras[min1]
        peakEdges[i,1] = peakNparas[i]+width/2#N_paras[min2]

        width = peakEdges[i,1]-peakEdges[i,0]
        
        if doPlot == 'both':
            axes[0].fill_between(N_paras[min1:min2], P_total[min1:min2], np.zeros(len(N_paras[min1:min2])), zorder = 10)  
            axes[0].plot([1,1],[10,10],label = r'$N_{||}$'+ f'= {N_paras[peak]:.2f}, {width:.2f} width \n Directivity = {directivities[i]:.2f}', lw = 5)
        if doPlot == 'spectrum':
            ax.fill_between(N_paras[min1:min2], P_total[min1:min2], np.zeros(len(N_paras[min1:min2])), zorder = 10)  
            #ax.plot([1,1],[10,10],label = r'$N_{||}$'+ f'= {N_paras[peak]:.2f}, {width:.2f} width,\n Directivity = {directivities[i]:.2f}', lw = 5)
            ax.plot([1,1],[10,10],label = r'$N_{||}$'+ f'= {N_paras[peak]:.2f}, \nDirectivity = {directivities[i]:.2f}', lw = 5)

        print(f'peak at {peakNparas[i]:.3f} with {directivities[i]} directivity and a width of {width}')

    if returnGenrayBestFit:

        def lobes(x, *params):
            num_lobes = len(params)//3
            y = np.zeros(len(x))
            for i in range(num_lobes):
                peak = params[i*3]
                amp = params[i*3+1]
                width = params[i*3+2]

                f = 2*np.pi*(x - peak)/width

                y += amp*np.sin(f)**2/f**2

            return y

        numDownSamplePoints = 1000
        downsampledNpara = np.linspace(np.min(N_paras), np.max(N_paras), numDownSamplePoints)
        downsampledP = interp1d(N_paras,P_total)(downsampledNpara)

        p0 = []
        bounds = []
        for i in range(numLobes):
            width = peakEdges[i,1]-peakEdges[i,0]
            p0.append(peakNparas[i])
            p0.append(1)
            p0.append(width)

            bounds.append([peakNparas[i]-width/2, peakNparas[i]+width/2])
            bounds.append([.05,1])
            bounds.append([width/2, width*2])
            
        popt,pcov = curve_fit(lobes, downsampledNpara, downsampledP, p0=p0, bounds = np.array(bounds).T) 

        genrayFit = lobes(N_paras, *popt)
        genrayArea = np.trapz(genrayFit, x=N_paras)

        genrayPeaksIndices, genrayPeakProps = find_peaks(genrayFit, height = 0.02)
        genrayPeakIndicesOfInterest = genrayPeaksIndices[np.flip(genrayPeakProps["peak_heights"].argsort())][:numLobes]
        genrayMinimaPeaks, _ = find_peaks(-genrayFit)

        #directivities = np.zeros(numLobes)
        peakNparas = np.zeros(numLobes)
        peakEdges = np.zeros((numLobes,2))
        for i in range(numLobes):
            if i >= len(genrayPeakIndicesOfInterest):
                continue
            peak = genrayPeakIndicesOfInterest[i]

            min1 = genrayMinimaPeaks[genrayMinimaPeaks < peak][-1] #min to the left of our peak
            min2 = genrayMinimaPeaks[genrayMinimaPeaks > peak][0] #min to the right of our peak
            
            #directivities[i] = np.trapz(genrayFit[min1:min2], x = N_paras[min1:min2])/genrayArea
            peakNparas[i] = N_paras[peak]
            width = N_paras[min2] - N_paras[min1]
            peakEdges[i,0] = peakNparas[i]-width/2#N_paras[min1]
            peakEdges[i,1] = peakNparas[i]+width/2
            print(f'genray peak at {peakNparas[i]:.3f} with width of {width}')


    if doPlot == 'both':
        if returnGenrayBestFit:
            axes[0].plot(N_paras, genrayFit, color = 'lime',lw = 4, zorder = 1, label = 'GENRAY best fit')
        axes[0].plot(N_paras, P_total, color = 'k',zorder = 1)
        axes[0].set_xlim([-10,10])
        axes[0].set_ylim([-.05,1.05])
        axes[0].legend()
        fig.tight_layout()
        plt.show()
    if doPlot == 'spectrum':
        #ax.plot(N_paras, genrayFit, color = 'r',lw = 1.5, zorder = 10, label = 'GENRAY spectrum', linestyle = 'dashed')
        ax.plot(N_paras, P_total, color = 'k',zorder = 10, lw = 1.25)#, label = 'Experimental spectrum')
        ax.set_xlim([-10,10])
        ax.set_ylim([-.05,1.05])
        ax.legend()
        fig.tight_layout()
        #plt.savefig('DIIID_203912.02700_exSpectrum.jpeg',dpi=300)

        plt.show()

    return peakNparas, peakEdges, directivities

def phaseMatrix():
    phases = np.linspace(0,2*np.pi,100)
    nparas = np.linspace(-50, 50, int(5e4+1))#needs to be an odd number of points

    powerMatrix = np.zeros((len(phases), len(nparas)))
    for i in range(len(phases)):
        phase = phases[i]
        peaks, edges, directs, power = generateSpectrum(modulePhaseShift = phase, analytic = False, powerRatio = [1,1,1,0,0,0,0,0])
        powerMatrix[i,:] = power

    fig,ax = plt.subplots()
    p2 = ax.pcolormesh(phases,nparas, powerMatrix.T,shading = 'nearest',cmap='viridis', vmin=np.nanmin(0), vmax = np.nanmax(powerMatrix))
    cbar = fig.colorbar(p2, ax = ax, shrink = .8, pad = .01)
    cbar.set_label(r"P(N$_{||})/P_{max}(N_{||})$")
    ax.set_ylim([-5,1])
    ax.set_ylabel(r'N$_{||}$')
    ax.set_xlabel(r'Module phase difference (rad)')
    ax.set_title('3 modules')
    fig.tight_layout()
    plt.show()

def phasePlot():
    phases = np.linspace(0, 2*np.pi,100)
    nparas = np.zeros(len(phases))
    fig,ax = plt.subplots()

    for i in range(len(phases)):
        phase = phases[i]
        peaks, _, _,_ = generateSpectrum(modulePhaseShift = phase)
        nparas[i] = peaks[0]

    #phases = 2*np.pi - phases
    phases = phases%(2*np.pi)
    phases[phases < 1.52] += 2*np.pi
    toFit_phases = np.copy(phases)
    toFit_npara = np.copy(nparas)

    m,b = np.polyfit(toFit_phases, toFit_npara, 1)
    print(f'slope: {m}, intercept: {b}')
    ax.scatter(toFit_phases,toFit_npara)
    ax.plot(phases, m*phases+b)
    plt.show()

#plots the directivity of the DIII-D launcher over a range of N||
def plotDirectivities():
    targets = np.arange(-3.5,-2,.1)

    fig,ax = plt.subplots()#figsize= (5.5,4.8))

    for target in targets:
        peakNparas, peakEdges, directivities= generateSpectrum(target_npara=target, analytic = True, 
                                            doPlot = False)
        ax.scatter([np.abs(peakNparas[0])], [directivities[0]], color ='k', marker = 'D')

    ax.set_ylabel('Directivity')
    ax.set_xlabel(r'Target |N$_{||}$|')
    ax.set_ylim([0,1])
    ax.set_xticks([3.5,3.3,3.1,2.9,2.7,2.5,2.3,2.1,1.9])
    fig.tight_layout()
    #plt.savefig('DIIID_8modDirect.jpeg',dpi=300)

    plt.show()

if __name__ == '__main__':

    generateSpectrum(-2.7, analytic = True, doPlot = 'spectrum', numLobes = 3)
    generateSpectrum(-2.7, analytic = True, doPlot = 'both', numLobes = 3, returnGenrayBestFit= True)
    #plotDirectivities()
