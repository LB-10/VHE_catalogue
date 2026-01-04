"""
Collection of VHE measurements of AGN

History of changes:
Version 1.0
- Created 9th February 2011

"""

__version__ = 1.0
__author__ = "M. Meyer // manuel.meyer@physik.uni-hamburg.de"

#===============================================================================================
import numpy as np
import glob
import warnings

class AGNcollection:
    def __init__(self):
	self.m = []

    def convert_collection_TeV_vs_TeVcm2s(self):
	new_m = []
	for m in self.m:
	    # Convert Energy to TeV 
	    if m['unit_x'] == 'GeV':
		for i,x in enumerate(m['data']['energy']):
		    m['data']['energy'][i] *= 1e-3
		m['data']['x_err'] *= 1e-3
		m['data']['x_err_sym'] *= 1e-3
	    # Convert Flux and Flux error to 1/(TeV cm^2 s)
	    if m['unit_f'] == '1/(TeVm^2s)' or m['unit_f'] == 'TeV/(m^2s)':
		m['data']['flux'] = map(lambda x: x*1e-4,m['data']['flux'])
		if len(m['data']['f_err_sym']):
		    m['data']['f_err_sym'] = map(lambda x: x*1e-4,m['data']['f_err_sym'])
		else:
		    m['data']['f_err_low'] = map(lambda x: x*1e-4,m['data']['f_err_low'])
		    m['data']['f_err_high'] = map(lambda x: x*1e-4,m['data']['f_err_high'])
	    if m['unit_f'] == 'TeV/(cm^2s)' or m['unit_f'] == 'erg/(cm^2s)' or m['unit_f'] == 'GeV/(cm^2s)':
		if m['unit_f'] == 'erg/(cm^2s)':
		    corr = 1.602
		elif m['unit_f'] == 'GeV/(cm^2s)':
		    corr = 1000.
		else:
		    corr = 1.
		for i,x in enumerate(m['data']['flux']):
			m['data']['flux'][i] /= (corr * m['data']['energy'][i]*m['data']['energy'][i])
			if len(m['data']['f_err_sym']):
			    m['data']['f_err_sym'][i] /= (corr*m['data']['energy'][i]*m['data']['energy'][i])
			else:
			    m['data']['f_err_low'][i] /= (corr*m['data']['energy'][i]*m['data']['energy'][i])
			    m['data']['f_err_high'][i] /= (corr*m['data']['energy'][i]*m['data']['energy'][i])
	    new_m.append(m)
	self.m = new_m
	return

    def get_agn_measurement_collection(self,PATH="./",All=False) :
	sources = glob.glob(PATH + "*/")
	for i,source in enumerate(sources):
	    #spectra = glob.glob("./" + source + "*.dat")
	    if All:
		spectra = glob.glob(source + "*.info")
	    else:
		spectra = glob.glob(source + "*.dat")
	    for spectrum in spectra:
		if spectrum == source + 'FERMI_2FGL.dat':
		    continue
		if All:
		    spectrum = spectrum.split(".info")[0] + ".dat"
		    try:
			data = np.loadtxt(spectrum)
		    except:
			data = np.zeros((1,1))
		else:
		    try:
			data = np.loadtxt(spectrum)
		    except ValueError:
			warnings.warn("Cannot handle data in File %s\n \
			Do all rows have the same lenght?" % spectrum, Warning)
			exit(-1)
		m = {}
		m['id'] = source.split("/")[-2] + "_" + spectrum.split("/")[-1].split(".dat")[0]
		m['data'] = {}
		if len(data[0,:]) == 4:
		# Check, how many energies there are in each row
		# if there are two then we have BinMin and BinMax
		# otherwise its FluxErrMin and FluxErrMax
		    EnergyIndex	= np.where(data[-1,:] > 1E-2)[0]
		    FluxIndex	= np.where(data[-1,:] < 1E-2)[0]
		    if len(EnergyIndex) == 2:
			m['data']['energy'] = (data[:,EnergyIndex[1]] + data[:,EnergyIndex[0]]) / 2.
			m['data']['flux'] = data[:,FluxIndex[0]]
			m['data']['f_err_sym'] = data[:,FluxIndex[1]]
		    elif len(FluxIndex) == 3:
			m['data']['energy'] = data[:,EnergyIndex[0]]
			m['data']['flux'] = data[:,FluxIndex[0]]

			m['data']['f_err_sym'] = np.zeros(len(m['data']['energy']))
			m['data']['f_err_sym'][np.where(data[:,FluxIndex[1]] == data[:,FluxIndex[2]])] = \
			    data[np.where(data[:,FluxIndex[1]] == data[:,FluxIndex[2]]),FluxIndex[1]]
			m['data']['f_err_sym'][np.where(data[:,FluxIndex[1]] > data[:,FluxIndex[2]])] = \
			    data[np.where(data[:,FluxIndex[1]] >= data[:,FluxIndex[2]]),FluxIndex[1]] \
				- data[np.where(data[:,FluxIndex[1]] >= data[:,FluxIndex[2]]),FluxIndex[0]]
			m['data']['f_err_sym'][np.where(data[:,FluxIndex[1]] < data[:,FluxIndex[2]])] = \
			    data[np.where(data[:,FluxIndex[1]] < data[:,FluxIndex[2]]),FluxIndex[0]] \
				- data[np.where(data[:,FluxIndex[1]] < data[:,FluxIndex[2]]),FluxIndex[1]]
			if len(np.where(m['data']['flux']-m['data']['f_err_sym'] < 0.)[0]):
			    warnings.warn("Error larger than flux in %s" % m['id'], Warning)
		    else:
			warnings.warn("Cannot handle data in Table in File %s" % spectrum, Warning)
		elif len(data[0,:]) == 3:
		    m['data']['energy'] = data[:,0]
		    m['data']['flux'] = data[:,1]
		    m['data']['f_err_sym'] = data[:,2]
		# The next condition is the case for Mkn501_HEGRA_1999 spectrum
		elif len(data[0,:]) == 5 and len(EnergyIndex) == 1:
		    m['data']['energy'] = data[:,0]
		    m['data']['flux'] = data[:,1]
		    m['data']['f_err_sym'] = data[:,2]
		    m['data']['f_err_sys'] = (data[:,3] + data[:,4]) / 2.
		elif len(data[0,:]) == 1 and data[0,0] == 0.:
		    ULIndex = []
		    m['data']['energy'] = []
		    m['data']['flux'] = []
		    m['data']['f_err_sym'] = []
		    warnings.warn("No Data in File %s" % spectrum, Warning)
		    continue
		else:
		    warnings.warn("Unknown number of columns in file %s" % spectrum, Warning)
		# Find Upper Limits

		if not len(data[0,:]) == 1:
		    ULIndex = np.where(m['data']['f_err_sym'] == 0.)[0]
		if len(ULIndex):
		    m['data']['UL_energy'] = m['data']['energy'][ULIndex[0]:]
		    m['data']['energy'] = m['data']['energy'][:ULIndex[0]]
		    m['data']['UL_flux'] = m['data']['flux'][ULIndex[0]:]
		    m['data']['flux'] = m['data']['flux'][:ULIndex[0]]
		    m['data']['f_err_sym'] = m['data']['f_err_sym'][:ULIndex[0]]
		# calculate bin widths

		m['data']['x_err'] = np.zeros((2,len(m['data']['energy'])))	# energy of lower and upper bin end
		for iE, E in enumerate(m['data']['energy'][:-1]):
		    m['data']['x_err'][0,iE] = E * (m['data']['energy'][iE + 1] / E) ** -0.5
		    m['data']['x_err'][1,iE] = E * (m['data']['energy'][iE + 1] / E) ** 0.5
		m['data']['x_err'][0,-1] = m['data']['x_err'][1,-2]
		m['data']['x_err'][1,-1] = m['data']['energy'][-1]**2. / m['data']['x_err'][0,-1]
		m['data']['x_err_sym'] = (m['data']['x_err'][1] - m['data']['x_err'][0])/2.

		# Read auxiliary Information
		try:
		    infofile = spectrum.split(".dat")[0] + ".info"
		    f = open(infofile, "r")
		    for l in f.readlines():
			try:
			# for the year
			    m[l.split(":")[0]] = int(l.split(":")[1])
			except:
			    try:
				# for the redshift 
				m[l.split(":")[0]] = float(l.split(":")[1])
			    except:
				# for anything else
				if len(l.split(":")) == 1:
				    m[l.split(":")[0].strip("\t\n ")] = ''
				else:
				    m[l.split(":")[0].strip("\t\n ")] = l.split(":")[1].strip("\t\n")
		    f.close()
		except:
		    warning = "No auxiliary information for %s" % infofile
		    warnings.warn(warning, Warning)
		    pass
		self.m.append(m)
		#for i in agn_measurements:
    		#print i 

		#for i,s in enumerate(sources):
    		#print i,s
	return
#===============================================================================================
class FermiSpectra:
    def __init__(self):
	self.m = []

    def GetFermiSpectra(self,PATH="."):
	sources = glob.glob(PATH + "*/FERMI_2FGL.dat")
	for source in sources:
	    m = {}
	# Energy bins in GeV
	    m['EBinMin'] = np.array([0.1,0.3,1.,3.,10.])
	    m['EBinMax'] = np.array([0.3,1.,3.,10.,100.])
	    m['ELogBinCenter'] = (np.log10(m['EBinMin'] * m['EBinMax'])/2.)
	    m['EBinCenter'] = 10.**m['ELogBinCenter']

	    data = np.loadtxt(source)
	    m['flux'] = data[:,0]
	    m['err'] = data[:,1]
	    f = open(source,'r')
	    m['object'] = f.readlines()[0].strip('#\n')
	    f.close()
	    self.m.append(m)
#===============================================================================================
