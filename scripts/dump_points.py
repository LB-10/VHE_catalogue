import vhe_spectra
import numpy as np
agn_meas =  vhe_spectra.AGNcollection()
agn_meas.get_agn_measurement_collection(PATH="../AGN/")
agn_meas.convert_collection_TeV_vs_TeVcm2s()
agn_m = agn_meas.m
agn_m = sorted(agn_m, key=lambda t: t['year'])
elist=np.array([])
flux=np.array([])
dflux=np.array([])
zlist=np.array([])
name=np.array([])

for m in agn_m:
# remove 
    if  m['z_uncertain']:
	continue
    if m['id']=='Mkn421_HEGRA_1999':   # not in the BV sample
	continue
    if m['id']=='Mkn421_HESS_2010':
	continue
    if m['id']=='Mkn421_MAGIC_2010_flare':
	continue
    if m['id']=='Mkn421_MAGIC_2011':
	continue
    if m['id']=='Mkn421_WHIPPLE_2002':
	continue
    if m['id']=='Mkn501_CAT_1999':
	continue
    if m['id']=='Mkn501_MAGIC_2011mw':
	continue
    if m['id']=='Mkn501_VERITAS_2009_Huang_a':
	continue
    if m['id']=='Mkn501_VERITAS_2009_Huang_b':
	continue
    if m['id']=='Mkn501_VERITAS_2009':
	continue
    if m['id']=='Mkn501_VERITAS_2011mw':
	continue
    if m['id']=='Mkn501_VERITAS_2011mwflare':
	continue
    if m['id']=='1ES1959_HEGRA_2003ls':
	continue
    if m['id']=='1ES1959_MAGIC_2006':
	continue
    if m['id']=='PKS2005_HESS_2005':
	continue
    if m['id']=='PKS2005_HESS_2011':
	continue
    if m['id']=='CenA_HESS_2009':
	continue
    if m['id']=='M87_MAGIC_2008':
	continue
    if m['id']=='M87_VERITAS_2007':
	continue
    if m['id']=='M87_VERITAS_2010':
	continue
    if m['id']=='M87_HEGRA_1999':
	continue
    if m['id']=='M87_HESS_2004':
	continue
    if m['id']=='M87_HESS_2005':
	continue
    if m['id']=='PKS0548_HESS_2010':
	continue
    if m['id']=='WCom_VERITAS_2009':
	continue
    if m['id']=='PKS2155_HESS_2005_stereo':
	continue
    if m['id']=='PKS2155_HESS_2007':
	continue
    if m['id']=='PKS2155_HESS_2007':
	continue
    if m['id']=='H1426_COMBINED_2003':
	continue
    if m['id']=='H2356_HESS_2006a':
	continue
    if m['id']=='H2356_HESS_2006b':
	continue
    if m['id']=='1ES1218_VERITAS_2013':
	continue
    if m['id']=='1ES1011_MAGIC_2012':
	continue
    if m['id']=='3C279_MAGIC_2011':
	continue
    if m['id']=='PKS1424_MAGIC_2011':
	continue
    if m['id']=='PKS1424_VERITAS_2009':
	continue
    if m['id']=='1ES0414_VERITAS_2013':
	continue
    if m['id']=='PKS1510_HESS_2013':
	continue
    if m['id']=='PKS1222_MAGIC_2011':
	continue
    if m['id']=='3C279_MAGIC_2008':
	continue
    if m['id']=='Mkn501_HEGRA_1999':
	continue
    

	
    
    elist = np.append(elist,np.array(m['data']['energy']))
    flux  = np.append(flux,np.array(m['data']['flux']))
    dflux = np.append(dflux,np.array(m['data']['f_err_sym']))
    for current in xrange(0,len(np.array(m['data']['energy']))):
	 name=np.append(name,m['id'])
	 zlist=np.append(zlist,m['z'])
    
np.savetxt('list.all',np.transpose([elist,flux,dflux,name,zlist]),
	fmt=['%s','%s','%s','%s','%s'])	
