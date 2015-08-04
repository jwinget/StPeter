#!/usr/bin/env python
#####################
# MS2 intensity-based label free quantification
# Algorithm by Griffin et al, Nat. Biotech 2010
# TPP integration by J. Winget
#####################

import argparse
import re
import sys
import pymzml
import bisect
import time
from lxml import etree
from numpy import interp
from math import log

start_time = time.time()

#-- Global dicts to hold some frequently used values --#
# Aminoacid_masses
md = { 'A': 71.03711,
   	'R': 156.1011,
   	'N': 114.04293,
   	'D': 115.02694,
	'C': 103.00919, # No carbam
#   	'C': 160.0306, # Incl. carbamidomethyl
   	'E': 129.04259,
   	'Q': 128.05858,
   	'G': 57.02146,
   	'H': 137.05891,
   	'I': 113.08406,
   	'L': 113.08406,
   	'K': 128.09496,
   	'M': 131.04049,
   	'F': 147.06841,
   	'P': 97.05276,
   	'S': 87.03203,
   	'T': 101.04768,
   	'W': 186.07931,
   	'Y': 163.06333,
   	'V': 99.06841,
   	'proton': 1.00728,
   	'H2O': 18.0155,
	'acetylation': 42.010565,
	'deamidation': 0.984016
   	}

# Namespaces, required by parser
NS = { 'ProtXML': 'http://regis-web.systemsbiology.net/protXML',
	'PepXML': 'http://regis-web.systemsbiology.net/pepXML',
	}
#-- End global dicts --#

def fast_iter(context, func, *args, **kwargs):
    """
    http://lxml.de/parsing.html#modifying-the-tree
    Based on Liza Daly's fast_iter
    http://www.ibm.com/developerworks/xml/library/x-hiperfparse/
    See also http://effbot.org/zone/element-iterparse.htm
    """
    for event, elem in context:
        func(elem, *args, **kwargs)
        # It's safe to call clear() here because no descendants will be
        # accessed
        elem.clear()
        # Also eliminate now-empty references from the root node to elem
        for ancestor in elem.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]
    del context

def calcFDR(protxml, fdr):
	''' Calculate the protein probability cutoff for a given FDR '''
	print('Processing '+protxml)
	probs = []
	fper = []
	def parseFDR(elem):
		probs.append(elem.get('min_probability'))
		fper.append(elem.get('false_positive_error_rate'))
	mytag = '{'+NS['ProtXML']+'}protein_summary_data_filter'
	context = etree.iterparse(protxml, tag=mytag)
	fast_iter(context, parseFDR)

	# Numpy interp requires these to be in reverse order
	probs.reverse()
	fper.reverse()
	pc = interp(fdr, fper, probs)
	print('Probability of '+str(pc)+' gives FDR of '+str(fdr))
	return pc

def getProts(protxml, pc):
	''' Parse proteins passing probability cutoff filter from protXML '''
	print('Parsing proteins')
	d = {} # This is our main data struct

	def parseProts(elem):
		prob = float(elem.get('probability'))
		if prob >= pc: # Keep it
			name = elem.get('protein_name')
			length = elem.xpath('ProtXML:parameter[@name="prot_length"]/@value', namespaces=NS)[0]
			try:
				desc = elem.xpath('ProtXML:annotation/@protein_description', namespaces=NS)[0]
			except:
				desc = ''
			d[name] = {}
			d[name]['length'] = int(length)
			d[name]['description'] = desc

			d[name]['peptides'] = {}
			peps = elem.xpath('ProtXML:peptide', namespaces=NS)
			for pep in peps:
				if pep.get('is_nondegenerate_evidence') == 'Y': # Keep it
					try: # Modified version if it exists
						seq = pep.xpath('ProtXML:indistinguishable_peptide/ProtXML:modification_info/@modified_peptide', namespaces=NS)[0]
					except: # Unmodified version
						seq = pep.get('peptide_sequence')
					seq = re.sub(r'C\[[0-9]{3}\]','C',seq) # Cys causes problems later
					d[name]['peptides'][seq] = {}

			if len(d[name]['peptides'].keys()) == 0: # Protein has no nondegenerate evidence
				del d[name]
				print(name+' had no unique peptides, discarding')

	mytag = '{'+NS['ProtXML']+'}protein'
	context = etree.iterparse(protxml, tag=mytag)
	fast_iter(context, parseProts)
	print('Parsed '+str(len(d.keys()))+' proteins')
	return d

def getDatabase(protxml):
	''' Parse the original search database from ProtXML '''
	dbs = []
	def parseDB(elem):
		thisdb = elem.get('reference_database')
		dbs.append(thisdb)
	mytag = '{'+NS['ProtXML']+'}protein_summary_header'
	context = etree.iterparse(protxml, tag=mytag)
	fast_iter(context, parseDB)
	db = dbs[0]
	return db

def calcLengths(db, nolengths):
	''' Parse a fasta file to generate protein lengths for proteins that need them'''
	lengths = {}
	return lengths

def checkLengths(d, protxml):
	''' Sometimes lengths 0 get written to ProtXML. Fun. '''
	nolength = []
	for protein in d.keys():
		if not d[protein]['length']:
			nolength.append(protein)
			del d[protein]
	if nolength: # We got some fixin' to do
		print("Could not find lengths for:")
		print(','.join(nolength))
		# Implement fix using database here
		#db = getDatabase(protxml)
		#print(db)
	return d

def getPepXML(protxml):
	''' Parse the source pepXML out of a protXML file '''

	sf = [] # This works but strings don't ::shrug::
	def parsePepXML(elem):
		sf.append(elem.get('source_files'))
	
	mytag = '{'+NS['ProtXML']+'}protein_summary_header'
	context = etree.iterparse(protxml, tag=mytag)
	fast_iter(context, parsePepXML)
	pepxml = sf[0]
	print('Found source file '+pepxml)
	
	return pepxml

def getSpectra(d, pepxml):
	''' Get scan events for the peptides of interest '''
	# Generate flat list of peptide sequences
	#all_peps = []
	#for protein, info in d.iteritems():
	#	all_peps += info['peptides'].keys()
	print('Parsing PSMs from pepXML')

	psms = {}
	counter = [0]
	def parsePSMs(elem):
		spectrum = elem.get('spectrum')
		mod_info = ''
		try: # Grab modified peptide sequence if present
			mod_info = elem.xpath('PepXML:search_result/PepXML:search_hit/PepXML:modification_info', namespaces=NS)
			seq = mod_info[0].get('modified_peptide')
		except: # Use unmodified sequence
			seq = elem.xpath('PepXML:search_result/PepXML:search_hit/@peptide', namespaces=NS)[0]

		if 'C[' in seq:
			# Some search engines write C[160] to PepXML, some don't. Isn't that fun.
			seq = re.sub(r'C\[[0-9]{3}\]','C',seq) # Cys is my least favorite amino acid

		# For speed, we first keep all PSMs and parse them into the struct later on
		counter[0] += 1
		try:
			psms[seq]['spectra'].append(spectrum)
		except: # Maybe spectral list doesn't exist
			try:
				psms[seq]['spectra'] = [spectrum]
			except: # Maybe this is a new sequence
				try:
					psms[seq] = {}
					psms[seq]['spectra'] = [spectrum]
				except:
					print('Error processing '+spectrum)
		if mod_info:
			mods = elem.xpath('PepXML:search_result/PepXML:search_hit/PepXML:modification_info/PepXML:mod_aminoacid_mass', namespaces=NS)
			modlist = []
			for m in mods:
				position = int(m.get('position'))
				mass = float(m.get('mass'))
				mt = (position, mass)
				modlist.append(mt)
			psms[seq]['modifications'] = modlist
		else: # No mods, just insert an empty list
			psms[seq]['modifications'] = []
		sys.stdout.write('\rParsed '+str(counter[0])+' total PSMs')


	mytag = '{'+NS['PepXML']+'}spectrum_query'
	context = etree.iterparse(pepxml, tag=mytag)
	fast_iter(context, parsePSMs)

	# Merge the matched PSMs with the main data struct
	print('\nMerging significant PSMs with protein data')
	for protein in d.keys():
		for peptide in d[protein]['peptides'].keys():
			try:
				spectra = psms[peptide]['spectra']
				mods = psms[peptide]['modifications']
				try:
					d[protein]['peptides'][peptide]['spectra'] += spectra
				except:
					d[protein]['peptides'][peptide]['spectra'] = spectra
				d[protein]['peptides'][peptide]['modifications'] = mods
			except:
				print('Did not match any PSMs to '+peptide)

	return d

def getRawData(pepxml):
	''' Parse raw data files from pepxml '''
	print('\nFinding raw data')

	raw_files = {}
	found_files = []
	def parseRawData(elem):
		base_name = elem.get('base_name')
		data_type = elem.get('raw_data')
		if data_type != '.mzML':
			if data_type == '':
				data_type = 'None'
			print('ERROR: Only mzML is currently supported, got '+data_type+' for '+base_name)
		else:
			raw_file = base_name + data_type
			raw_short = base_name.split('/')[-1]
			if raw_file not in found_files:
				raw_files[raw_short] = raw_file
				found_files.append(raw_file)

	mytag = '{'+NS['PepXML']+'}msms_run_summary'
	context = etree.iterparse(pepxml, tag=mytag)
	fast_iter(context, parseRawData)

	print('Found '+str(len(raw_files.keys()))+' raw files')
	return raw_files

def compileScanEvents(d, rdf):
	''' Determine which scans need to be extracted from each raw file '''
	scans = {}
	for rdf_short in rdf.keys():
		scans[rdf_short] = []
	
	for protein in d.keys():
		for peptide in d[protein]['peptides'].keys():
			for spectrum in d[protein]['peptides'][peptide]['spectra']:
				parts = spectrum.split('.')
				short = '.'.join(parts[0:-3])
				scan_num = int(parts[-3])
				scans[short].append(scan_num)

	for rdf_short, scanlist in scans.iteritems():
		scanlist.sort()
	return scans

def extractPeaks(rdf, scans):
	''' Extract peak information from raw files '''
	peak_info = {}
	for rdf_short, fname in rdf.iteritems():
		print('\nExtracting peak intensities from '+rdf_short)
		count = 0
		peak_info[rdf_short] = {}
		needed_scans = scans[rdf_short]
		msrun = pymzml.run.Reader(fname)
		for s in msrun:
			if s['id'] in needed_scans:
				peak_info[rdf_short][s['id']] = []
				for mz, i in s.peaks:
					peak_info[rdf_short][s['id']].append((mz, i))
				count += 1
				pct_complete = round(100*(float(count) / len(needed_scans)), 1)
				sys.stdout.write('\r'+str(pct_complete)+'% complete')

	return peak_info

def calcIons(peptide, mods):
	''' Calculate ion series for a given modified peptide '''
	theoretical_ions = []
	# Convert peptide to stripped form, handling N-terminal modifications
	if re.match(r'^n', peptide): # N-term mods are special
		parts = peptide.split('[')[1].split(']')
		val = parts[0]
		naa = parts[1][0]
		if val == '43':
			mass = md['acetylation'] + md[naa]
			mods.append((1, mass))
		elif val == '1':
			mass = md['deamidation'] + md[naa]
			mods.append((1, mass))
		else:
			print ("MAJOR ERROR: Could not decipher N-terminal mod. B-ion series will be corrupted")

	stripped = re.sub(r'[a-z]*\[[0-9]*\]','',peptide) # Strip out other mod info in string
	if not mods:
		aa_masses = [md[aa] for aa in stripped] # pull masses from global dict
	else: # Handle mods
		mod_positions = [mod[0] for mod in mods]
		mod_masses = [mod[1] for mod in mods]
		aa_masses = []
		for i in range(0,len(stripped)):
			if i+1 not in mod_positions:
				try:
					aa_masses.append(md[stripped[i]])
				except:
					print('Error assigning mass to '+peptide+' at position '+str(i))
			else:
				aa_masses.append(mod_masses[mod_positions.index(i+1)])

	series_list = ['b1','b2','y1','y2'] # Calculate +1 and +2 charged b and y ions
	
	for series in series_list:
		scharge = int(series[1])
		if series[0] == 'b':
			for i in range(1, len(aa_masses)):
				if i+1 < len(aa_masses):
					ion_mass = (sum(aa_masses[0:i+1]) + md['proton']) / scharge
					theoretical_ions.append(ion_mass)
		elif series[0] == 'y':
			y_masses = aa_masses[::-1]
			for i in range(0, len(y_masses)):
				if i+1 < len(y_masses):
					ion_mass = (sum(y_masses[0:i+1]) + md['H2O'] + md['proton']) / scharge
					theoretical_ions.append(ion_mass)
	return list(sorted(theoretical_ions)) # list

def extractIntensities(d, raw_files, scans, tol):
	''' Extract raw data, peak match, and record intensities in one shot '''
	print('Extracting intensity data')
	# Remap peptide:spectrum information as spectrum:peptide
	# Also need peptide:protein
	spec2pep = {}
	pep2prot = {}
	for protein in d.keys():
		for peptide in d[protein]['peptides'].keys():
			pep2prot[peptide] = protein
			spectra = d[protein]['peptides'][peptide]['spectra']
			for spectrum in spectra:
				parts = spectrum.split('.')
				short = '.'.join(parts[0:-3])
				scan_num = int(parts[-3])
				#short = parts[0]
				#scan_num = int(parts[1])
				try:
					spec2pep[short][scan_num] = peptide
				except:
					spec2pep[short] = {}
					spec2pep[short][scan_num] = peptide

	gi = 0 # Global matched intensity
	for short in spec2pep.keys():
		count = 0
		print('\nExtracting from '+short)
		rf = raw_files[short]
		needed_scans = scans[short]
		msrun = pymzml.run.Reader(rf)
		for s in msrun:
			if s['id'] in needed_scans:
				# Info from the scan
				peak_masses = [mz for mz, i in s.peaks]
				peak_intensities = [i for mz, i in s.peaks]

				# Now map back to the peptide + mods
				thispep = spec2pep[short][s['id']]
				prot = pep2prot[thispep]
				pepmods = d[prot]['peptides'][thispep]['modifications']

				# And generate the ions to match against
				theoretical_ions = calcIons(thispep, pepmods)

				# Now match theoretical to observed ions, get intensities and write to struct
				intensities = []
				for ti in theoretical_ions:
					pos = bisect.bisect(peak_masses, ti)
					matched_intensities = []
					left_boundary = 0
					right_boundary = 0
					if pos == 0:
						left_boundary = 1
					elif pos == len(peak_masses):
						right_boundary = 1
					i = 1
					j = 0
					while not left_boundary:
						if abs(ti - peak_masses[pos - i]) <= tol:
							matched_intensities.append(peak_intensities[pos-i])
							if pos - i == 0:
								left_boundary = 1
							i += 1
						else:
							left_boundary = 1
					while not right_boundary:
						if abs(ti - peak_masses[pos + j]) <= tol:
							matched_intensities.append(peak_intensities[pos+j])
							j += 1
							if pos + j == len(peak_masses):
								right_boundary = 1
						else:
							right_boundary = 1
					if matched_intensities:
						top_intensity = list(sorted(matched_intensities))[-1]
						if top_intensity not in intensities: # Don't double-count peaks inside tolerance
							intensities.append(top_intensity)
				try:
					d[prot]['peptides'][thispep]['intensity'] += sum(intensities)
				except:
					d[prot]['peptides'][thispep]['intensity'] = sum(intensities)
				gi += sum(intensities)
				count += 1
				pct_complete = round(100*(float(count) / len(needed_scans)), 1)
				sys.stdout.write('\r'+str(pct_complete)+'% complete')

	return d, gi

def compileIntensities(d, peak_info, tol):
	''' Here we go, time to merge the protein data structure with the raw peak information '''
	print('\nCalculating theoretical ion series for all peptides')
	gi = 0
	for protein in d.keys():
		intensities = []
		for peptide in d[protein]['peptides'].keys():
			mods = d[protein]['peptides'][peptide]['modifications']
			theoretical_ions = calcIons(peptide, mods)

			# Now get peak intensities
			spectra = d[protein]['peptides'][peptide]['spectra']
			for spectrum in spectra:
				parts = spectrum.split('.')
				short = parts[0]
				scan_num = int(parts[1])
				peaks = peak_info[short][scan_num]
				peak_masses = [peak[0] for peak in peaks]
				peak_intensities = [peak[1] for peak in peaks]

				for ti in theoretical_ions:
					pos = bisect.bisect(peak_masses, ti)
					matched_intensities = []
					left_boundary = 0
					right_boundary = 0
					if pos == 0:
						left_boundary = 1
					elif pos == len(peak_masses):
						right_boundary = 1
					i = 1
					j = 0
					while not left_boundary:
						if abs(ti - peak_masses[pos - i]) <= tol:
							matched_intensities.append(peak_intensities[pos-i])
							if pos - i == 0:
								left_boundary = 1
							i += 1
						else:
							left_boundary = 1
					while not right_boundary:
						if abs(ti - peak_masses[pos + j]) <= tol:
							matched_intensities.append(peak_intensities[pos+j])
							j += 1
							if pos + j == len(peak_masses):
								right_boundary = 1
						else:
							right_boundary = 1
					if matched_intensities:
						top_intensity = list(sorted(matched_intensities))[-1]
						if top_intensity not in intensities: # Don't double-count peaks inside tolerance
							intensities.append(top_intensity)
			d[protein]['peptides'][peptide]['intensity'] = sum(intensities)
			gi += sum(intensities)

	return d, gi

def calcNSI(d, gi):
	''' Produce the normalized spectral index '''
	print('\nCalculating Normalized Spectral Indices')
	for protein in d.keys():
		length = d[protein]['length']
		pep_intensities = []
		for peptide in d[protein]['peptides'].keys():
			try:
				pep_intensities.append(d[protein]['peptides'][peptide]['intensity'])
			except:
				print('Did not extract intensity for '+peptide+' in '+protein)
		prot_intensity = sum(pep_intensities)

		try:
			nsi = log((prot_intensity / gi / length), 2)
			d[protein]['nsi'] = nsi
		except:
			print("Error calculating SIn for "+protein)
			print(prot_intensity, gi, length)
	
	return d

def outputNSI(protxml, d):
	''' Write output '''
	outroot = protxml.split('.')[0]
	outfn = outroot + '_nsi.csv'
	print('Writing output to '+outfn)
	with open(outfn, 'wb') as f:
		headers = ['Name','Description','SIn']
		f.write(','.join(headers)+'\n')
		for protein in d.keys():
			outline = [protein, '"'+d[protein]['description']+'"', str(d[protein]['nsi'])]
			f.write(','.join(outline)+'\n')
	return True

if __name__ == '__main__':
	# Initiate and set up argument parser
	parser = argparse.ArgumentParser(
		description = 'MS2 intensity-based label-free quantification'
		)
	parser.add_argument('protXML',
		help = 'A protXML file containing identifications, scored with iProphet or PeptideProphet')
	parser.add_argument('-f', '--fdr', type=float, default=0.01,
		help = 'An FDR cutoff value, e.g. 0.01.  Default is 0.01')
	parser.add_argument('-t', '--tolerance', type=float, default=0.4,
		help = 'Mass tolerance for matching MS2 peaks (Daltons). Default is 0.4')

	args = parser.parse_args()
	protxml = args.protXML
	fdr = args.fdr
	tol = args.tolerance
	# End argument parsing

	pc = calcFDR(protxml, fdr)
	d = getProts(protxml, pc)
	d = checkLengths(d, protxml)
	pepxml = getPepXML(protxml)
	d = getSpectra(d, pepxml)
	raw_files = getRawData(pepxml)
	scans = compileScanEvents(d, raw_files)
	d, gi = extractIntensities(d, raw_files, scans, tol)
	d = calcNSI(d, gi)
	outputNSI(protxml, d)
	print("--- %s seconds ---" % (time.time() - start_time))
