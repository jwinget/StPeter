import csv
import glob

infiles = glob.glob('10*_nsi.csv')
#infiles = ['Healthy_nsi.csv','SLS_nsi.csv','Nonlesional_nsi.csv','Lesional_nsi.csv']
print(infiles)

data = {}
descs = {}
i = 0
for item in infiles:
	exp = item.split('_')[0]
	with open(item, 'rb') as f:
		reader = csv.reader(f)
		reader.next() # Skip header
		for row in reader:
			accession = row[0]
			desc = row[1]
			nsi = row[2]

			if accession not in descs.keys():
				descs[accession] = desc
				data[accession] = ['' for item in infiles]
			data[accession][i] = nsi
	i += 1

outfile = 'nsi_merged.csv'

with open(outfile, 'wb') as f:
	headers = ['Accession', 'Description']
	for item in infiles:
		headers.append(item)
	f.write(','.join(headers)+'\n')
	for d, nsilist in data.iteritems():
		outrow = []
		outrow.append(d)
		outrow.append('"'+descs[d]+'"')
		for v in nsilist:
			outrow.append(v)
		f.write(','.join(outrow)+'\n')

print('done')

