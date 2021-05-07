import pickle
import pandas as pd
import sys,re


# Will have to define cdr3 and contignt as def *char (perhaps str)
def check_flanking_cdr3(cdr3, contignt):
    '''
    function check flanking conserved C and F(W)
    '''
    try:
        start = re.search(cdr3, contignt).start()-3
        end   = re.search(cdr3, contignt).end()
    except:
        sys.stderr.write("{} not flanked by C and F in {}".format(cdr3, contignt))
        return False

    if contignt[start:start+3] in ['TGT','TGC'] and contignt[end:end+3] in ['TTT','TTC','TGG']:
        return True    
    
    return False


def parse_header(line):
    '''
    Functon parses line where e.g. "Length=437"
    is returned as a dict {'Length':'437'}
    line should NOT include "\n"
    '''

    pos = [0] + [i.start()+1 for i in re.finditer("[|=]", line)]
    lpos = len(pos)-1

    eqs = [line[pos[idx]:pos[idx+1]-1] for idx in range(0,lpos)] + [line[pos[lpos]:]]
    results = {eqs[n].lower():eqs[n+1] for n in range(1,len(eqs)-1,2)}
    results['UMI'] = eqs[0][1:]
    
    return  results


def igblast_parse(fin):
    # results table
    d={}

    # define all column names THAT ARE NOT IN THE HEADER!
    column_names = ['Chain type', 'D region', 'D-J junction','E_value', 'J start', 'Productive', 'Strand',
                    'Top D gene match', 'Top J gene match', 'Top V gene match', 'V end', 'V-D junction',
                    'V-J frame', 'V-J junction', 'cdr3', 'cdr3aa', 'cdr3end', 'cdr3start','stop codon', 
                    'conserved_cfw']
    column_names = [i.replace(" ","_").replace("-","_").lower() for i in column_names]

    # This is needed to include extra information in separate lines
    rearrangement_flag, E_flag, header_flag, junction_flag, cdr3_flag, summary_flag = [False]*6
    
    for line in fin:
        # skip empty rows
        if len(line) < 2: continue
        try:
            # start new item
            if line[:5] == "Query":
                header= line.strip()[7:]
                header_flag=True
                # process the last one!

                if 'Umi' in locals():
                    if 'cdr3' not in d[Umi].keys() or 'contignt' not in d[Umi].keys():
                        d[Umi]['conserved_cfw'] = False 
                    else:
                        d[Umi]['conserved_cfw'] = check_flanking_cdr3(d[Umi]['cdr3'], d[Umi]['contignt'])


                print(Umi, "*************************")

       # this finishes the header including the length
            elif line[:6] == "Length":
                header += '|' + line.strip()
                header_flag=False

                # extract info from header
                tmp_d = parse_header(header)
                Umi = tmp_d.pop('UMI')
                # Depending on the fastq used, UMI could be repeated
                # HERE we assume that UMIs should not be repeated
                if Umi in d:
                    sys.stderr.write('UMI={} is repeated. stick with the first one found in the data\n'.format(Umi))
                    continue

                # initialize this item of d to avoid d-items with different number of key-value items                
                d[Umi] = {i:"N/A" for i in column_names}

                # line is now Length
                for tmpk, tmpv in tmp_d.items():
                    d[Umi][tmpk.replace(" ","_").replace("-","_").lower()] = tmpv
                #d[UMI] = tmp_d
                d[Umi]['contignt'] = ''
                continue

            # the header inifo is split into 2 lines
            elif header_flag:
                header += line.strip().lower()
            
            # stop here if not done with header
            #if 'header' not in locals() or header_flag: continue

            #------------------------------------------------------
            # all extra components that are not in the header 
            if line[:8] == "Sequence":
                E_flag=True
                continue

            # for now we only include the info that was included in regeneron analysis and later we could 
            # decide on including extra information if usefill for anuthing => Only taking the first V-(D)-J line
            elif line[:21] == 'V-(D)-J rearrangement': 
                rearrangement_flag = True
                Keys = re.search("query sequence \((.*)\).",line).group(1).split(',')
                Keys = [i[1:].replace(" ","_").replace("-","_") if i[0]==" " else i for i in Keys]
                continue

            elif line[:16] == 'V-(D)-J junction':
                junction_flag = True
                Keys = re.search("matches \((.*)\)\.",line).group(1).split(',')
                Keys = [i[1:].replace(" ","_").replace("-","_") if i[0]==" " else i for i in Keys]
                continue

            elif line[:19] == 'Sub-region sequence':
                cdr3_flag = True
                continue

            elif line[:17] == 'Alignment summary':
                summary_flag = True
                continue

            # knowing the flags, proceed with updating the values in d[UMI]
            elif re.search('Query_', line):
                d[Umi]['contignt'] +=  re.search('[ACTGRYSWKMBDHVN]{1,}',line).group()
                continue

            if E_flag:
                d[Umi]['e_value'] = line.strip().split(' ')[-1]
                E_flag=False

            elif rearrangement_flag:
                # replace "," for ";" to use "," in ngs_agg to split line in some values
                for i,j in zip(map(lambda x:x.lower().replace(" ","_").replace("-","_"), Keys), 
                               line.strip().replace(",",";").split('\t')): 
                
                    d[Umi][i] = j
                rearrangement_flag=False

            if junction_flag:
                for i,j in zip(map(lambda x:x.lower().replace(" ","_").replace("-","_"), Keys), line.strip().split('\t')):
                    d[Umi][i]=j
                junction_flag=False

            elif cdr3_flag:
                for i,j in zip(['cdr3', 'cdr3aa', 'cdr3start', 'cdr3end'],line.strip().split('\t')[1:]):
                    d[Umi][i]=j
                cdr3_flag=False

        except Exception as e:
            sys.stderr.write('ERROR > {}\n'.format(str(e)))

    fin.close()
    
    return d


def main():
	from signal import signal, SIGPIPE, SIG_DFL
	signal(SIGPIPE,SIG_DFL) 


	if len({'-help','--help'}.intersection(set(sys.argv)))>0:
		sys.exit('''\n 
						igblast_output | python <this script> <options> --out <out_filename_prefix>\n
						python <this script <options> --in <igblast_output> --out <out_filename_prefix>\n
						options: --seq-len (if the aligned sequence is longer than N nucleotides, trim it to "NNN"\n
						and let the user know about that in the std err.)\n
				 ''')
	seq_len_filter = 1000

	# initialize prefix filename or get from user 
	out_prefix = 'igblast_output'
	fin = sys.stdin
	for n,i in enumerate(sys.argv):
		if i in ['--out', '-out']:
			out_prefix = sys.argv[n+1]
		if i in ['--in', '-in']:
			fin = open(sys.argv[n+1])
		if i in ['--seq-len', 'filter-len']:
			seq_len_filter = int(sys.argv[n+1])

	d = igblast_parse(fin)

	# if contignt is too long, apply the filter. Otherwise, galaxy will throw and error
	for umi in d.keys(): 
		seq = d[umi]['contignt']
		if len(seq) > seq_len_filter: 
			d[umi]['contignt'] = 'NNN' 
			sys.stderr.write('sequence starting with {} is too long ({} bp) it will be replaced in the output with NNN\n'.format(\
							 seq[:20], len(seq)))

	pd.DataFrame(d).T.to_csv(path_or_buf=out_prefix + '.csv', index_label="umi")


if __name__ == '__main__': 
	main()
