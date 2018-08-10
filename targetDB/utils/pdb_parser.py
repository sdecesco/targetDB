#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import re
from Bio.PDB import *
from collections import OrderedDict

Bioparser = PDBParser(PERMISSIVE=1, QUIET=True)


class parsed_pdb:
	def __init__(self, pdb_code):
		self.code = pdb_code
		self.biomolecules = []


def parse_header(pdb_code, pdb_file):
	pdb_parsed = parsed_pdb(pdb_code)
	with open(pdb_file, 'r') as pdb:
		biomolecule = {}
		for line in pdb:
			if 'REMARK 350' in line:
				if line.startswith('REMARK 350'):
					current_line = line.replace("REMARK 350 ", '').strip('\n').strip(' ')
					if current_line.startswith('BIOMOLECULE:'):
						biomol_id = current_line.split(':')[1].strip(' ')
					if current_line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
						biomol_chains = current_line.split(':')[1].replace(' ', '').rstrip(',').split(',')
						biomolecule[biomol_id] = biomol_chains
					if current_line.startswith('AND CHAINS: '):
						biomol_chains = current_line.split(':')[1].replace(' ', '').rstrip(',').split(',')
						biomolecule[biomol_id].extend(biomol_chains)


		pdb_parsed.biomolecules = OrderedDict(sorted(biomolecule.items(),key=lambda t:t[0]))
	return pdb_parsed


def get_sequence(pdb_code,pdb_file, chains_to_do, domains):
	three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
					'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
					'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
					'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
					'SEC': 'U'}
	modified_residue = {'MSE': 'MET'}

	with open(pdb_file, 'r') as pdb:
		sequence = {}
		list_of_chain = []
		chain_stat = {}
		for line in pdb:
			if 'ATOM' in line:
				if line.startswith('ATOM'):
					break
			elif 'DBREF' in line:
				if line.startswith('DBREF'):
					current_line = list(filter(None, line.strip('\n').replace(' ', '_').split('_')))
					if current_line[5] == 'PDB':
						continue
					if current_line[2] not in chain_stat.keys():
						chain_stat[current_line[2]] = {'start': 0,'stop': 0,'start_stop':[]}
					if current_line[0] == 'DBREF1':
						continue
					elif current_line[0] =='DBREF2':
						chain_stat[current_line[2]]['start_stop'].append((int(re.sub("\D", "", current_line[-2])),int(re.sub("\D", "", current_line[-1]))))
					else:
						chain_stat[current_line[2]]['start_stop'].append((int(re.sub("\D", "", current_line[-2])), int(re.sub("\D", "", current_line[-1]))))
			elif 'SEQRES' in line:
				if line.startswith('SEQRES'):
					current_line = list(
						filter(None, line.strip('\n').replace('   ', ' ').replace('  ', ' ').replace(' ', '_').split('_')))
					chain = current_line[2]
					if chain not in list_of_chain:
						list_of_chain.append(chain)
						sequence[chain] = ''
					seq = current_line[4:]
					for res in seq:
						try:
							r = three_to_one[res]
						except KeyError:
							try:
								r = three_to_one[modified_residue[res]]
							except KeyError:
								try:
									with open(pdb_file, 'r') as pdb_2:
										for modres in pdb_2:
											if 'MODRES' in modres:
												if modres.startswith('MODRES'):
													mod_line = list(filter(None, modres.strip('\n').split(' ')))
													modified_residue[mod_line[2]] = mod_line[5]
											elif 'ATOM' in modres:
												if modres.startswith('ATOM'):
													break
									r = three_to_one[modified_residue[res]]
								except KeyError:
									r = ''
						sequence[chain] += r
		for k,v in chain_stat.items():
			v['start'] = v['start_stop'][0][0]
			v['stop'] = v['start_stop'][-1][1]
		list_seq_df = pd.DataFrame(columns=['PDB_code','chain_name', 'sequence', 'equal','start', 'stop','length','start_stop_pairs','domain','domain_id','seq_list'])
		for key, values in sequence.items():
			if key in chains_to_do:
				if key in chain_stat.keys():
					list_seq_df.loc[pdb_code+'_'+key] = {'PDB_code':pdb_code,'chain_name': key, 'sequence': values, 'equal': [],'start': chain_stat[key]['start'], 'stop': chain_stat[key]['stop'],'length': len(values),'start_stop_pairs':chain_stat[key]['start_stop'],'domain':[None],'domain_id':[None],'seq_list':[None]}
				else:
					list_seq_df.loc[pdb_code+'_'+key] = {'PDB_code':pdb_code,'chain_name': key, 'sequence': values, 'equal': [],'start': 1, 'stop': len(values),'length': len(values),'start_stop_pairs':[],'domain':[None],'domain_id':[None],'seq_list':[None]}

		for keyA in list_seq_df.index:
			for keyB in list_seq_df.index:
				if keyA == keyB:
					continue
				elif list_seq_df.loc[keyA].sequence == list_seq_df.loc[keyB].sequence:
					list_seq_df.loc[keyA]['equal'].append(list_seq_df.loc[keyB].chain_name)
					list_seq_df.loc[keyB]['equal'].append(list_seq_df.loc[keyA].chain_name)
		for key in list_seq_df.index:
			domain_name = []
			domain_id = []
			for domain in domains.index:
				for pair in list_seq_df.loc[key]['start_stop_pairs']:
					try:
						if pair[0] >= domains.loc[domain]['Stop']:
							continue

						elif pair[0] <= domains.loc[domain]['Start']:
							if pair[1] >= domains.loc[domain]['Stop']:
								domain_name.append(domains.loc[domain]['name'])
								try:
									domain_id.append(domains.loc[domain]['domain_id'])
								except KeyError:
									pass
								continue
							elif pair[1] <= domains.loc[domain]['Start']:
								continue
							else:
								ratio = (float(pair[1]) - float(domains.loc[domain]['Start'])) / (
									float(domains.loc[domain]['Stop']) - float(domains.loc[domain]['Start']))
								if ratio > 0.3:
									domain_name.append(domains.loc[domain]['name'])
									try:
										domain_id.append(domains.loc[domain]['domain_id'])
									except KeyError:
										pass
									continue
						else:
							if pair[1] < domains.loc[domain]['Stop']:
								ratio = (float(pair[1]) - float(pair[0])) / (
									float(domains.loc[domain]['Stop']) - float(domains.loc[domain]['Start']))
								if ratio > 0.3:
									domain_name.append(domains.loc[domain]['name'])
									try:
										domain_id.append(domains.loc[domain]['domain_id'])
									except KeyError:
										pass
									continue
							else:
								ratio = (float(domains.loc[domain]['Stop']) - float(pair[0])) / (
									float(domains.loc[domain]['Stop']) - float(domains.loc[domain]['Start']))
								if ratio > 0.3:
									domain_name.append(domains.loc[domain]['name'])
									try:
										domain_id.append(domains.loc[domain]['domain_id'])
									except KeyError:
										pass
									continue
					except TypeError:
						pass

			if len(domain_name) == 0:
				domain_name = ['']
			seq_list = list(list_seq_df.sequence.loc[key])
			list_seq_df.at[key,'domain'] = domain_name
			list_seq_df.at[key,'domain_id'] = domain_id
			list_seq_df.at[key,'seq_list'] = seq_list
	return list_seq_df


if __name__ == "__main__":
	pass
