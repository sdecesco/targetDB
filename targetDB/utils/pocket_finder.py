#!/usr/bin/env python
import os
import subprocess

import pandas as pd
from operator import itemgetter
from targetDB.utils import pdb_parser as parser
from Bio.PDB import *
from pathlib import Path

Bioparser = PDBParser(PERMISSIVE=1, QUIET=True)
io = PDBIO()


def fpocket_launcher(pdb, pdb_file, sphere_size_min=3.0,verbose=False):
	if verbose:
		print('[POCKET SEARCH]: ' + str(pdb_file))
	subprocess.check_output([fpocket_path, '-f', str(pdb)])


def parse_pockets(out_file):
	with open(out_file, 'r') as pockets_list:
		pockets_dict = {}
		pocket_parameters = {}
		pocket_number = None
		for line in pockets_list:
			line = line.replace('\t', '').replace('\n', '')
			if line.startswith('Pocket'):
				line = line.replace('Pocket ', 'p').replace(' ', '').split(':')
				pocket_number = line[0]
			else:
				line = line.split(':')
				if line[0] == '':
					pockets_dict[pocket_number] = pocket_parameters
					pocket_parameters = {}
				else:
					param = line[0]
					param_value = line[1]
					param = param.strip(' ').replace('-', '').replace('.', '').replace(' ', '_').lower()
					param_value = param_value.strip(' ').replace(' ', '_')
					pocket_parameters[param] = param_value
	return pockets_dict


def get_object(p_dict, pdb_name, out_pockets_dir_path, pdb_info, domain):
	pockets_objects = {}
	for key, values in p_dict.items():
		pockets_objects[key] = Pockets(values, pdb_name, key, out_pockets_dir_path, pdb_info, domain)
	return pockets_objects


class ChainSelect(Select):
	def __init__(self, list_of_chains):
		self.list = list_of_chains

	def accept_chain(self, chain):
		if chain.id in self.list:
			return True
		else:
			return False


def get_pockets(path, sphere_size=3.0, pdb_info=pd.DataFrame(), domain=pd.DataFrame(), alternate=False, alternate_pdb=pd.DataFrame(),
				uniprot_id=None,fpocket_exe=None,verbose=False):
	global fpocket_path
	fpocket_path = fpocket_exe


	pock = pd.DataFrame(columns=['PDB_code','Target_id','Pocket_number','Pocket_id','Score','DrugScore','apolar_sasa','polar_sasa','total_sasa','volume','blast','druggable'])
	pock_chain = pd.DataFrame(columns=['Pocket_id','Chain_id','List_of_contacts'])
	pock_domain = pd.DataFrame(columns=['Pocket_id','Domain_id','Coverage'])

	if os.name == 'nt':
		return {'pockets': pock, 'pockets_chain': pock_chain, 'pockets_domain': pock_domain}


	files = []
	if not pdb_info.empty:
		for pdb_code in pdb_info.index:
			if Path(pdb_info.path.loc[pdb_code]).is_file():
				files.append(pdb_info.loc[pdb_code]['path'])
	if not alternate_pdb.empty:
		group = alternate_pdb.groupby('PDB_code')
		c_list = group.chain_letter.apply(list)
		p_list = group.path.apply(set).apply(list).apply(lambda x: ''.join(x))
		alternate_pdb = pd.DataFrame({'chain_letter': c_list, 'path': p_list})
		for pdb_code in alternate_pdb.index:
			if Path(alternate_pdb.path.loc[pdb_code]).is_file():
				files.append(alternate_pdb.loc[pdb_code]['path'])
	results = {}
	for f in files:

		pdb_file = Path(f).name
		pdb_code = str(pdb_file).rstrip('.pdb')

		pocket_path = path.joinpath('POCKETS')
		if not pocket_path.is_dir():
			pocket_path.mkdir(parents=True, exist_ok=True)

		pdb_parsed = parser.parse_header(pdb_code, f)

		pdb_strip_path = pocket_path.joinpath(pdb_code + '_ALL.pdb')
		out_path = pocket_path.joinpath(pdb_code + '_ALL_out')
		out_file_path = out_path.joinpath(pdb_code + '_ALL_info.txt')
		out_pockets_dir_path = out_path.joinpath('pockets')

		if pdb_parsed.biomolecules:
			chain_to_keep = []
			if not pdb_info.empty:
				for biomolecule,chains in pdb_parsed.biomolecules.items():
					if not chain_to_keep:
						for c in pdb_info.loc[pdb_code]['Chain']:
							if c in chains:
								chain_to_keep = pdb_parsed.biomolecules[biomolecule]
								pdb_strip_path = pocket_path.joinpath(pdb_code +'_BIO'+biomolecule+'.pdb')
								out_path = pocket_path.joinpath(pdb_code +'_BIO'+biomolecule+'_out')
								out_file_path = out_path.joinpath(pdb_code +'_BIO'+biomolecule+ '_info.txt')
								out_pockets_dir_path = out_path.joinpath('pockets')
								break
			elif not alternate_pdb.empty:
				for biomolecule,chains in pdb_parsed.biomolecules.items():
					if not chain_to_keep:
						for c in alternate_pdb.loc[pdb_code]['chain_letter']:
							if c in chains:
								chain_to_keep = pdb_parsed.biomolecules[biomolecule]
								pdb_strip_path = pocket_path.joinpath(pdb_code + '_BIO' + biomolecule + '.pdb')
								out_path = pocket_path.joinpath(pdb_code + '_BIO' + biomolecule + '_out')
								out_file_path = out_path.joinpath(pdb_code + '_BIO' + biomolecule + '_info.txt')
								out_pockets_dir_path = out_path.joinpath('pockets')
								break
			else:
				chain_to_keep = []

			structure = Bioparser.get_structure(pdb_code, f)
			io.set_structure(structure[0])
			if chain_to_keep:
				chain_select = ChainSelect(chain_to_keep)
				io.save(str(pdb_strip_path), chain_select)
			else:
				io.save(str(pdb_strip_path))
		if alternate:
			info = pd.DataFrame()
		elif not pdb_info.empty:
			info = pdb_info.loc[pdb_code]
		else:
			info = pd.DataFrame()

		if out_path.is_dir():
			if out_file_path.is_file():
				if out_pockets_dir_path.is_dir():
					p_out = len([x for x in out_pockets_dir_path.iterdir()])
					if p_out!=0:
						if verbose:
							print('[POCKET ALREADY EXISTS]: ' + pdb_code)
						pockets = parse_pockets(str(out_file_path))
						results[pdb_code] = get_object(pockets, pdb_code, str(out_pockets_dir_path), info,
													   domain)
						continue
					else:
						fpocket_launcher(pdb_strip_path, pdb_file, sphere_size_min=sphere_size,verbose=verbose)
				else:
					fpocket_launcher(pdb_strip_path, pdb_file, sphere_size_min=sphere_size,verbose=verbose)
			else:
				fpocket_launcher(pdb_strip_path, pdb_file, sphere_size_min=sphere_size,verbose=verbose)
		else:
			fpocket_launcher(pdb_strip_path, pdb_file, sphere_size_min=sphere_size,verbose=verbose)
		if out_path.is_dir():
			if out_file_path.is_file():
				if out_pockets_dir_path.is_dir():
					p_out = len([x for x in out_pockets_dir_path.iterdir()])
					if p_out != 0:
						if verbose:
							print('[POCKET SEARCH DONE]: ' + pdb_code)
						pockets = parse_pockets(str(out_file_path))
						results[pdb_code] = get_object(pockets, pdb_code, str(out_pockets_dir_path), info,
													   domain)
					else:
						print(
							'[ERROR]: No output file detected (problem with PDB file or no pocket found): ' + pdb_file)
						print('[ERROR]: Script will now continue ...')
				else:
					print('[ERROR]: No output file detected (problem with PDB file or no pocket found): ' + pdb_file)
					print('[ERROR]: Script will now continue ...')
			else:
				print('[ERROR]: No output file detected (problem with PDB file or no pocket found): ' + pdb_file)
				print('[ERROR]: Script will now continue ...')
		else:
			print('[ERROR]: No output file detected (problem with PDB file or no pocket found): ' + pdb_file)
			print('[ERROR]: Script will now continue ...')



	if results:
		pock = pd.DataFrame(columns=['PDB_code','Target_id','Pocket_number','Pocket_id','Score','DrugScore','apolar_sasa','polar_sasa','total_sasa','volume','blast','druggable'])
		pock_chain = pd.DataFrame(columns=['Pocket_id','Chain_id','List_of_contacts'])
		pock_domain = pd.DataFrame(columns=['Pocket_id','Domain_id','Coverage'])
		for pdb_code,pockets in results.items():
			for p in pockets.keys():
				pocket_id = uniprot_id+'_'+pdb_code+'_'+pockets[p].pocket_number
				if alternate:
					blast = 'TRUE'
				else:
					blast = 'FALSE'
				pock.loc[pocket_id] = {'PDB_code':pdb_code,'Target_id':uniprot_id,'Pocket_number':pockets[p].pocket_number,'Pocket_id':pocket_id,'Score':pockets[p].score,'DrugScore':pockets[p].druggability_score,'apolar_sasa':pockets[p].apolar_sasa,'polar_sasa':pockets[p].polar_sasa,'total_sasa':pockets[p].total_sasa,'volume':pockets[p].volume,'blast':blast,'druggable':pockets[p].druggable}

				if pockets[p].chain_coverage:
					for c in pockets[p].chain_coverage.keys():
						pock_chain.loc[pocket_id+'_'+c] = {'Pocket_id':pocket_id,'Chain_id':pdb_code + '_' + c,'List_of_contacts':pockets[p].chain_coverage[c]}
				if pockets[p].part_of_domain:
					for d in pockets[p].part_of_domain:
						pock_domain.loc[pocket_id + '_' + d['domain_id']] = {'Pocket_id':pocket_id,'Domain_id':d['domain_id'],'Coverage':d['coverage']}

	return {'pockets':pock,'pockets_chain':pock_chain,'pockets_domain':pock_domain}


def get_bests(dict_pdbs, n_pockets=4, rank_by='druggability_score'):
	rank_list = {}
	best_pockets = {}
	rank_list['other'] = []
	best_pockets['other'] = []
	for k in dict_pdbs.keys():
		for p in dict_pdbs[k].keys():
			for d in dict_pdbs[k][p].part_of_domain:
				rank_list[d['domain']] = []
				best_pockets[d['domain']] = []

	for k in dict_pdbs.keys():
		for p in dict_pdbs[k].keys():
			in_list = False
			for d in dict_pdbs[k][p].part_of_domain:
				if d['coverage'] > 5.0:
					in_list = True
					rank_list[d['domain']].append((k, p, float(dict_pdbs[k][p].param[rank_by])))
			if not in_list:
				rank_list['other'].append((k, p, float(dict_pdbs[k][p].param[rank_by])))

	for key in rank_list.keys():
		data = sorted(rank_list[key], key=itemgetter(2), reverse=True)
		best = data[:n_pockets]
		for entry in best:
			best_pockets[key].append(dict_pdbs[entry[0]][entry[1]])
	return best_pockets


def get_druggable_pockets(list_of_pockets):
	"""Get a dict of PDB with a dict of Pockets object as entry, output a list of druggable pockets"""
	new_dict = {}
	for k in list_of_pockets.keys():
		new_dict[k] = {}
		for p in list_of_pockets[k].keys():
			if list_of_pockets[k][p].druggable == 'TRUE':
				new_dict[k][p] = list_of_pockets[k][p]
	key_to_del = []
	for k in new_dict.keys():
		if not new_dict[k]:
			key_to_del.append(k)
	for k in key_to_del:
		new_dict.pop(k)
	return new_dict


def show_pockets(pocket_dict):
	"""get a pocket_dict from the get_pockets() function"""
	for v in pocket_dict.values():
		for obj in v.values():
			obj.show()


def show_pockets_sorted(pocket_dict):
	"""get a pocket_dict from the get_pockets() function"""
	p = get_bests(pocket_dict)
	for pocket in p:
		pocket.show()


def open_pymol(path_to_pml):
	subprocess.call(['/usr/local/scripts/pymol-x64-COS7', str(path_to_pml)])


class Pockets:
	"""
	Attributes available for pocket objects:

	self.score
	self.druggability_score
	self.number_of_alpha_spheres
	self.total_sasa
	self.polar_sasa
	self.apolar_sasa
	self.volume
	self.mean_local_hydrophobic_density
	self.mean_alpha_sphere_radius
	self.mean_alp_sph_solvent_access
	self.apolar_alpha_sphere_proportion
	self.hydrophobicity_score
	self.volume_score
	self.polarity_score
	self.charge_score
	self.proportion_of_polar_atoms
	self.alpha_sphere_density
	self.cent_of_mass__alpha_sphere_max_dist
	self.flexibility
	self.pdb_name
	self.pocket_number

	"""

	def __init__(self, values, pdb_name, pnumber, out_pockets_dir_path, pdb_info, domain):
		self.druggable = 'FALSE'
		self.param = {}

		for k, v in values.items():
			self.param[k] = v
			setattr(self, k, v)
		self.pdb_name = pdb_name
		pnumber = 'p' + str(int(str(pnumber).lstrip('p')))
		self.pocket_number = pnumber
		self.pocket_path = out_pockets_dir_path
		if not pdb_info.empty and not domain.empty:
			self.get_residues(domain, pdb_info)
		else:
			self.chain_coverage = {}
			self.part_of_domain = []
		if float(self.param['druggability_score']) >= 0.5 and float(self.param['volume']) >= 250.0:
			self.druggable = 'TRUE'

	def show(self):
		print('==================================================')
		print('pdb_name: ' + str(self.pdb_name))
		print('==================================================')
		print('Name: ' + str(self.pocket_number))
		print('Score: ' + str(self.score))
		print('Druggability score: ' + str(self.druggability_score))
		print('Volume: ' + str(self.volume))
		print('Area: ' + str(self.total_sasa))
		print('==================================================')

	def get_residues(self, domain, pdb_info):
		pocket_number = str(self.pocket_number).lstrip('p')
		file = Path(self.pocket_path).joinpath('pocket' + pocket_number + '_atm.pdb')

		pocket = Bioparser.get_structure(self.pocket_number, str(file))
		data = {}
		for chain in pocket[0]:
			data[str(chain.id).strip(' ')] = []
			for residue in chain:
				data[str(chain.id).strip(' ')].append(str(residue.id[1]).strip(' '))
		self.chain_coverage = data
		self.part_of_domain = []
		domain_coverage = {}
		for d in domain.index:
			domain_coverage[domain.loc[d]['name']] = 0
		for k in self.chain_coverage.keys():
			if k in pdb_info['Chain']:
				for res_number in self.chain_coverage[k]:
					for d in domain.index:
						if int(domain.loc[d]['Start']) <= int(res_number) <= int(domain.loc[d]['Stop']):
							domain_coverage[domain.loc[d]['name']] += 1
		for key in domain_coverage.keys():
			if domain_coverage[key] != 0:
				for d in domain.index:
					if key == domain.loc[d]['name']:
						try:
							self.part_of_domain.append(
								{'domain': key, 'coverage': round((domain_coverage[key] / domain.loc[d]['length']) * 100, 0),
								 'domain_id': domain.loc[d]['domain_id']})
						except KeyError:
							self.part_of_domain.append(
								{'domain': key, 'coverage': round((domain_coverage[key] / domain.loc[d]['length']) * 100, 0)})
		if not self.part_of_domain:
			count = 0
			for k in self.chain_coverage.keys():
				if k in pdb_info['Chain']:
					count += len(self.chain_coverage[k])
			try:
				if pdb_info['length'] == 0:
					self.part_of_domain.append(
						{'domain': 'other', 'coverage': '', 'domain_id': 'other'})
				else:
					self.part_of_domain.append(
						{'domain': 'other', 'coverage': round((count / pdb_info['length']) * 100, 0),
						 'domain_id': 'other'})
			except KeyError:
				self.part_of_domain.append(
					{'domain': 'other', 'coverage': '', 'domain_id': 'other'})
