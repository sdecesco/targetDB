#!/usr/bin/env python
import glob
import os
import sys
import subprocess
import shutil
import argparse
from operator import itemgetter
from druggability3 import pdb_parser as parser
from Bio.PDB import *

Bioparser = PDBParser(PERMISSIVE=1, QUIET=True)
io = PDBIO()


def fpocket_launcher(pdb, pdb_file, sphere_size_min=3.0):
    print('[POCKET SEARCH]:  (' + str(pdb_file) + ')')
    subprocess.check_output(
        ['/ssddata/sdecesco/fpocket2/bin/fpocket', '-f', str(pdb), '-m', str(sphere_size_min), '-i', str(35), '-r',
         str(4.5), '-M', str(5)])


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


def get_pockets(path, sphere_size=3.0, pdb_info=None, domain=None, alternate=False, alternate_pdb=None,
                uniprot_id=None):
    files = []
    if pdb_info:
        for key in pdb_info.keys():
            if 'path' in pdb_info[key].keys():
                files.append(pdb_info[key]['path'])
    if alternate_pdb:
        for key in alternate_pdb.keys():
            if 'path' in alternate_pdb[key].keys():
                files.append(alternate_pdb[key]['path'])
    results = {}
    for f in files:

        names = str(f.split('/')[
                        -1])  # names[0] = outputs / names[1] = gene name / name [2] = PDB/PDB_BLAST or BLAST_XXsimilarity_GENE / name[3] = PDB/PDB_BLAST or pdb file name / name[4] = pdb_file_name
        pdb_file = names
        pdb_code = str(pdb_file).rstrip('.pdb')

        if uniprot_id:
            pocket_path = path + 'POCKETS'
            if not os.path.exists(pocket_path):
                os.makedirs(pocket_path)

            out_path = pocket_path + '/' + str(pdb_code) + '_' + str(uniprot_id) + '_out'
            out_file_path = out_path + '/' + pdb_code + '_' + str(uniprot_id) + '_info.txt'
            out_pockets_dir_path = out_path + '/pockets'
            pdb_parsed = parser.parse_header(pdb_code, f)
            pdb_strip_path = pocket_path + '/' + str(pdb_code) + '_' + str(uniprot_id) + '.pdb'
        else:
            out_path = str(f).rstrip('.pdb') + '_strip_out'
            out_file_path = out_path + '/' + pdb_code + '_strip_info.txt'
            out_pockets_dir_path = out_path + '/pockets'
            pdb_parsed = parser.parse_header(pdb_code, f)
            pdb_strip_path = str(f).rstrip('.pdb') + '_strip.pdb'
        if pdb_parsed.biomolecules:
            chain_to_keep = pdb_parsed.biomolecules['1']
            structure = Bioparser.get_structure(pdb_code, f)
            io.set_structure(structure[0])
            chain_select = ChainSelect(chain_to_keep)
            io.save(pdb_strip_path, chain_select)
        if alternate:
            info = None
        elif pdb_info:
            info = pdb_info[pdb_code]
        else:
            info = None

        if os.path.exists(out_path):
            if os.path.isfile(out_file_path):
                if os.path.exists(out_pockets_dir_path):
                    if os.listdir(out_pockets_dir_path):
                        print('[POCKET ALREADY EXISTS]: ' + pdb_code)
                        pockets = parse_pockets(out_file_path)
                        results[pdb_code] = get_object(pockets, pdb_code, out_pockets_dir_path, info,
                                                       domain)
                        continue
                    else:
                        fpocket_launcher(pdb_strip_path, pdb_file, sphere_size_min=sphere_size)
                else:
                    fpocket_launcher(pdb_strip_path, pdb_file, sphere_size_min=sphere_size)
            else:
                fpocket_launcher(pdb_strip_path, pdb_file, sphere_size_min=sphere_size)
        else:
            fpocket_launcher(pdb_strip_path, pdb_file, sphere_size_min=sphere_size)
        if os.path.exists(out_path):
            if os.path.isfile(out_file_path):
                if os.path.exists(out_pockets_dir_path):
                    if os.listdir(out_pockets_dir_path):
                        print('[POCKET SEARCH DONE]: ' + pdb_code)
                        pockets = parse_pockets(out_file_path)
                        results[pdb_code] = get_object(pockets, pdb_code, out_pockets_dir_path, info,
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
    return results


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
        pnumber = 'p' + str(int(str(pnumber).lstrip('p')) - 1)
        self.pocket_number = pnumber
        self.pocket_path = out_pockets_dir_path
        if pdb_info and domain:
            self.get_residues(domain, pdb_info)
        else:
            self.chain_coverage = {}
            self.part_of_domain = []
        if float(self.param['druggability_score']) >= 0.5 and (
                float(self.param['apolar_sasa']) / float(self.param['total_sasa'])) > 0.5 and float(
            self.param['volume']) >= 250.0 and float(self.param['total_sasa']) >= 50.0:
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
        file = self.pocket_path + '/pocket' + pocket_number + '_atm.pdb'

        pocket = Bioparser.get_structure(self.pocket_number, file)
        data = {}
        for chain in pocket[0]:
            data[str(chain.id).strip(' ')] = []
            for residue in chain:
                data[str(chain.id).strip(' ')].append(str(residue.id[1]).strip(' '))
        self.chain_coverage = data
        self.part_of_domain = []
        domain_coverage = {}
        for d in domain:
            domain_coverage[d['name']] = 0
        for k in self.chain_coverage.keys():
            if k in pdb_info['Chain']:
                for res_number in self.chain_coverage[k]:
                    for d in domain:
                        if int(d['Start']) <= int(res_number) <= int(d['Stop']):
                            domain_coverage[d['name']] += 1
        for key in domain_coverage.keys():
            if domain_coverage[key] != 0:
                for d in domain:
                    if key == d['name']:
                        try:
                            self.part_of_domain.append(
                                {'domain': key, 'coverage': round((domain_coverage[key] / d['length']) * 100, 0),
                                 'domain_id': d['domain_id']})
                        except KeyError:
                            self.part_of_domain.append(
                                {'domain': key, 'coverage': round((domain_coverage[key] / d['length']) * 100, 0)})
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


if __name__ == "__main__":
    sys.path.append('/ssddata/sdecesco/data/Scripts/Druggability_3')
    from druggability3 import druggability_assesment as drug

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str, help='output folder name', metavar='')
    parser.add_argument('-i', '--in_file', help='Name of the input pdb file', metavar='')
    parser.add_argument('-p', '--get_pdb',
                        help='enter a single PDB code (Note: if -p option used --output option is required', metavar='')
    parser.add_argument('-S', '--show', help='display results in the console', action='store_true', default=False)
    parser.add_argument('-D', '--show_druggable', help='display druggable pockets', action='store_true', default=False)
    parser.add_argument('-P', '--open_pymol', help='Open at the end', action='store_true', default=False)
    args = parser.parse_args()

    if args.output:
        if os.path.exists(args.output):
            pass
        else:
            os.makedirs(args.output)

    if args.in_file:
        if os.path.isfile(args.in_file):
            if args.output:
                shutil.rmtree(args.output, ignore_errors=True)
            p = get_pockets(args.in_file)
            if args.show:
                show_pockets(p)
            if args.show_druggable:
                show_pockets_sorted(p)
            if args.open_pymol:
                path = str(args.in_file).rstrip('.pdb')
                prev_dir = os.getcwd()
                os.chdir(path + '_out/')
                path_to_pml = path + '_out/' + path[-4:] + '.pml'
                open_pymol(path_to_pml)
                os.chdir(prev_dir)
        else:
            print('No such file exists')
            sys.exit()
    if args.get_pdb:
        if len(args.get_pdb) != 4:
            print("Not a valid PDB code")
            sys.exit()
        else:
            if not args.output:
                print('No output name given')
                sys.exit()
            path = drug.get_pdb([str(args.get_pdb).upper()], args.output + '/')
            p = get_pockets(path)
            if args.show:
                show_pockets(p)
            if args.show_druggable:
                show_pockets_sorted(p)
            if args.open_pymol:
                path_to_pml = args.get_pdb + '.pml'
                prev_dir = os.getcwd()
                out_path = str(path).rstrip('.pdb')
                os.chdir(out_path + args.get_pdb + '_out/')
                open_pymol(path_to_pml)
                os.chdir(prev_dir)
