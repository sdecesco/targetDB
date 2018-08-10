#!/usr/bin/env python

import re
from tkinter.filedialog import askopenfilename,askdirectory
from tkinter.simpledialog import askstring
from pathlib import Path


def get_config_from_user(config, todo=None, new=False):
	print('========================================================')
	print('=================== CONFIGURATION ======================')
	print('========================================================')
	if new:
		config['database_path'] = {}
		config['output_path'] = {}
		config['pubmed_email'] = {}
		config['executables'] = {}
		config['database_path']['targetdb'] = ''
		config['database_path']['chembl'] = ''
		config['pubmed_email']['email'] = ''
		config['output_path']['db_files'] = ''
		config['output_path']['single'] = ''
		config['output_path']['list'] = ''
		config['executables']['blast'] = ''
		config['executables']['blastdb_path'] = ''
		config['executables']['fpocket'] = ''

	if 'chembl' in todo:
		path_exist = False
		while not path_exist:
			print('========================================================')
			print('================ CHEMBL SQLITE FILE ====================')
			print('========================================================\n')
			chembldb_path = Path(askopenfilename(title='Select ChEMBL sqlite database',
			                           filetypes=[("sqliteDB", "*.db")]))
			if chembldb_path.is_file():
				path_exist = True
			else:
				print('[ERROR]: The file you have entered does not exists \n\n')
		config['database_path']['chembl'] = str(chembldb_path)

	if 'targetdb' in todo:
		path_exist = False
		while not path_exist:
			print('========================================================')
			print('=============== TargetDB SQLITE FILE ===================')
			print('========================================================\n')
			targetDB_path = Path(askopenfilename(title='Select TargetDB sqlite database',
						                           filetypes=[("sqliteDB", "*.db")]))
			if targetDB_path.is_file():
				path_exist = True
			else:
				print('[ERROR]: The file you have entered does not exists \n\n')
		config['database_path']['targetdb'] = str(targetDB_path)

	if 'blast' in todo:
		path_exist = False
		while not path_exist:
			print('========================================================')
			print('================ BLASTP EXECUTABLE  =====================')
			print('========================================================\n')
			blast = Path(askopenfilename(title='Select BlastP executable location'))
			blast_db = Path(askdirectory(title='Select directory where blast database are located'))
			if blast.is_file() and blast_db.is_dir():
				path_exist = True
			else:
				print('[ERROR]: The file you have entered does not exists \n\n')
		config['executables']['blast'] = str(blast)
		config['executables']['blastdb_path'] = str(blast_db)

	if 'fpocket' in todo:
		path_exist = False
		while not path_exist:
			print('========================================================')
			print('================ FPOCKET EXECUTABLE  ===================')
			print('========================================================\n')
			fpocket = Path(askopenfilename(title='Select fPocket executable location'))
			if fpocket.is_file():
				path_exist = True
			else:
				print('[ERROR]: The excutable you have entered does not exists \n\n')
		config['executables']['fpocket'] = str(fpocket)

	if 'single' in todo:
		path_exist = False
		while not path_exist:
			print('========================================================')
			print('=========== SINGLE TARGET OUTPUT FOLDER ================')
			print('========================================================\n')
			single_output = Path(askdirectory(title='Select directory to save single output files'))
			if single_output.is_dir():
				path_exist = True
			else:
				print('[ERROR]: The folder you have entered does not exists \n\n')
		config['output_path']['single'] = str(single_output)

	if 'list' in todo:
		path_exist = False
		while not path_exist:
			print('========================================================')
			print('============= LIST TARGET OUTPUT FOLDER ================')
			print('========================================================\n')
			list_output = Path(askdirectory(title='Select directory to save list output files'))
			if list_output.is_dir():
				path_exist = True
			else:
				print('[ERROR]: The folder you have entered does not exists \n\n')
		config['output_path']['list'] = str(list_output)

	if 'db_files' in todo:
		path_exist = False
		while not path_exist:
			print('========================================================')
			print('========= DATABASE WORKING FILES OUTPUT  ===============')
			print('=====/!\ WARNING A LOT OF SPACE CAN BE REQUIRED /!\=====')
			print('========================================================\n')
			db_output = Path(askdirectory(title='Select directory to save working database files (large files)'))
			if db_output.is_dir():
				path_exist = True
			else:
				print('[ERROR]: The folder you have entered does not exists \n\n')
		config['output_path']['db_files'] = str(db_output)

	if 'email' in todo:
		email_correct = False
		while not email_correct:
			print('========================================================')
			print('============= EMAIL FOR PUBMED SEARCH ==================')
			print('========================================================\n')
			pubmed_email = askstring("Enter your email",'Enter your email address (used for pubmed searches - pubmed api requires an email for batch requests)')
			if is_email(pubmed_email):
				email_correct = True
			else:
				print('[ERROR]: The email address you have entered is invalid (correct format: youraddress@domain.com) \n')

		config['pubmed_email']['email'] = pubmed_email

	return config


def is_email(email):
	match = re.match('^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,4})$', email)
	if match is None:
		return False
	return True