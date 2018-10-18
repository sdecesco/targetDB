#!/usr/bin/env python

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import pkg_resources
from pathlib import Path


def generate_model():
    # Model generated with optimised parameters (see .ipynb file)
    rf_model = RandomForestClassifier(n_estimators=1000, max_depth=21, max_features=3, min_samples_leaf=2,
                                      min_samples_split=5)

    # recover the training data from the data file

    pck_path = Path(str(pkg_resources.resource_filename('targetDB.utils', ''))).parent
    ml_data = pck_path.joinpath('ml_data')
    ml_data = ml_data.joinpath('ml_training_data.zip')

    training_df = pd.read_json(ml_data, compression='zip')

    col_to_keep = ['OT_max_association_score', 'dis_AScore', 'gen_AScore', 'safe_GScore', 'dis_TScore', 'bio_PScore',
                   'bio_GScore', 'dis_UScore', 'alt_drug_score', '%_sequence_covered', 'mean_druggability_score',
                   '%_domain_covered', 'pdb_count_score', 'bio_MScore', 'chembl_potent', 'domain_drug_score',
                   'gen_GScore', 'bio_AScore', 'safe_EScore', 'druggability']
    training_df = training_df[col_to_keep]

    training_set, training_labels = training_df.drop('druggability', axis=1), training_df['druggability']

    return rf_model.fit(training_set, training_labels)


def predict(model, data):
    df = data.copy()
    df.index = df.Target_id
    df.drop(columns=['Target_id'], inplace=True)
    df.replace({True: 1, False: 0}, inplace=True)
    df = df.fillna(0)
    # Only columns to consider in the model (see .ipynb file for selection of the columns)

    col_to_keep = ['OT_max_association_score', 'dis_AScore', 'gen_AScore', 'safe_GScore', 'dis_TScore', 'bio_PScore',
                   'bio_GScore', 'dis_UScore', 'alt_drug_score', '%_sequence_covered', 'mean_druggability_score',
                   '%_domain_covered', 'pdb_count_score', 'bio_MScore', 'chembl_potent', 'domain_drug_score',
                   'gen_GScore', 'bio_AScore', 'safe_EScore']
    df = df[col_to_keep]
    return model.predict(df)


def predict_prob(model, data):
    df = data.copy()
    df.index = df.Target_id
    df.drop(columns=['Target_id'], inplace=True)
    df.replace({True: 1, False: 0}, inplace=True)
    df = df.fillna(0)
    # Only columns to consider in the model (see .ipynb file for selection of the columns)

    col_to_keep = ['OT_max_association_score', 'dis_AScore', 'gen_AScore', 'safe_GScore', 'dis_TScore', 'bio_PScore',
                   'bio_GScore', 'dis_UScore', 'alt_drug_score', '%_sequence_covered', 'mean_druggability_score',
                   '%_domain_covered', 'pdb_count_score', 'bio_MScore', 'chembl_potent', 'domain_drug_score',
                   'gen_GScore', 'bio_AScore', 'safe_EScore']
    df = df[col_to_keep]
    return model.predict_proba(df)


def in_training_set(data):

    pck_path = Path(str(pkg_resources.resource_filename('targetDB.utils', ''))).parent
    ml_data = pck_path.joinpath('ml_data')
    ml_data = ml_data.joinpath('ml_training_data.zip')

    training_df = pd.read_json(ml_data, compression='zip')

    df = data.copy()
    df.index = df.Target_id

    df['Is_in_training_set'] = 'No'
    df.loc[df.index.isin(training_df.index),['Is_in_training_set']] = 'Yes'
    return df['Is_in_training_set'].values
