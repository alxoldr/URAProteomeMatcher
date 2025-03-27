#!/usr/bin/python3.13

import os
import json
import pandas as pd
import logging
#import collections

from tabulate import tabulate

class URAProteomeMatcher():
    """ Class used to identify upstream regulators for given protein or gene names """

    def __init__(self, logger=logging.getLogger(__name__)):
        self.logger = logger

    def load_proteome_data(self, proteome_data_file: str):
        """
        Creates a pandas dataframe called proteome_data from a provided CSV file, e.g.

        Accession,Protein_Names,Genes,Mean_Intensity
        Q9NRM6,I17RB_HUMAN,IL17RB,693300
        P30939;P21918,5HT1F_HUMAN;DRD5_HUMAN,HTR1F;DRD5,418100
        Q9UMX1,SUFU_HUMAN,SUFU,15700
        P00734,THRB_HUMAN,F2,122800
        P10828,THB_HUMAN,THRB,134700
        Q9Y6Y9,LY96_HUMAN,LY96,170300
        P01189,COLI_HUMAN,POMC,189900


        ╒════╤═══════════════╤════════════════════════╤════════════╤═════════════╕
        │    │   Accession   │          Protein_Names │   Genes    │   Intensity │
        ╞════╪═══════════════╪════════════════════════╪════════════╪═════════════╡
        │ 0  │    Q9NRM6     │            I17RB_HUMAN │   IL17RB   │      693300 │
        ├────┼───────────────┼────────────────────────┼────────────┼─────────────┤
        │ 1  │ P30939;P21918 │ 5HT1F_HUMAN;DRD5_HUMAN │ HTR1F;DRD5 │      418100 │
        ├────┼───────────────┼────────────────────────┼────────────┼─────────────┤
        │ 2  │    Q9UMX1     │             SUFU_HUMAN │    SUFU    │       15700 │
        ├────┼───────────────┼────────────────────────┼────────────┼─────────────┤
        │ 3  │    P00734     │             THRB_HUMAN │     F2     │      122800 │
        ├────┼───────────────┼────────────────────────┼────────────┼─────────────┤
        │ 4  │    P10828     │              THB_HUMAN │    THRB    │      134700 │
        ├────┼───────────────┼────────────────────────┼────────────┼─────────────┤
        │ 5  │ O00206;A8K1Y8 │             TLR4_HUMAN │    TLR4    │      170300 │
        ├────┼───────────────┼────────────────────────┼────────────┼─────────────┤
        │ 6  │    P01189     │             COLI_HUMAN │    POMC    │      189900 │
        ╘════╧═══════════════╧════════════════════════╧════════════╧═════════════╛

        """
        file_path = os.path.abspath(proteome_data_file)
        proteome_data = pd.read_csv(file_path).replace('\xa0', '', regex=True)
        self.logger.info(f'LOADED PROTEOME: {file_path}')
        self.logger.debug('\n' + tabulate(proteome_data, proteome_data.columns, tablefmt='fancy_grid') + '\n')
        return proteome_data

    def load_upstream_regulators(self, upstream_regulator_file: str):
        """
        Creates a list of upstream regulators to search from from a given file (containing URs seperated by newline), e.g.

        ITAE
        IL-17 receptor family
        CD5
        THRB
        TLR4 complex (LPS-binding)
        GRPR
        DRD5
        FOS
        TLR4
        AP-2 complex


        """
        file_path = os.path.abspath(upstream_regulator_file)

        # # do we want to keep this as a csv file? seems like it doesn't have to be
        upstream_regulators = []
        with open(file_path) as file:
            upstream_regulators = list(set(filter(None, file.read().split("\n"))))
            self.logger.info(f'LOADED UPSTREAM REGULATORS: {file_path}')
            self.logger.debug('\n' + tabulate({ None : upstream_regulators }, ['UPSTREAM_REGULATORS'], tablefmt='fancy_grid') + '\n')
            file.close()

        return upstream_regulators

    def load_group_data(self, group_file: str):
        """
        Creates a dictionary of gene groups based on a provided JSON file, e.g.

        {
            "TLR4 complex (LPS-binding)": ["TLR4","LY96"],
            "IL-17R family": ["IL17RA","IL17RB","IL17RC","IL17RD"],
            "AP-2 complex": ["AP2A1","AP2A2"]
        }

        """
        if not group_file:
            self.logger.warning('NO GROUP FILE PROVIDED')
            return

        group_data = {}
        file_path = os.path.abspath(group_file)
        with open(file_path, 'r') as file:
            group_data = json.load(file)
            self.logger.info(f'LOADED GROUP DATA: {file_path}')
            self.logger.debug('\n' + tabulate(group_data, group_data.keys(), tablefmt='fancy_grid'))
            file.close()

        return group_data

    def validate_proteome_data(self, proteome_data: pd.DataFrame, expected_columns: list):
        """ Confirms that expected columns are present in proteome_data """
        expected_columns = set(expected_columns)
        if expected_columns.issubset(set(proteome_data.columns)):
            self.logger.info(f'EXPECTED COLUMNS FOUND: {expected_columns}')
        else:
            raise ValueError(f'MISSING COLUMNS IN proteome DATA: {expected_columns}')

    def drop_columns(self, dataframe: pd.DataFrame, columns: list):
        """ Drop passed columns if present """
        for column in columns:
            if column in dataframe.columns:
                self.logger.info(f'DROPPED COLUMN: {column}')
                dataframe.drop(column, axis=1)
        self.logger.info('DONE')

    def get_ur_exact_matches(self, index: int, upstream_regulators: list, genes: list, protein_names: list, ignore_protein_name_matches: bool):
        """
        Looks for upstream regulator matches in a list of genes, and optionally in protein names
        For gene matches, the upstream regulator must be an exact match, e.g.

            DRD5 (upstream regulator) in HTR1F;DRD5 (genes) = match

        For protein name matches, the the upstream regulator with a suffix of _MOUSE or _HUMAN must be present in protein names, e.g.

            THRB (upstream regulator) in THRB_HUMAN (protein names) = match

        Matches will get written in the UR_EXACT_MATCH column as a ; separated list
        """
        ur_exact_matches = set()
        for regulator in upstream_regulators:
            if regulator in genes:
                self.logger.info(f'ROW {index} UR_EXACT_MATCH {regulator} FOUND IN Genes: {";".join(genes)}')
                ur_exact_matches.add(regulator)
            elif f'{regulator}_HUMAN' in protein_names or f'{regulator}_MOUSE' in protein_names:
                if not ignore_protein_name_matches:
                    self.logger.info(f'ROW {index} UR_EXACT_MATCH {regulator} FOUND IN Protein_Names: {";".join(protein_names)}')
                    ur_exact_matches.add(regulator)
            else:
                pass
        return ur_exact_matches

    def get_ur_group_matches(self, index: int, group_data: dict, genes: list):
        """
        Looks for group matches based on a provided JSON file containing groups and the gene names within them, e.g.
        {
            "TLR4 complex (LPS-binding)": ["TLR4","LY96"],
            "IL-17R family": ["IL17RA","IL17RB","IL17RC","IL17RD"],
            "AP-2 complex": ["AP2A1","AP2A2"]
        }

        For each gene in a row, this will search for a match in the dict values
        If a match is found, the group key followed by the matched gene will get recorded in the UR_GROUP_MATCH column, e.g.

            IL-17 receptor family [IL17RB]
        
        """
        group_matches = {}
        if group_data:
            for group, group_genes in group_data.items():
                for gene in genes:
                    if gene in group_genes:
                        if group not in group_matches.keys():
                            group_matches[group] = []
                        self.logger.info(f'ROW {index} UR_GROUP_MATCH {group} [{gene}] FOUND')
                        group_matches[group].append(gene)

        ur_group_matches = set([f"{k} [{','.join(v)}]" for k, v in group_matches.items()])
        return ur_group_matches

    def get_ur_best_matches(self, index: int, ur_exact_matches: set, ur_group_matches: set):
        """
        Writes a match value to the UR_BEST_MATCH column in the following priority:

        1) UR_EXACT_MATCH
        2) UR_GROUP_MATCH
        3) NULL
        """
        ur_best_matches = set()
        if ur_exact_matches:
            self.logger.info(f'ROW {index} UR_BEST_MATCH EXACT: {";".join(ur_exact_matches)}')
            ur_best_matches.update(ur_exact_matches)
        elif ur_group_matches:
            self.logger.info(f'ROW {index} UR_BEST_MATCH GROUP: {";".join(ur_group_matches)}')
            ur_best_matches.update(ur_group_matches)
        else:
            pass
        return ur_best_matches

    def update_proteome_data(self, proteome_data: pd.DataFrame, match_columns: list, match_data: list):
        """ Merges the proteome_data dataframe with the match data """
        match_df = pd.DataFrame(columns=match_columns, data=match_data)
        updated_proteome_data = pd.merge(proteome_data, match_df, left_index=True, right_index=True, how='outer').fillna('')
        self.logger.info('UPDATED PROTEOME WITH URS')
        self.logger.debug('\n' + tabulate(updated_proteome_data, updated_proteome_data.columns, tablefmt='fancy_grid') + '\n')
        return updated_proteome_data

    def write_dataframe_file(self, dataframe: pd.DataFrame, out_file: str):
        OUT_FILE = os.path.abspath(out_file)
        dataframe.to_csv(OUT_FILE, index=False)
        self.logger.info(f'WROTE FILE TO: {OUT_FILE}')
