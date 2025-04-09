#!/usr/bin/python3.13

import os
import sys
import argparse
import traceback
from datetime import datetime
from classes.Matcher import URAProteomeMatcher
from classes.Log import Logger
from classes.Config import Config
from classes.UI import *

def arg_parser():
      """ Defines named arguments, provides the --help/-h flag, and allows the option to use required arguments """
      parser = argparse.ArgumentParser(description = 'arguments')
      parser.add_argument('-pd', '--proteome-data-file', help = 'REQUIRED: Local file path to proteome data', nargs = '?', type = str, required = True)
      parser.add_argument('-ur', '--upstream-regulator-file', help = 'REQUIRED: Local file path to upstream regulator data', nargs = '?', type = str, required = True)
      parser.add_argument('-ug', '--upstream-regulator-group-file', help = 'OPTIONAL: Local file path to json gene group file', nargs = '?', type = str, required = False)
      parser.add_argument('-ipn', '--ignore-protein-name-matches', help = 'OPTIONAL: Only match Genes, ignore Protein_Names matches', action = 'store_true')
      parser.add_argument('-ll', '--log-level', help = 'OPTIONAL: Log level, defaults to INFO', nargs = '?', type = str, required = False, default='INFO')
      parser.add_argument('-lf', '--log-to-file', help = 'OPTIONAL: Log to file, defaults to False', action = 'store_true')
      parser.add_argument('-of', '--out-file', help = 'OPTIONAL: Local file to output results to', nargs = '?', type = str, required = False)
      return parser.parse_args()

def main(args):
    """ Main function, parses proteome data file and adds upstream regulator matches/groups """
    # instantiate the logger and matcher
    current_time = datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
    logger = Logger.get_logger(args.log_level, args.log_to_file, current_time)
    matcher = URAProteomeMatcher(logger)

    # load input files
    print('\n--- load files ---\n')
    proteome_data = matcher.load_proteome_data(args.proteome_data_file)
    upstream_regulators = matcher.load_upstream_regulators(args.upstream_regulator_file)
    group_data = matcher.load_group_data(args.upstream_regulator_group_file)

    # drop match columns from source dataframe if they are present
    print('\n--- drop old match column results (if present) ---\n')
    match_columns = ['UR_EXACT_MATCH', 'UR_GROUP_MATCH', 'UR_BEST_MATCH']
    matcher.drop_columns(proteome_data, match_columns)

    # validate proteome data, confirm expected columns are present
    print('\n--- validate proteome data has expected columns ---\n')
    if args.ignore_protein_name_matches:
        expected_columns = ['Genes']
    else:
        expected_columns = ['Protein_Names', 'Genes']
    matcher.validate_proteome_data(proteome_data, expected_columns)

    # find URs for exact, group, and best match columns
    print('\n--- look for matches, generate match data ---\n')
    match_data = []
    ## error handle this in case protein names isnt present
    for index, row in proteome_data.iterrows():
        genes = row.get('Genes', '').split(';')
        protein_names = row.get('Protein_Names', '').split(';')

        ur_exact_matches = matcher.get_ur_exact_matches(index, upstream_regulators, genes, protein_names, args.ignore_protein_name_matches)
        ur_group_matches = matcher.get_ur_group_matches(index,upstream_regulators, group_data, genes)
        ur_best_matches = matcher.get_ur_best_matches(index, ur_exact_matches, ur_group_matches)

        match_data.append([';'.join(ur_exact_matches), ';'.join(ur_group_matches), ';'.join(ur_best_matches)])

    # write match data to the proteome data
    print('\n--- update proteome data with matches---\n')
    updated_proteome_data = matcher.update_proteome_data(proteome_data, match_columns, match_data)    
    

    # find the predicted activation state for your UR EXACT_MATCH and BEST_MATCH columns
    prediction_data = []
    for index, row in updated_proteome_data.iterrows():
        exact_match_urs = row.get('UR_EXACT_MATCH', '').split(';')
        ur_em_predictions = matcher.get_em_predictions(index, upstream_regulators, exact_match_urs)
        
        best_match_urs = row.get('UR_BEST_MATCH', '').split(';')
        ur_bm_predictions = matcher.get_bm_predictions(index, upstream_regulators, best_match_urs)
        
        prediction_data.append([';'.join(ur_em_predictions), ';'.join(ur_bm_predictions)])
    
    predic_columns = ['EXACT_MATCH_PREDIC', 'BEST_MATCH_PREDIC']
    proteome_data_with_predicts = matcher.update_proteome_data(updated_proteome_data, predic_columns, prediction_data)  

        
    # write the updated proteome data to a file
    print('\n--- write updated proteome data to file ---\n')
    out_file = os.path.abspath(args.out_file or f'URAProteomeMatcher_OUT_{current_time}.csv')
    matcher.write_dataframe_file(proteome_data_with_predicts, out_file) 
    return out_file

def ui_main():
    root = Tk()
    root.title('URAProteomeMatcher')
    frame = Frame(root)
    frame.pack(fill=X)

    pd_entry = create_entry(frame, '')
    pd_button = create_button(frame, 'Select Proteome CSV File', get_file_path, pd_entry)

    ur_entry = create_entry(frame, '')
    ur_button = create_button(frame, 'Select Upstream Regulator Data', get_file_path, ur_entry)

    urg_entry = create_entry(frame, '')
    urg_button = create_button(frame, 'Select Upstream Regulator Group Reference', get_file_path, urg_entry)

    of_entry = create_entry(frame, '')
    of_label = create_label(frame, 'Output CSV File Name (Optional)')

    ll = create_optionmenu(frame, ['info', 'debug'])
    ipnm = create_checkbox(frame, 'Ignore Protein Name Matches')

    def on_submit():
        ui_config = {
            'proteome_data_file': pd_entry.get(),
            'upstream_regulator_file': ur_entry.get(),
            'upstream_regulator_group_file': urg_entry.get(),
            'out_file': of_entry.get() or None,
            'log_level': ll.get() or 'info',
            'ignore_protein_name_matches': ipnm.get(),
            'log_to_file': False # Turned off due to error on Windows
        }
        config = Config()
        for k, v in ui_config.items():
            setattr(config, k, v)
        try:
            out_file = main(config)
            open_popup(root, 'Done!', f'Output File: {out_file}')
            open(out_file)
        except Exception as e:
            open_popup(root, 'ERROR', f'ERROR: {e}:\n{traceback.format_exc()}')


    submit_button = Button(frame, text='Submit', command=on_submit)
    submit_button.pack(fill=X)

    frame.mainloop()

if __name__ == "__main__":
    try:
        if len(sys.argv) == 1:
            # run UI if no args are provided
            ui_main()
        else:
            # parse arguments, generate a logger at specified level, instantiate matcher class, and call main
            args = arg_parser()
            main(args)
    except Exception as e:
        # wrap main in a try/except to raise an Exception as soon as an error occurs
        raise Exception('Error', e) 

