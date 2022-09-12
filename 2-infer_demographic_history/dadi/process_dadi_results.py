import sys
import pandas
import argparse

# read parameters from optimization table
def get_params_dict_n_best_opti(opti_table, chosen_model, n):
    # read ordered optimization table
    df = pandas.read_csv(opti_table, sep='\t')
    # extract the results from the chosen model
    df_model = df[df.Model == chosen_model].reset_index()
    # create a dictionary nth scoring results 
    params_dict = df_model.iloc[n].to_dict()
    return params_dict

# calculate Nref
def get_Nref(theta, mu, L):
    # theta = 4 * Nref * mu * L
    Nref = theta / (4 * mu * L)
    return Nref

# calculate time
def get_time(dadi_time, Nref):
    # time is given in units of 2 Nref
    real_time = dadi_time * 2 * Nref
    return real_time

# calculate population size
def get_nu(dadi_nu, Nref):
    real_nu = dadi_nu * Nref
    return real_nu

# calculate migration
def get_mig(dadi_mig, Nref):
    # migration is given in units of Mij = 2 Nref mij
    # proportion of individuals that are migrants from another population in a given generation
    real_mig = dadi_mig / (2 * Nref)
    return real_mig

# main function
def print_processed_table(opti_table):
    # mutation rate
    mu = 6e-9
    # table header
    print("nth_best", "model", "pop_pair", "log-likelihood", "Nref", "Tsplit", "Tbot1", "Tbot2", "N_iber_c", "N_iber_i", "N_iber_a", "N_eura_c", "N_eura_a", "ancient_mig", "inter_mig_12", "inter_mig_21", "recent_mig_12", "recent_mig_21", sep = '\t')
    # for each model
    for model in ['model_1_a', 'model_1_b', 'model_2_a', 'model_2_b', 'model_2_c']:
        # extract results of first 4 optimizations
        for n in range(4):
            # create dictionary
            params_dict = get_params_dict_n_best_opti(opti_table, model, n)
            # segment length
            if params_dict['pop_pair'] == "lpa-wel":
                L = 612311182
            elif params_dict['pop_pair'] == "lpa-eel":
                L = 612178519
            elif params_dict['pop_pair'] == "lpa-sel":
                L = 612085036
            # calculate Nref
            Nref = round(get_Nref(params_dict['theta'], mu, L), 2)
            # calculate split time (tsplit + tbot1 + tbot2) * 5 (gen time)
            Tsplit = round((get_time(params_dict['Tsplit'], Nref) + get_time(params_dict['Tbot2'], Nref) + get_time(params_dict['Tbot1'], Nref)) * 5, 2)
            # calculate first bottleneck time (tbot1 + tbot2) * 5 (gen time)
            Tbot1 = round((get_time(params_dict['Tbot2'], Nref) + get_time(params_dict['Tbot1'], Nref)) * 5, 2)
            # calculate second bottleneck time (tbot2) * 5 (gen time)
            Tbot2 = round((get_time(params_dict['Tbot2'], Nref)) * 5, 2)
            # calculate Ne of contemporary (after bot1 and bot2) iberian
            N_iber_c = round(get_time(params_dict['iber_a'], Nref) * params_dict['iber_pr_a'] * params_dict['iber_pr'], 2)
            # calculate Ne of intermediate (after bot1) iberian
            N_iber_i = round(get_time(params_dict['iber_a'], Nref) * params_dict['iber_pr_a'], 2)
            # calculate Ne of ancient (before bottlenecks) iberian
            N_iber_a = round(get_time(params_dict['iber_a'], Nref), 2)
            # calculate Ne of contemporary (after bottleneck) eurasian
            N_eura_c = round(get_time(params_dict['eura_a'], Nref) * params_dict['eura_pr'], 2)
            # calculate Ne of ancient (before bottleneck) eurasian
            N_eura_a = round(get_time(params_dict['eura_a'], Nref), 2)
            # calculate after split symmetrical migration rate
            ancient_mig = get_mig(params_dict['m'], Nref)
            # calculate intermediate (after bot1) migration rate into iberian (1) from eurasian (2)
            inter_mig_12 = get_mig(params_dict['ma_12'], Nref)
            # calculate intermediate (after bot1) migration rate into eurasian (2) from iberian (1)
            inter_mig_21 = get_mig(params_dict['ma_21'], Nref)
            # calculate recent (after bot2) migration rate into iberian (1) from eurasian (2)
            recent_mig_12 = get_mig(params_dict['m_12'], Nref)
            # calculate recent (after bot2) migration rate into eurasian (2) from iberian (1)
            recent_mig_21 = get_mig(params_dict['m_21'], Nref)
            
            # create a list of the calculated parameters
            real_params = [n, model, params_dict['pop_pair'], params_dict['log-likelihood'], Nref, Tsplit, Tbot1, Tbot2, N_iber_c, N_iber_i, N_iber_a, N_eura_c, N_eura_a, ancient_mig, inter_mig_12, inter_mig_21, recent_mig_12, recent_mig_21]
            
            # print parameters
            params_line = '\t'.join(map(str, real_params))
            print(params_line)

# run main function if executing script
if __name__ == '__main__':
    
    # parse arguments for script
    parser = argparse.ArgumentParser()
    
    parser.add_argument("opti_table", type=str,
                    help="file containing ordered optimization results")

    args = parser.parse_args()
    opti_table = args.opti_table

    print_processed_table(opti_table)