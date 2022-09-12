import dadi
import pandas
import my_models
import pylab

"""
This script is a collection of functions needed in order to plot dadi results
"""

def get_pop_pair_projs(pops):
    """
    this function will select the correct projections for each population pair (see best_projections.py and dadi.md)
    """
    if pops == "lpa-wel":
        proj_pop1 = 34
        proj_pop2 = 36

    if pops == "lpa-eel":
        proj_pop1 = 34
        proj_pop2 = 34
    
    if pops == "lpa-sel":
        proj_pop1 = 34
        proj_pop2 = 22
    
    proj = [proj_pop1, proj_pop2]
    
    return proj

def get_data_fs(vcf_file, pop_file, pops, proj):
    """
    this function takes a VCF file, a dadi population file and the name of the population pair (pop1-pop2)
    to build an unpolarized frequency spectrum (FS) object of the data for dadi.
    """
    # build data dictionary
    dd = dadi.Misc.make_data_dict_vcf(vcf_file, pop_file)
    
    # population id's from the pop pair    
    pop_ids = pops.split("-")

    # build the FS object
    fs_data = dadi.Spectrum.from_data_dict(dd, pop_ids = pop_ids, projections = proj, polarized = False)
    
    return fs_data
    
def get_models_params_n_best_opti_results(opti_table, chosen_model, n):
    """
    this function will take an ordered optimization table (see order_optimization_results.R and dadi.md) and
    extract the Nth scoring (0 based index) parameters for the chosen model
    """
    # read ordered optimization table
    df = pandas.read_csv(opti_table, sep='\t')
    
    # extract the results from the chosen model
    df_model = df[df.Model == chosen_model].reset_index()
    
    # create a dictionary nth scoring results 
    row_dict = df_model.iloc[n].to_dict()
    
    # select appropriate parameters for chosen model
    if chosen_model == 'model_1_a' or chosen_model == 'model_1_b':
        params = [row_dict["Tsplit"], row_dict["Tbot1"], row_dict["iber_a"], row_dict["iber_pr"],\
                  row_dict["eura_a"], row_dict["eura_pr"], row_dict["m"], row_dict["m_12"], row_dict["m_21"]]
    
    if chosen_model == 'model_2_a' or chosen_model == 'model_2_b' or chosen_model == 'model_2_c':
        params = [row_dict["Tsplit"], row_dict["Tbot2"], row_dict["Tbot1"], row_dict["iber_a"],\
                  row_dict["iber_pr_a"], row_dict["iber_pr"], row_dict["eura_a"], row_dict["eura_pr"],\
                  row_dict["m"], row_dict["ma_12"], row_dict["ma_21"], row_dict["m_12"], row_dict["m_21"]]
    
    return params


def get_model_obj(chosen_model, params, proj):
    """
    this function will create a dadi model object for the chosen model, using a list of appropriate parameters
    and selected projections
    """
    # define pts here so I can change it more easily
    pts = 60

    if chosen_model == "model_1_a":
    	model = my_models.model_1_a(params=params, ns=proj, pts=pts)
    
    elif chosen_model == "model_1_b":
    	model = my_models.model_1_b(params=params, ns=proj, pts=pts)
    
    elif chosen_model == "model_2_a":
    	model = my_models.model_2_a(params=params, ns=proj, pts=pts)
    
    elif chosen_model == "model_2_b":
    	model = my_models.model_2_b(params=params, ns=proj, pts=pts)
    
    elif chosen_model == "model_2_c":
    	model = my_models.model_2_c(params=params, ns=proj, pts=pts)
    
    elif chosen_model == "model_2_d":
    	model = my_models.model_2_d(params=params, ns=proj, pts=pts)
    
    model_folded = model.fold()
    
    return model_folded

def plot_save_fs_vs_model_2D(plot_file_name, model_folded, fs_data):
    """
    this function will generate a 2D SFS multinom plot with the selected file name with data from
    two dadi FS objects, the folded model (from get_model_obj) and the data (from get_data_fs)
    """
    # get population ids
    pop_ids = fs_data.pop_ids

    # build and save plot
    fig = pylab.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model_folded, fs_data, pop_ids = pop_ids, residual = 'Anscombe')
    fig.savefig(plot_file_name)