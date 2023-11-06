import numpy as np
from modules.import_elife import DataImporter, Import_Theune
from colors import Colors
from import_model import import_model, split_rubisco
from zavrel_FIT import get_kl, pd_fitting, pd_fitting_w_resp
from zavrel_FIT import plus_blue, fit_haldane  # , proteincontent
from zavrel_FIT import photodamage_helper, dark_resp
from add_glycogen import update_glycogen, update_prot_glyc
from fit_glycogen_funcs import fit as fit_glyc
import multiprocessing as mp
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import optimize as op