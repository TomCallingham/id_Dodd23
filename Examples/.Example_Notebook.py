# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import astropy
import numpy as np
from astropy.table import QTable
from id_Dodd23.coordinates import add_units

# %%
from id_Dodd23.load_data import data_folder, test_data_file

# %%
astro_tab = QTable.read(test_data_file)
astro_tab = add_units(astro_tab)

# %%
astro_tab

# %%
astro_tab["distance"]

# %%
observables = ["ra","dec","distance","pmra","pmdec","radial_velocity"]
obs_data = {key:astro_tab[key]for key in observables}

# %%
from id_Dodd23.coordinates import coord_transform_icrs_Galacto

# %%
xyz_observable = coord_transform_icrs_Galacto(obs_data)

# %%
print(xyz_observable[:,0])

# %%
print(np.asarray(astro_tab["X"]))

# %%
from id_Dodd23.dynamics import dynamics_calc_H99

# %%
vec_ELzLp  = dynamics_calc_H99(xyz_observable)

# %%
print(vec_ELzLp[:5,:])

# %%
print(astro_tab[:5][["E","Lz","Lp"]])

# %%
from id_Dodd23.id_stars import  add_maha_members_to_groups

# %%
from id_Dodd23.load_data import named_Groups, group_covar,group_mean 

# %%
print(group_mean)

# %%
labels = add_maha_members_to_groups(vec_ELzLp,group_mean,group_covar,max_dis=2.13)

# %%
print(labels)

# %%
original_labels = np.asarray(astro_tab["derived_labels_group"]).astype(int)

# %%
print(original_labels)

# %%
(original_labels==labels).all()

# %%
original_labels==labels

# %%
