# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: py13
#     language: python
#     name: py13
# ---

# %% [markdown]
# # Example Notebook
# How to:
# - Calculate Galactocentric Position Velocity from Observables
# - Calculate Energy, Lz, and Lperp in the potential of Dodd23
# - Identify the grouping of the star from its position in IoM space

# %% [markdown]
# ## Load Test Data
# Test data is a small subset of the GaiaDR3 Halo dataset of Dodd23.

# %%
import id_Dodd23

# %%
import numpy as np
import matplotlib.pyplot as plt

import astropy
from astropy.table import QTable
from id_Dodd23.coordinates import add_units

# %%
from id_Dodd23.load_data import test_data_file

# %%
astro_tab = QTable.read(test_data_file)
astro_tab = add_units(astro_tab)

# %%
astro_tab

# %% [markdown]
# ## Coordinantes
# Observables to Galactocentric XYZ

# %%
from id_Dodd23.coordinates import coord_transform_icrs_Galacto

# %%
observables = ["ra","dec","distance","pmra","pmdec","radial_velocity"]
obs_data = {key:astro_tab[key]for key in observables}

# %%
xyz = coord_transform_icrs_Galacto(obs_data)

# %% [markdown]
# ## Calculate Dynamics
# From galacto PosVel to \[En,Lz,Lperp\]

# %%
from id_Dodd23.dynamics import dynamics_calc_H99, calc_vtoomre, calc_vT_vP, vlsr
vec_ELzLp  = dynamics_calc_H99(xyz)

# %%
dist = obs_data["distance"].value

# %%
print(vec_ELzLp[:5,:])

# %%
print(astro_tab[:5][["E","Lz","Lp"]])

# %%
v_toomre = calc_vtoomre(xyz)
vT,vP = calc_vT_vP(xyz)

# %%
plt.figure()
plt.hist(v_toomre,bins=100)
plt.xlabel(r"$V_{\mathrm{toomre}}$, km/s")
plt.axvline(210,c="r")
plt.show()

# %% [markdown]
# ## Identify Groups
# Using the position in IoM space, identify stars that are close to known groups

# %%
from id_Dodd23.id_stars import  groups_from_dynamics
from id_Dodd23.load_data import named_Groups, group_covar,group_mean  , g_name_to_index
from id_Dodd23.plotting import draw_ellipse

# %%
groups = groups_from_dynamics(vec_ELzLp,v_toomre)

# %%
original_labels = np.asarray(astro_tab["derived_labels_group"]).astype(int)
original_groups = named_Groups[original_labels]
original_groups[original_labels==-1] = "Other"
original_groups[original_labels==-2] = "Disc"

# %%
(original_groups==groups).all()

# %%
plotting_groups = np.unique(groups)

# %%
for g in plotting_groups:
    pop = (groups==g).sum()
    print(g, pop)

# %% [markdown]
# ## Selection Plots

# %%
from id_Dodd23.plotting import g_colours,draw_ellipse,label_dic,lims

# %%
vT,vP = calc_vT_vP(xyz)

# %%
fig, axs = plt.subplots(ncols=3,figsize=(24,8))
for g in plotting_groups:
    if g in ["Disc","Other"]:
        continue
    i = g_name_to_index[g]
    g_filt =(groups==g)
    c = g_colours[g]
    draw_ellipse(group_mean[i,:2][::-1],group_covar[i,:2,:2][::-1,::-1],ax=axs[0],c=c)
    draw_ellipse(group_mean[i,1:],group_covar[i,1:,1:],ax=axs[1],c=c)
    pop = g_filt.sum()
    if pop==0:
        axs[1].scatter(None,None,c=c,label=compact_label.get(g,g),marker="+")
        continue
    print(g)
    scatter_args =  {"color":c, "ec":"k"}
    axs[0].scatter( vec_ELzLp[g_filt,1],vec_ELzLp[g_filt,0],label=f"{g} : {pop}",**scatter_args)
    axs[1].scatter( vec_ELzLp[g_filt,1],vec_ELzLp[g_filt,2],**scatter_args)
    axs[2].scatter( vT[g_filt], vP[g_filt],**scatter_args)



other_args = {"zorder":-2,"c":"k","s":80,"marker":".","ec":"none"}
g_filt = groups=="Other"
pop = g_filt.sum()
label = f"Not Grouped: {pop}"
axs[0].scatter( vec_ELzLp[g_filt,1],vec_ELzLp[g_filt,0],label=label,**other_args)
axs[1].scatter( vec_ELzLp[g_filt,1],vec_ELzLp[g_filt,2],**other_args)
axs[2].scatter( vT[g_filt],vP[g_filt],**other_args)

disc_args = {"zorder":-1,"c":"red","s":80,"marker":".","ec":"none"}
g_filt = groups=="Disc"
pop = g_filt.sum()
label =  f"Disc (Vtoomre<210) :{pop}"
axs[0].scatter( vec_ELzLp[g_filt,1],vec_ELzLp[g_filt,0],label=label,**disc_args)
axs[1].scatter( vec_ELzLp[g_filt,1],vec_ELzLp[g_filt,2],**disc_args)
axs[2].scatter( vT[g_filt],vP[g_filt],**disc_args)

axs[0].set_xlabel(label_dic["Lz"])
axs[0].set_ylabel(label_dic["En"])
axs[0].set_xlim(lims["Lz"])
axs[0].set_ylim(lims["En"])

axs[1].set_xlim(lims["Lz"])
axs[1].set_ylim(lims["Lp"])
axs[1].set_xlabel(label_dic["Lz"])
axs[1].set_ylabel(label_dic["Lp"])

v_toomre_circ = plt.Circle([vlsr.value,0],210,color="r",alpha=0.2)
axs[2].add_artist(v_toomre_circ)

axs[2].set_xlabel(label_dic["vT"])
axs[2].set_ylabel(label_dic["vP"])

handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.2))

plt.suptitle("Stars in Dodd23 Groups")
plt.show()
