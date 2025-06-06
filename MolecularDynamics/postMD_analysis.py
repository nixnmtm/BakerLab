import os
import os.path as path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
import plotly.graph_objs as go
import plotly.io as pio

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis import distances

import logging

# Figure rendering style for publication
def publication_style(fontsize):
   plt.rc('font', family='Arial', serif="Arial", weight='normal')
   #plt.rc('text', usetex=False)
   plt.rc('xtick', labelsize=fontsize)
   plt.rc('ytick', labelsize=fontsize)
   plt.rc('axes', labelsize=fontsize)
   plt.rc('axes', labelweight='normal')
   plt.rc('lines', markersize=1.2)
   plt.rc('lines', linewidth=0.8)
   plt.rc('legend', fontsize=fontsize)
    
def RMSD(ref, mobile, **kwargs):

    """ Get the RMSD between two structures 
    
    ref: reference struture
    mobile: the structure to check for RMSD and fit
    returns:
    rotation matrix and rmsd
    
    kwargs:
    
    selection: atom selection
    write: save new fitted mobile structure to pdb format
    
    """
    
    selection = kwargs.get('selection', None)
    writefile = kwargs.get('writefile', None)
    if selection is None:
        selection = "protein and name CA"
        
    ref = mda.Universe(ref)
    mobile = mda.Universe(mobile)
    mobile_sel = mobile.select_atoms(selection)
    ref_sel = ref.select_atoms(selection)
    print(f"# of atoms selected in mobile: {mobile_sel.atoms.n_atoms}")
    print(f"# of atoms selected in reference: {ref_sel.atoms.n_atoms}")
    
    ref0 =  ref_sel.positions - ref_sel.center_of_mass()
    mobile0 =  mobile_sel.positions - mobile_sel.center_of_mass()
    R, rmsd = align.rotation_matrix(mobile0, ref0)
    if writefile is not None:
        mobile.atoms.translate(-mobile_sel.center_of_mass())
        mobile.atoms.rotate(R)
        mobile.atoms.translate(ref_sel.center_of_mass())
        mobile.atoms.write(writefile)
    return R, rmsd
    

def plot_rmsd(xvglist, xticklabel=None, yticklabel=None, 
              xlimit=None, ylimit=None, fontsize=12, figsize=None, slide=1000, label=None,
             col=None, title="RMSD", savepath=None, outname=None, leg_loc="best", ticklabel_rotate=90):
    
    """
    Plot the RMSd of two systems
    :param slide, to convert the time points to ns scale
    :param: xvglist: list of xvg files to plot
    """
    #fig, ax = plt.subplots(1,figsize=(3.33,2.06))
    fig, ax = plt.subplots(1,figsize=figsize)
    for n,data in enumerate(xvglist):
        publication_style(fontsize=fontsize)
        #x1 = data[0][1:][::slide]/1000  # converting in us scale
        x1 = data[0][1:][::slide]
        y1 = data[1][1:][::slide]*10
        #print(len(x1), len(y1)) 
        
        if col is not None:
            ax.plot(x1, y1, label=label[n], color=col[n])
        else:
            ax.plot(x1, y1, label=label[n])
        x_formatter = FixedFormatter([str(i) for i in xticklabel])
        y_formatter = FixedFormatter([str(i) for i in yticklabel])
        x_locator = FixedLocator(xticklabel)
        y_locator = FixedLocator(yticklabel)
        ax.xaxis.set_major_formatter(x_formatter)
        ax.yaxis.set_major_formatter(y_formatter)
        ax.xaxis.set_major_locator(x_locator)
        ax.yaxis.set_major_locator(y_locator)
        ax.tick_params(axis='x', labelrotation=ticklabel_rotate)
        
        ax.set_xlabel(r"Time (ns)", fontsize=fontsize+2)
        ax.set_ylabel(r"RMSD ($\mathrm{\AA}$)", fontsize=fontsize+2)
        ax.set_title(title)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        if xlimit is not None:
            ax.set_xlim(xlimit[0], xlimit[1])
        else:
            ax.set_xlim(x1.min(), x1.max())
        if ylimit is not None:
            ax.set_ylim(ylimit[0], ylimit[1])
        else:
            ax.set_xlim(y1.min(), y1.max())
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #leg = ax.legend(loc=leg_loc, bbox_to_anchor=(0.5, 0.5), frameon=False)
        leg = ax.legend(loc=leg_loc, frameon=False)
        for line in leg.get_lines():
            line.set_linewidth(3.0)
        #ax.legend(frameon=False)
    if len(xvglist) != len(label):
        logging.WARNING("Number of labels is not equal to the number of line plots given")
    if figsize is None:
        fig.set_size_inches(width, height)
    fig.tight_layout()
    if outname is not None:
        fig.savefig(path.join(savepath,f"{outname}.pdf"))
    plt.show()

def read_rmsd_xvg(file, skiplines=18):
    _read = pd.read_csv(file, skiprows=skiplines, header=None, delim_whitespace=True, engine='python')
    return _read

def read_rmsd_xvg_list(filelist, skiplines=18):
    xvglist = []
    for file in filelist:
        xvglist.append(read_rmsd_xvg(file, skiplines=skiplines))
    return xvglist

def read_rmsf_xvg(filename):
    _read = pd.read_csv(filename, skiprows=17, header=None, delim_whitespace=True, engine='python')
    return _read

def plot_rmsf(filename,
              col='blue',
              outname=None,
              xlimit=None,
              ylimit=None,
              fontsize=12,
              label=None,
              ax=None,
              xtick_label_spacing=1,
              savepath='.'):
    """
    Plot the RMSF of one system from an XVG file with properly aligned major and minor ticks.

    Parameters:
        filename (str): Path to the XVG file.
        col (str): Line color.
        outname (str, optional): Base name for saved figure (without extension).
        xlimit (tuple, optional): Limits for the x-axis.
        ylimit (tuple, optional): Limits for the y-axis.
        fontsize (int): Base font size for publication style.
        label (str, optional): Legend label.
        ax (matplotlib.axes.Axes, optional): Axis to plot on.
        xtick_label_spacing (int): Interval for major tick labels (every Nth point).
        savepath (str): Directory to save the figure.
    Returns:
        matplotlib.axes.Axes: The axis with the plot.
    """
    # Apply publication style
    publication_style(fontsize=fontsize)

    # Read data
    data = read_rmsf_xvg(filename)
    # Convert x-values to numeric array (residue IDs)
    x1 = pd.to_numeric(data[0], errors='coerce').to_numpy()
    # RMSF values converted to Angstroms
    y1 = data[1].to_numpy() * 10

    # Create figure/axis if not supplied
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.figure

    # Plot RMSF
    ax.plot(x1, y1, label=label, linewidth=2, color=col)

    # Labels and formatting
    ax.set_xlabel("Residue", fontsize=fontsize+2)
    ax.set_ylabel(r"RMSF ($\text{\AA}$)", fontsize=fontsize+2)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Apply axis limits if provided
    if xlimit is not None:
        ax.set_xlim(*xlimit)
    if ylimit is not None:
        ax.set_ylim(*ylimit)

    # Legend
    if label:
        ax.legend(frameon=False)

    # Save figure
    if outname is not None:
        fig.tight_layout()
        fig.savefig(os.path.join(savepath, f"{outname}.pdf"))

    return ax

def get_filepath(basepath, protein, filename):
    xvgpath = path.join(basepath, protein)
    filepath = path.join(xvgpath, filename)
    return filepath

def compute_com_distance(universe, sel1_str, sel2_str, stride=1):
    sel1 = universe.select_atoms(sel1_str)
    sel2 = universe.select_atoms(sel2_str)

    times = []
    distances = []

    for ts in universe.trajectory[::stride]:
        com1 = sel1.center_of_mass()
        com2 = sel2.center_of_mass()
        dist = np.linalg.norm(com1 - com2)

        times.append(ts.time)
        distances.append(dist)

    return np.array(times), np.array(distances)

def angle_between_domains(universe, sel_ig20, sel_hinge, sel_ig21, stride=1):
    sel1 = universe.select_atoms(sel_ig20)
    sel2 = universe.select_atoms(sel_ig21)
    sel_hinge = universe.select_atoms(sel_hinge)

    times = []
    angles = []

    for ts in universe.trajectory[::stride]:
        com1 = sel1.center_of_mass()
        com2 = sel2.center_of_mass()
        com_hinge = sel_hinge.center_of_mass()

        vec1 = com1 - com_hinge  # vector from hinge to Ig20
        vec2 = com2 - com_hinge  # vector from hinge to Ig21

        # Normalize
        vec1 /= np.linalg.norm(vec1)
        vec2 /= np.linalg.norm(vec2)

        # Angle in degrees
        angle_rad = np.arccos(np.clip(np.dot(vec1, vec2), -1.0, 1.0))
        angle_deg = np.degrees(angle_rad)

        times.append(ts.time)
        angles.append(angle_deg)

    return np.array(times), np.array(angles)

def track_com_distance(universe, resid_target, core_resids, use_sidechain=True):
    """
    Track COM distance between a residue and a hydrophobic core.

    Parameters:
    - universe: MDAnalysis Universe
    - resid_target: int
    - core_resids: list of ints
    - use_sidechain: if True, excludes backbone atoms (N, C, CA, O)

    Returns:
    - times, distances
    """
    atom_filter = "not name N C CA O" if use_sidechain else "all"

    sel_target = universe.select_atoms(f"resid {resid_target} and {atom_filter}")
    sel_core = universe.select_atoms(" or ".join([f"resid {r} and {atom_filter}" for r in core_resids]))

    if sel_target.n_atoms == 0 or sel_core.n_atoms == 0:
        raise ValueError("Selection returned 0 atoms. Check residue IDs or atom names.")

    times = []
    distances = []

    for ts in universe.trajectory:
        com_target = sel_target.center_of_mass()
        com_core = sel_core.center_of_mass()
        dist = np.linalg.norm(com_target - com_core)
        times.append(universe.trajectory.time)
        distances.append(dist)

    return np.array(times), np.array(distances)

def plot_contact_map(
    u,
    sel1,
    sel2,
    label="",
    stride=10,
    distance_cutoff=5.0,
    contact_threshold=0.0,
    highlight_threshold=None,
    ax=None,
    cmap="coolwarm",
    vmin=0,
    vmax=None,
    annotate=True,
    tick_label_spacing=5
):
    """
    Plot contact map between two residue selections in an MDAnalysis Universe.

    Parameters:
    - u: MDAnalysis.Universe
    - sel1, sel2: selection strings (e.g., "resid 2140:2149 and backbone")
    - label: subplot title label
    - stride: frame stride for analysis
    - distance_cutoff: atom distance (Å) to define a contact
    - contact_threshold: minimum frequency to annotate a contact
    - highlight_threshold: frequency to draw highlight (e.g. 0.5 highlights persistent contact >= 0.5 )
    - ax: matplotlib axis object to plot on
    - cmap: color map
    - vmax: maximum for heatmap normalization
    """

    group1 = u.select_atoms(sel1).residues
    group2 = u.select_atoms(sel2).residues

    n1 = len(group1)
    n2 = len(group2)
    contact_matrix = np.zeros((n1, n2))
    print(contact_matrix.shape)
    for ts in u.trajectory[::stride]:
        for i, res_i in enumerate(group1):
            for j, res_j in enumerate(group2):
                if res_i.atoms.n_atoms > 0 and res_j.atoms.n_atoms > 0:
                    d = distances.distance_array(res_i.atoms.positions, res_j.atoms.positions)
                    if np.any(d < distance_cutoff):
                        contact_matrix[i, j] += 1

    contact_matrix /= len(u.trajectory[::stride])

    # Annotation labels
    annot_labels = np.where(
        contact_matrix > contact_threshold,
        np.round(contact_matrix, 2).astype(str),
        ""
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    if annotate == True:
        sns.heatmap(
            contact_matrix,
            annot=annot_labels,
            fmt="",
            cmap=cmap,
            xticklabels=group2.resids,
            yticklabels=group1.resids,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={"label": "Contact Frequency"},
            ax=ax,
            linewidths=0.3,
            linecolor="gray")
    else:
        sns.heatmap(
            contact_matrix,
            cmap=cmap,
            xticklabels=group2.resids,
            yticklabels=group1.resids,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={"label": "Contact Frequency"},
            ax=ax,
            linewidths=0.3,
            linecolor="gray")

    # Highlight strong contacts
    if highlight_threshold:
        for (i, j), val in np.ndenumerate(contact_matrix):
            if val >= highlight_threshold:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='yellow', lw=2))


    # Reduce label density: show every 5th residue
    xticks = np.arange(0, len(group2.resids), tick_label_spacing)
    yticks = np.arange(0, len(group1.resids), tick_label_spacing)
    ax.set_xticks(xticks)
    ax.set_xticklabels([group2.resids[i] for i in xticks], rotation=90, ha='right', fontsize=8)
    ax.set_yticks(yticks)
    ax.set_yticklabels([group1.resids[i] for i in yticks], fontsize=8)
    
    ax.set_title(label)
    ax.set_xlabel("Residues: " + sel2)
    ax.set_ylabel("Residues: " + sel1)
    ax.invert_yaxis()
    return contact_matrix

def plot_contact_difference_heatmap(contact_WT, contact_MUT, resids_x, resids_y, 
                                    annotate=False,    
                                    vmin=-1,
                                    vmax=1,
                                    output_name="contact_diff_map.png", tick_label_spacing=5):
    """
    Plot the difference in contact frequency between MUT and WT.
    
    Args:
        contact_WT (np.ndarray): WT contact matrix (2D array)
        contact_MUT (np.ndarray): MUT contact matrix (same shape as WT)
        resids_y (list or array): Y-axis (ex: Ig20 β-strand) residue numbers
        resids_x (list or array): X-axis (ex: Ig21 CD-face) residue numbers
        output_name (str): Output filename for the figure
    """

    diff_matrix = contact_MUT - contact_WT

    # Only annotate non-zero differences
    annot_labels = np.where(diff_matrix != 0, np.round(diff_matrix, 2).astype(str), "")

    plt.figure(figsize=(12, 8))
    if annotate == True:
        sns.heatmap(diff_matrix, annot=annot_labels, fmt="", cmap="bwr", center=0,
                    xticklabels=resids_x, yticklabels=resids_y, vmin=vmin,vmax=vmax, 
                    cbar_kws={'label': 'Δ Contact Frequency (MUT - WT)'})
    else:
        sns.heatmap(diff_matrix, cmap="bwr", center=0,
                    xticklabels=resids_x, yticklabels=resids_y,vmin=vmin,vmax=vmax, 
                    cbar_kws={'label': 'Δ Contact Frequency (MUT - WT)'})

    # Reduce label density: show every 5th residue
    xticks = np.arange(0, len(resids_x), tick_label_spacing)
    yticks = np.arange(0, len(resids_y), tick_label_spacing)
    
    plt.xticks(ticks=xticks, labels=[resids_x[i] for i in xticks], rotation=90, ha='right', fontsize=8)
    
    plt.yticks(ticks=yticks, labels=[resids_y[i] for i in yticks], fontsize=8)
    
    plt.gca().invert_yaxis()
    plt.xlabel("Ig21 residues (CD-face)")
    plt.ylabel("Ig20 residues (β-strand)")
    plt.title("Contact Frequency Difference Map (MUT - WT)")
    plt.tight_layout()
    plt.savefig(output_name, dpi=300)
    # Ensure plot is shown in Jupyter
    from IPython.display import display
    display(plt.gcf())
    plt.close()

    return diff_matrix

def plot_contact_difference_interactive(contact_WT, contact_MUT, resids_x, resids_y,
                                        tick_label_spacing=5,
                                        vmin=-1, vmax=1,
                                        xlabel="Ig21 residues",ylabel="Ig20 residues",
                                        output_name="contact_diff_map.html"):
    """
    Interactive Plotly heatmap of contact frequency difference between MUT and WT.

    Args:
        contact_WT (np.ndarray): WT contact matrix
        contact_MUT (np.ndarray): MUT contact matrix
        resids_x (list): X-axis residues (e.g., Ig21 CD-face)
        resids_y (list): Y-axis residues (e.g., Ig20 β-strand)
        tick_label_spacing (int): Residue label spacing
        vmin (float): Min color scale
        vmax (float): Max color scale
        output_name (str): Output HTML file
    """

    diff_matrix = contact_MUT - contact_WT

    # Build tick labels
    xticks = [str(r) if i % tick_label_spacing == 0 else '' for i, r in enumerate(resids_x)]
    yticks = [str(r) if i % tick_label_spacing == 0 else '' for i, r in enumerate(resids_y)]

    # Reverse Y-axis for biological orientation
    #diff_matrix = diff_matrix[::-1]
    #yticks = yticks[::-1]

    fig = go.Figure(data=go.Heatmap(
        z=diff_matrix,
        x=resids_x,
        y=resids_y[::-1],
        text=np.round(diff_matrix[::-1], 2),
        hovertemplate="Ig20 resid: %{y}<br>Ig21 resid: %{x}<br>ΔFreq: %{z}<extra></extra>",
        colorscale="RdBu",
        zmin=vmin,
        zmax=vmax,
        colorbar=dict(title="Δ Contact Freq")
    ))

    fig.update_layout(
        title="Interactive Contact Frequency Difference Map (MUT - WT)",
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        xaxis=dict(tickmode='array', tickvals=resids_x[::tick_label_spacing], ticktext=xticks[::tick_label_spacing]),
        yaxis=dict(tickmode='array', tickvals=resids_y[::-1][::tick_label_spacing], ticktext=yticks[::tick_label_spacing]),
        autosize=False,
        width=800,
        height=600,
        margin=dict(l=80, r=80, t=100, b=80)
    )

    # Save to HTML
    pio.write_html(fig, file=output_name, auto_open=True)

    fig.show()

    return diff_matrix
