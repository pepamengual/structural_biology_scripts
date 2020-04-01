import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def load_traj(trajectory_file, topology_file):
    traj = md.load(trajectory_file, top=topology_file)
    return traj

def selecting_atoms(traj, selection):
    atom_selection = traj.topology.select(selection)
    return atom_selection

def defining_traj_atom_selection(traj, atom_selection):
    traj_selection = traj.atom_slice(atom_selection)
    return traj_selection

def compute_rmsd(traj_selection):
    rmsd = md.rmsd(traj_selection, traj_selection)
    return rmsd

def compute_rmsf(traj_selection):
    traj_selection.center_coordinates() #docs says its faster in this way
    rmsf = md.rmsf(traj_selection, traj_selection, precentered=True) #precentered only if center_coodinates
    return rmsf

def mean_rmsf_by_residue(traj_selection, rmsf):
    topology_data = traj_selection.topology

    rmsf_all_atoms_residue = {}
    for i, index in enumerate(topology_data.atoms):
        residue = str(index).split("-")[0]
        rmsf_all_atoms_residue.setdefault(residue, []).append(rmsf[i])
    
    rmsf_mean_residue = {}
    for residue, data in rmsf_all_atoms_residue.items():
        rmsf_mean_residue.setdefault(residue, np.mean(data)*10)
    return rmsf_mean_residue

def plot_trimer(rmsf_mean_residue_1, rmsf_mean_residue_2, position, mutation, trajectory_number):
    split_values_1 = np.array_split(list(rmsf_mean_residue_1.values()), 3)
    split_values_2 = np.array_split(list(rmsf_mean_residue_2.values()), 3)
    r = list(range(position-30,position+30))
    max_value = 0

    fig, axs = plt.subplots(3)   
    for i, array in enumerate(split_values_1):
        rmsf_subset_1 = list(array[r[0]:r[-1]])
        rmsf_subset_2 = list(split_values_2[i][r[0]:r[-1]])
        max_subset = max(rmsf_subset_1 + rmsf_subset_2)
        if max_subset > max_value:
            max_value = max_subset
        axs[i].plot(r[0:-1], rmsf_subset_1, label="Native Spike")
        axs[i].plot(r[0:-1], rmsf_subset_2, label="{} Spike".format(mutation))

    for i, ax in enumerate(axs.flat):
        ax.set_ylim(0, max_value)
        if i == 0:
            ax.legend()
        if i == 1:
            ax.set(ylabel="RMSF (Å)")
        if i == 2:
            ax.set(xlabel="Positions")

    plt.tight_layout()
    plt.savefig("{}_{}_all".format(mutation, trajectory_number), dpi=300)

def iterator(trajectory_file, topology_file, selection):
    traj = load_traj(trajectory_file, topology_file)
    atom_selection = selecting_atoms(traj, selection)
    traj_selection = defining_traj_atom_selection(traj, atom_selection)
    rmsf = compute_rmsf(traj_selection)
    rmsf_mean_residue = mean_rmsf_by_residue(traj_selection, rmsf)
    return rmsf_mean_residue

def main():
    #selection = "protein and name CA"
    selection = "protein"
    for i in range(1,5):
        trajectory_file = "covid19_trimer/trajectory_{}.xtc".format(i)
        topology_file = "covid19_trimer/topology_0.pdb"
        rmsf_mean_residue_1 = iterator(trajectory_file, topology_file, selection)
        
        mutation = "G387A"
        trajectory_file = "covid19_trimer_{}/trajectory_{}.xtc".format(mutation, i)
        topology_file = "covid19_trimer_{}/topology_0.pdb".format(mutation)
        position = int(mutation[1:-1])
        rmsf_mean_residue_2 = iterator(trajectory_file, topology_file, selection)
        plot_trimer(rmsf_mean_residue_1, rmsf_mean_residue_2, position, mutation, i)
        print("--> Trajectory {}_{} finished".format(mutation, i))
main()
