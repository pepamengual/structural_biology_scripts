import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def load_traj(pdb_file):
    pdb = md.load(pdb_file)
    return pdb

def compute_sasa(pdb):
    sasa = md.shrake_rupley(pdb, mode="residue")
    return sasa[0]

def plot_sasa(sasa_cov, sasa_cov_2):
    plt.plot(sasa_cov, label="SARS-CoV")
    plt.plot(sasa_cov_2, label="SARS-CoV-2")
    plt.legend()
    plt.show()

def chunker(sasa, chunk_size):
    len_sasa = int(len(sasa)/3)
    sasa_monomer = list(sasa)[0:len_sasa]
    chunks = int(len_sasa/chunk_size)
    chunks_sasa_monomer = np.array_split(sasa_monomer, chunks)
    sasa_each = []
    for array in chunks_sasa_monomer:
        mean = np.mean(array)
        sasa_each.append(mean)
    return sasa_each

def plot_sasa_monomer(sasa_each, name, chunk_size):
    plt.plot(sasa_each, label=name, linestyle='-', marker='D')
    plt.ylabel('Total SASA (nm)^2', size=12)
    plt.xlabel('Every {} residues in closed Spike conformation'.format(chunk_size), size=12)
    plt.tight_layout()
    plt.legend()
    plt.savefig("sasa_{}.png".format(name), dpi=300)
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure
    plt.close() # Close a figure window

def plot_trimer(sasa_cov, sasa_cov_2, chunk_size):
    plt.plot(sasa_cov, label="SARS-CoV", linestyle='-', marker='D')
    plt.plot(sasa_cov_2, label="SARS-CoV-2", linestyle='-', marker='D')
    plt.ylabel('Total SASA (nm)^2', size=12)
    plt.xlabel('Every {} residues in closed Spike conformation'.format(chunk_size), size=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig("sasa_trimer.png", dpi=300)

def main():
    cov_closed  = "sars_cov_trimmer_6acc_hbond-opt.pdb"
    cov_2_closed = "sars_cov2_trimmer_swiss_6vxx_preprocessed_hbond-opt.pdb"
    
    cov_closed_pdb = load_traj(cov_closed)
    cov_2_closed_pdb = load_traj(cov_2_closed)
    
    sasa_cov_closed = compute_sasa(cov_closed_pdb)
    sasa_cov_2_closed = compute_sasa(cov_2_closed_pdb)
    
    chunk_size = 10
    sasa_cov = chunker(sasa_cov_closed, chunk_size)
    sasa_cov_2 = chunker(sasa_cov_2_closed, chunk_size)

    plot_sasa_monomer(sasa_cov, "SARS-CoV", chunk_size)
    plot_sasa_monomer(sasa_cov_2, "SARS-CoV-2", chunk_size)

    plot_trimer(sasa_cov, sasa_cov_2, chunk_size)
main()
