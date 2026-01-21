import pandas as pd
import matplotlib.pyplot as plt

def plot_final():
    # Use classic style
    plt.style.use('classic')
    plt.figure(figsize=(12, 7))
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Load Data
    try:
        df_base = pd.read_csv('baseline_results.csv')
        df_hyb = pd.read_csv('hybrid_stable.csv')
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    # Filter valid numeric data
    df_base = df_base[pd.to_numeric(df_base['Step'], errors='coerce').notnull()]
    df_hyb = df_hyb[pd.to_numeric(df_hyb['Step'], errors='coerce').notnull()]
    
    df_base['Step'] = df_base['Step'].astype(float)
    df_hyb['Step'] = df_hyb['Step'].astype(float)
    df_base['AvgEnergy'] = df_base['AvgEnergy'].astype(float)
    df_hyb['AvgEnergy'] = df_hyb['AvgEnergy'].astype(float)

    # Truncate to common length
    max_step = min(df_base['Step'].max(), df_hyb['Step'].max())
    df_base = df_base[df_base['Step'] <= max_step]
    df_hyb = df_hyb[df_hyb['Step'] <= max_step]

    # Plot Lines
    plt.plot(df_base['Step'], df_base['AvgEnergy'], label='Baseline (Local MCMC)', color='blue', alpha=0.8, linewidth=2)
    plt.plot(df_hyb['Step'], df_hyb['AvgEnergy'], label='Hybrid (NF + MCMC)', color='red', alpha=0.8, linewidth=2)
    
    # Identify Accepted Jumps
    # Check for both int 1 and string 'True'
    jumps = df_hyb[(df_hyb['Accepted'] == 1) | (df_hyb['Accepted'] == 'True') | (df_hyb['Accepted'] == 'ACC')]
    
    if len(jumps) > 0:
        plt.scatter(jumps['Step'], jumps['AvgEnergy'], color='green', marker='*', s=150, label='Accepted NF Jump', zorder=10, edgecolors='black')

    # Formatting
    plt.xlabel('Simulation Step', fontsize=12)
    plt.ylabel('Total Energy (Lower is Better)', fontsize=12)
    plt.title('Convergence Comparison: Baseline vs Hybrid NF-MCMC', fontsize=14, fontweight='bold')
    plt.legend(fontsize=10)
    plt.tight_layout()
    
    # Save
    plt.savefig('energy_comparison_final.png', dpi=100)
    print("Saved energy_comparison_final.png")

if __name__ == "__main__":
    plot_final()
