import pandas as pd
import matplotlib.pyplot as plt

def load_data(filename):
    try:
        df = pd.read_csv(filename)
        # Ensure Step is numeric
        df = df[pd.to_numeric(df['Step'], errors='coerce').notnull()]
        for col in ['Step', 'AvgEnergy', 'AvgShift', 'AvgAcc', 'AvgMSD']:
            if col in df.columns:
                df[col] = df[col].astype(float)
        return df
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def get_accepted_steps(df):
    if df is None or 'Accepted' not in df.columns:
        return []
    # Check for various True formats
    mask = (df['Accepted'] == 1) | \
           (df['Accepted'] == 'True') | \
           (df['Accepted'] == 'ACC')
    return df[mask]

def plot_separate():
    df_base = load_data('baseline_results.csv')
    df_pre = load_data('pretrain_v4.csv')
    df_post = load_data('hybrid_test_results.csv')
    
    datasets = [
        ('Baseline', df_base, 'blue'),
        ('Pre-train', df_pre, 'green'),
        ('Post-train', df_post, 'red')
    ]
    
    metrics = [
        ('AvgAcc', 'Average Acceptance Ratio', 'plot_acc.png'),
        ('AvgEnergy', 'Average Energy', 'plot_energy.png'),
        ('AvgShift', 'Average Shift', 'plot_shift.png'),
        ('AvgMSD', 'Mean Squared Displacement', 'plot_msd.png')
    ]
    
    plt.style.use('classic')
    
    for metric, title, filename in metrics:
        plt.figure(figsize=(10, 6))
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.title(title, fontweight='bold', fontsize=14)
        plt.xlabel('Step', fontsize=12)
        plt.ylabel(metric, fontsize=12)
        
        for name, df, color in datasets:
            if df is not None and metric in df.columns:
                # Plot Line
                plt.plot(df['Step'], df[metric], label=name, color=color, alpha=0.8, linewidth=2)
                
                # Plot Stars for Acceptance (Only for Pre-train and Post-train to avoid clutter)
                if name in ['Pre-train', 'Post-train']:
                    accepted = get_accepted_steps(df)
                    if len(accepted) > 0:
                        plt.scatter(accepted['Step'], accepted[metric], 
                                    color=color, marker='*', s=200, 
                                    edgecolors='black', zorder=10, 
                                    label=f'{name} Accepted')
        
        plt.legend(fontsize=10)
        plt.tight_layout()
        plt.savefig(filename, dpi=150)
        print(f"Saved {filename}")
        plt.close()

if __name__ == "__main__":
    plot_separate()
