import pandas as pd
import matplotlib.pyplot as plt

def load_data(filename):
    try:
        df = pd.read_csv(filename)
        # Ensure Step is numeric
        df = df[pd.to_numeric(df['Step'], errors='coerce').notnull()]
        for col in ['Step', 'AvgEnergy', 'AvgShift', 'AvgAcc', 'AvgMSD', 'FlowLoss']:
            if col in df.columns:
                df[col] = df[col].astype(float)
        return df
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def plot_metrics():
    df_base = load_data('baseline_results.csv')
    df_pre = load_data('pretrain_v4.csv')
    df_post = load_data('hybrid_test_results.csv')
    
    # Define datasets for loop
    datasets = [
        ('Baseline', df_base, 'blue'),
        ('Pre-train', df_pre, 'green'),
        ('Post-train', df_post, 'red')
    ]
    
    # 1. Main Metrics Plot (2x2)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    # Style
    plt.style.use('classic')
    
    metrics = [
        ('AvgAcc', 'Average Acceptance Ratio', axes[0,0]),
        ('AvgEnergy', 'Average Energy', axes[0,1]),
        ('AvgShift', 'Average Shift', axes[1,0]),
        ('AvgMSD', 'Mean Squared Displacement', axes[1,1])
    ]
    
    for metric, title, ax in metrics:
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel('Step')
        
        for name, df, color in datasets:
            if df is not None and metric in df.columns:
                ax.plot(df['Step'], df[metric], label=name, color=color, alpha=0.8, linewidth=1.5)
        
        ax.legend()

    plt.tight_layout()
    plt.savefig('comparison_metrics.png', dpi=150)
    print("Saved comparison_metrics.png")
    
    # 2. Loss Plot (Separate)
    plt.figure(figsize=(10, 6))
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.title('Flow Loss Comparison (Pre vs Post)', fontweight='bold')
    plt.xlabel('Step')
    plt.ylabel('Loss (-log p)')
    
    # Only Pre and Post have Loss (Baseline has 0 or NaN)
    if df_pre is not None and 'FlowLoss' in df_pre.columns:
        plt.plot(df_pre['Step'], df_pre['FlowLoss'], label='Pre-train Loss', color='green', linewidth=1.5)
    
    if df_post is not None and 'FlowLoss' in df_post.columns:
        plt.plot(df_post['Step'], df_post['FlowLoss'], label='Post-train Loss', color='red', linewidth=1.5)

    plt.legend()
    plt.tight_layout()
    plt.savefig('comparison_loss.png', dpi=150)
    print("Saved comparison_loss.png")

if __name__ == "__main__":
    plot_metrics()
