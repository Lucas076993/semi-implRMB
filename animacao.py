import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Configuração geral
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'  # ajuste se necessário

# Arquivos de entrada
files = {
    'R': 'R.csv',
    'M': 'M.csv',
    'B': 'B.csv'
}

# Leitura e pré-processamento
print("Lendo arquivos CSV...")
dfs = {}
for var, fname in files.items():
    df = pd.read_csv(fname)
    col_name = df.columns[3]  # assume que a quarta coluna é o valor
    dfs[var] = df

# Determinar a malha espacial (supõe-se que seja a mesma para todos os tempos)
df_ex = dfs['R']
t_vals = sorted(df_ex['t'].unique())
x_vals = sorted(df_ex['x'].unique())
y_vals = sorted(df_ex['y'].unique())
nx, ny = len(x_vals), len(y_vals)

# Função para construir matriz 2D para um dado t e variável
def get_frame(df, t, nx, ny):
    """Retorna matriz 2D (ny, nx) com os valores no instante t."""
    subset = df[df['t'] == t]
    grid = np.zeros((ny, nx))
    for _, row in subset.iterrows():
        i = x_vals.index(row['x'])
        j = y_vals.index(row['y'])
        # Queremos grid[j, i] se y é linha, x é coluna
        grid[j, i] = row.iloc[3]  # quarta coluna é o valor
    return grid

# Preparar dados iniciais
t0 = t_vals[0]
frames_data = []
for var in ['R', 'M', 'B']:
    frames_data.append(get_frame(dfs[var], t0, nx, ny))

vmins = {}
vmaxs = {}
for var in ['R', 'M', 'B']:
    all_vals = dfs[var].iloc[:, 3]  # quarta coluna
    vmins[var] = all_vals.min()
    vmaxs[var] = all_vals.max()

# Criar figura e eixos
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Simulação NSFD – Modelo Espacial', fontsize=14)

# Listas para armazenar os objetos de imagem e as barras de cor
ims = []
cbs = []
titles = ['Rês', 'Moscas', 'Besouros']
cmaps = ['Blues', 'Oranges', 'Greens']  # cores sugestivas

for ax, data, title, cmap, var in zip(axes, frames_data, titles, cmaps, ['R','M','B']):
    im = ax.imshow(data, origin='lower', extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]],
                   cmap=cmap, vmin=vmins[var], vmax=vmaxs[var], aspect='auto')
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    # Barra de cor
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cb = fig.colorbar(im, cax=cax)
    ims.append(im)
    cbs.append(cb)

# Texto para exibir o tempo
time_text = fig.text(0.5, 0.02, f't = {t0:.3f}', ha='center', fontsize=12)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# Função de atualização para animação
def update(frame_idx):
    t = t_vals[frame_idx]
    for i, var in enumerate(['R', 'M', 'B']):
        grid = get_frame(dfs[var], t, nx, ny)
        ims[i].set_array(grid)
    time_text.set_text(f't = {t:.3f}')
    return ims + [time_text]

# Criar animação
ani = animation.FuncAnimation(fig, update, frames=len(t_vals), interval=100, blit=True)

# Salvar como vídeo
print("Salvando animação em 'simulacao.mp4'...")
ani.save('simulacao.mp4', writer='ffmpeg', fps=10)

