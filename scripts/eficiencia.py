import matplotlib.pyplot as plt

# =========================
# Datos
# =========================
procesadores = [1, 2, 4, 8, 16, 32]
tamanios = [100000, 1000000, 10000000, 100000000, 1000000000]

ef_block = {
    100000:     [1.000, 0.845, 0.270, 0.143, 0.029, 0.003],
    1000000:    [1.000, 0.530, 0.282, 0.124, 0.042, 0.022],
    10000000:   [1.000, 0.519, 0.269, 0.131, 0.065, 0.031],
    100000000:  [1.000, 0.526, 0.265, 0.130, 0.067, 0.033],
    1000000000: [1.000, 0.531, 0.265, 0.130, 0.065, 0.032],
}

ef_nonblock = {
    100000:     [1.000, 0.586, 0.336, 0.136, 0.057, 0.026],
    1000000:    [1.000, 0.578, 0.318, 0.159, 0.052, 0.026],
    10000000:   [1.000, 0.577, 0.315, 0.160, 0.081, 0.038],
    100000000:  [1.000, 0.593, 0.316, 0.157, 0.079, 0.040],
    1000000000: [1.000, 0.606, 0.318, 0.158, 0.080, 0.041],
}

# =========================
# Visualización
# =========================
estilos = ['-o', '-s', '-^', '-D', '-v']
colores = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']

# --- Eficiencia Bloqueante ---
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios[:-1]):
    plt.plot(procesadores, ef_block[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'Block - N={N}', linewidth=2, markersize=8)

plt.title('Eficiencia Bloqueante vs Número de Procesadores', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Eficiencia', fontsize=12)
plt.ylim(0, 1.05)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tamaño del problema (N)", fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.savefig('eficiencia_block.png')
plt.show()

# --- Eficiencia No Bloqueante ---
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios[:-1]):
    plt.plot(procesadores, ef_nonblock[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'NonBlock - N={N}', linewidth=2, markersize=8)

plt.title('Eficiencia No Bloqueante vs Número de Procesadores', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Eficiencia', fontsize=12)
plt.ylim(0, 1.05)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tamaño del problema (N)", fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.savefig('eficiencia_nonblock.png')
plt.show()
