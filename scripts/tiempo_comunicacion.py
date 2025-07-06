import matplotlib.pyplot as plt
# ================================
# Datos de tiempo de cómputo puro
# ================================
procesadores = [1, 2, 4, 8, 16, 32]
tamanios = [100000, 1000000, 10000000, 100000000, 1000000000]

# =============================
# Estilos visuales del gráfico
# =============================
estilos = ['-o', '-s', '-^', '-D', '-v', '-x']
colores = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

# Datos de tiempo de comunicación bloqueante
tiempo_comunicacion_block = {
    100000:     [0.000, 0.057, 0.393, 0.373, 1.380, 7.091],
    1000000:    [0.000, 2.397, 3.386, 4.794, 8.903, 9.304],
    10000000:   [0.000, 24.279, 34.429, 42.588, 48.897, 59.498],
    100000000:  [0.000, 247.034, 353.477, 430.062, 478.711, 603.618],
    1000000000: [0.000, 2495.966, 3613.272, 4434.876, 4887.489, 6246.657],
}

# Datos de tiempo de comunicación no bloqueante
tiempo_comunicacion_nonblock = {
    100000:     [0.000, 0.192, 0.292, 0.396, 0.690, 0.769],
    1000000:    [0.000, 1.996, 2.872, 3.586, 7.105, 7.614],
    10000000:   [0.000, 19.439, 27.686, 33.874, 38.524, 49.510],
    100000000:  [0.000, 192.657, 277.593, 347.216, 401.292, 489.215],
    1000000000: [0.000, 1901.970, 2811.093, 3522.987, 3893.973, 4895.487],
}

# Gráfico: Tiempo de Comunicación Bloqueante
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_comunicacion_block[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'N = {N}', linewidth=2, markersize=8)

plt.title('Tiempo de Comunicación Bloqueante vs Número de Procesadores', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Tiempo de Comunicación (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tamaño del problema (N)", fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.savefig('tiempo_comunicacion_block.png')
plt.show()

# Gráfico: Tiempo de Comunicación No Bloqueante
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_comunicacion_nonblock[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'N = {N}', linewidth=2, markersize=8)

plt.title('Tiempo de Comunicación No Bloqueante vs Número de Procesadores', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Tiempo de Comunicación (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tamaño del problema (N)", fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.savefig('tiempo_comunicacion_nonblock.png')
plt.show()

# Gráfico comparativo: Comunicación Bloqueante vs No Bloqueante
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_comunicacion_block[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'Block - N={N}', linewidth=2, markersize=6)
    plt.plot(procesadores, tiempo_comunicacion_nonblock[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], linestyle='--', label=f'NonBlock - N={N}', linewidth=2, markersize=6)

plt.title('Comparación de Tiempo de Comunicación Bloqueante vs No Bloqueante', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Tiempo de Comunicación (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tipo y Tamaño del problema", fontsize=9, title_fontsize=10, loc='upper left')
plt.tight_layout()
plt.savefig('tiempo_comunicacion_vs.png')
plt.show()
