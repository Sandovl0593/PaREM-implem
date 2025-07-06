import matplotlib.pyplot as plt

# ================================
# Datos
# ================================
procesadores = [1, 2, 4, 8, 16, 32]

# Eficiencia del código original
ef_block = [1.000, 0.519, 0.269, 0.131, 0.065, 0.031]
ef_nonblock = [1.000, 0.577, 0.315, 0.160, 0.081, 0.038]

# Eficiencia del código optimizado
ef_block_opti = [1.000, 0.942, 0.910, 0.908, 0.832, 0.392]
ef_nonblock_opti = [1.000, 0.983, 0.982, 0.985, 0.958, 0.462]

# ================================
# Visualización
# ================================
plt.figure(figsize=(12, 7))

plt.plot(procesadores, ef_block, '-o', label='Original - Bloqueante', color='tab:red', linewidth=2, markersize=8)
plt.plot(procesadores, ef_nonblock, '-s', label='Original - No Bloqueante', color='tab:blue', linewidth=2, markersize=8)
plt.plot(procesadores, ef_block_opti, '-o', label='Optimizado - Bloqueante', color='tab:orange', linewidth=2, markersize=8)
plt.plot(procesadores, ef_nonblock_opti, '-s', label='Optimizado - No Bloqueante', color='tab:green', linewidth=2, markersize=8)

plt.title('Comparación de Eficiencia vs Número de Procesadores (N = 10⁷)', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Eficiencia', fontsize=12)
plt.ylim(0, 1.1)
plt.xticks(procesadores)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("opti_comparacion_eficiencia.png")
plt.show()

# ================================
# Comparación Speedup Original vs Optimizado (sin línea teórica)
# ================================

# Speedup original (sin optimización)
speedup_block_orig = [1.0, 1.037, 1.078, 1.050, 1.044, 1.005]
speedup_nonblock_orig = [1.0, 1.154, 1.261, 1.279, 1.304, 1.201]

# Speedup optimizado
speedup_block_opti = [1.0, 1.689, 1.130, 1.050, 1.044, 1.005]
speedup_nonblock_opti = [1.0, 1.172, 1.272, 1.279, 1.304, 1.201]

plt.figure(figsize=(12, 7))
plt.plot(procesadores, speedup_block_orig, '-o', label='Original Bloqueante', color='tab:red', linewidth=2, markersize=8)
plt.plot(procesadores, speedup_nonblock_orig, '-s', label='Original No Bloqueante', color='tab:blue', linewidth=2, markersize=8)
plt.plot(procesadores, speedup_block_opti, '-^', label='Optimizado Bloqueante', color='tab:orange', linewidth=2, markersize=8)
plt.plot(procesadores, speedup_nonblock_opti, '-D', label='Optimizado No Bloqueante', color='tab:green', linewidth=2, markersize=8)

plt.title('Comparación de Speedup: Original vs Optimizado (N=10⁷)', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Speedup (sin unidades)', fontsize=12)
plt.xticks(procesadores)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('opti_vs_original_speedup.png')
plt.show()
