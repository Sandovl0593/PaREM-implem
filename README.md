<h1 align="center">PaREM-implem</h1>

#### Implementación de _Parallel Regular Expression Matching_ con MPI y OpenMP en C++.

Descripción de la implementación y comparación de mediciones con la complejidad teórica y experimental de la implementación bajo memoria distribuida y compartida.

Se toma en cuenta que no implementamos de Expresiones Regulares a AFD's según el paper por complejidad de código, sino directamente solo AFD como clase resultante a partir de una expresión regular y que se instancie directamente para test del proyecto.

## Compilación y Entorno

- En secuencial

```bash
g++ seciencial.cpp -o secuencial
```

- Con MPI
- Con OMP   

## Integrantes

- Arleth Ivhy Lapa Carhuamaca
- Adrian Sandoval Huamani
- Salvador Eliot Hilares Barrios
