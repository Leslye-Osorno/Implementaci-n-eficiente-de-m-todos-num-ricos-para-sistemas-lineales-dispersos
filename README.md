# Implementación eficiente de métodos numéricos para sistemas lineales dispersos

Este repositorio contiene el código fuente y las librerías desarrolladas como parte de la tesis para obtener el título de Matemática en el Centro Universitario de Ciencias Exactas e Ingenierías (CUCEI) de la Universidad de Guadalajara.

El proyecto consiste en una librería escrita de forma modular en el lenguaje C, enfocada en resolver sistemas de ecuaciones lineales explotando de forma eficiente la estructura de las matrices dispersas para reducir tanto el consumo de memoria como el tiempo de cómputo.

## Características Principales

* **Formatos de Almacenamiento Disperso**: Implementación y gestión de estructuras de datos híbridas y nativas, abarcando Lista de Coordenadas (COO), Lista de Listas (LIL) y su conversión a *Compressed Sparse Row* (CSR) para cómputo de alto rendimiento.
* **Métodos Numéricos Implementados**:
    * **Factorización PLU**: Implementación directa con técnicas para matrices dispersas.
    * **Descomposición de Cholesky**: Aplicada a sistemas simétricos.
    * **Gradiente Conjugado**: Enfoque iterativo optimizado, evaluando su comportamiento frente al fenómeno de relleno (*fill-in*).
* **Telemetría y Rendimiento**: Capacidad para medir de forma automatizada tiempos de ejecución, consumo de memoria estática y dinámica en el heap (huella de memoria), y cálculo de las normas de residuales algebraicos.

## Arquitectura del Código

La librería está dividida en los siguientes módulos funcionales, garantizando independencia y escalabilidad:

* `main.c`: Interfaz principal mediante menús interactivos, inicialización de datos, cronometraje de ejecución y llamado a rutinas de exportación.
* `metodos_lineales.c` / `.h`: Contiene la lógica de los resolvedores lineales (PLU, Cholesky, Gradiente Conjugado) abstraídos mediante estructuras intermedias optimizadas.
* `operaciones_basicas.c` / `.h`: Núcleo matemático para el álgebra vectorial y el producto matriz-vector de alto rendimiento en formato CSR, así como rutinas de conversión.
* `leer_matriz.c` / `.h`: Mapeo dinámico y carga estructurada de matrices en formatos base desde archivos en disco.
* `tools.c` y `io_tools.c` / `.h`: Funciones de auditoría analítica del rendimiento, persistencia de variables métricas y volcado de resultados a bitácoras estructuradas (`.csv`).
* `errores.c` / `.h`: Prevención y manejo de fallos matemáticos (matrices asimétricas, diagonales críticas) o computacionales (falta de memoria dinámica).

## Autoría y Contexto Académico

* **Autora**: Leslye Valeria Osorno Amaya
* **Director de Tesis**: Abel Palafox González
* **Institución**: Universidad de Guadalajara, División de Ciencias Básicas
