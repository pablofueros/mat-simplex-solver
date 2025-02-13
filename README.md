# Solver del Método Simplex en MATLAB

Este proyecto implementa el **método simplex** para la resolución de problemas de **programación lineal** en MATLAB. Fue desarrollado como **trabajo final de la asignatura de Optimización 1** en el grado de Matemáticas, en colaboración con **Ángel**, un excompañero de clase. Además, representa uno de mis primeros proyectos completos en los que diseñé un sistema modular que integra varios módulos y gestiona la ejecución en función de un input.

## 🚀 Características
- Implementación del **método simplex revisado** con actualización eficiente mediante **factorización QR**.
- Resolución de problemas de programación lineal con variables artificiales mediante un enfoque de **dos fases**.
- Código modular y estructurado para facilitar su uso y modificación.
- Ejemplos prácticos de aplicación a problemas reales.

## 📁 Estructura del Proyecto
```
mat-simplex-solver/
│── src/
│   ├── qr_add_column.m      % Añadir una columna a la factorización QR
│   ├── qr_remove_column.m   % Eliminar una columna de la factorización QR
│   ├── simplex_phase1.m     % Implementación de la Fase I del método simplex
│   ├── simplex_phase2.m     % Implementación de la Fase II del método simplex
│   ├── simplex_solver.m     % Función principal para resolver problemas LP
│── examples/
│   ├── Dieta.m              % Ejemplo de optimización de dieta
│   ├── Prueba_1.m           % Prueba básica del método simplex
│   ├── Prueba_2.m           % Otra prueba con distintos parámetros
│   ├── Transporte.m         % Ejemplo de optimización de transporte
│── README.md                % Este archivo de documentación
```

## 🛠️ Instalación y Uso
1. Clona este repositorio o descarga los archivos.
   ```bash
   git clone https://github.com/pablofueros/mat-simplex-solver.git
   ```
2. Abre MATLAB y navega hasta la carpeta `src`.
3. Llama a `simplex_solver.m` con los parámetros adecuados para resolver un problema de programación lineal.
4. Para ver ejemplos de uso, revisa los archivos en la carpeta `examples` y ejecútalos en MATLAB.

## 📜 Descripción de los Archivos

### 🔹 Factorización QR
- **`qr_add_column.m`**: Añade una columna a la factorización QR y actualiza las matrices correspondientes.
- **`qr_remove_column.m`**: Elimina una columna de la factorización QR y ajusta las matrices \( Q \) y \( R \).

### 🔹 Método Simplex
- **`simplex_phase1.m`**: Encuentra una solución factible inicial si el problema lo requiere.
- **`simplex_phase2.m`**: Optimiza la función objetivo a partir de una solución factible.
- **`simplex_solver.m`**: Función principal que coordina la ejecución del método simplex.

### 🔹 Ejemplos de Aplicación
- **`Dieta.m`**: Modelo de optimización de una dieta equilibrada.
- **`Prueba_1.m`** y **`Prueba_2.m`**: Casos de prueba para validar el método simplex.
- **`Transporte.m`**: Ejemplo de optimización de costos de transporte.

## 📌 Notas
- Este código fue desarrollado en el contexto académico y puede requerir ajustes.
- Se recomienda tener conocimientos previos sobre el método simplex y factorización QR para entender mejor su funcionamiento.

## 📖 Créditos
Desarrollado por **[Pablo Garcia](https://github.com/pablofueros)** y **Ángel de Hoyos** como parte del trabajo final de la asignatura de **Optimización 1** en el grado de Matemáticas.

## 🏆 Licencia
Este proyecto se distribuye bajo la licencia MIT. ¡Siéntete libre de usarlo, modificarlo y compartirlo! 🎯

