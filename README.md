# Estimación del contenido de agua precipitable a partir de simulaciones climáticas a alta resolución en la Península Ibérica
Resumen del trabajo: (EN PROCESO)

## Scripts

**read_cdfs.ipynb**: Notebook auxiliar para leer información de archivos NETCDFs.

### Generación de los datos de rocío recolectable
**netcdfs_dewyield_CERRA.py**: Script que genera un archivo NETCDF para cada año con los datos de rocío recolectable por hora obtenidos de los datos de CERRA cada seis horas entre 1991 y 2020.

**netcdfs_dewyield_WRF**: Script igual al anterior con los datos de WRF. Incluye la interpolación de los datos WRF para adaptarlos a la rejilla de CERRA (de 5km a 5,5km)

(faltan los atributos en los NETCDFs, si son necesarios se añadirán en el futuro)

### Estudio de los campos medios
**CERRA_WRF_comp.py**: Script con varias funciones:
- Genera un plot que compara la media de rocío recolectable al año en el periodo 1991-2020 de los datos CERRA y el modelo WRF.
- Calcula y genera un plot de la diferencia relativa de las medias de los modelos. Realiza también un test U de Mann-Whitney.
- Calcula las desviaciones típicas y genera un plot del cociente de las desviaciones de ambos modelos. Realiza también un test F de Snedecor.

## Figuras
**CERRA_WRF_dewyiled.png** Comparación lado a lado de la media de producción de rocío en el periodo 1991-2020 en el modelo de WRF y el reanálisis de CERRA.

**dif_relativas.png** Diferencias relativas de la media de producción de rocío en la PI entre el modelo WRF y los datos CERRA.

**coc_var.png** Mapa de cociente de varianzas