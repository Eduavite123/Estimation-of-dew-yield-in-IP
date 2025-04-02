# Estimación del contenido de agua precipitable a partir de simulaciones climáticas a alta resolución en la Península Ibérica
Resumen del trabajo: (EN PROCESO)

## Scripts

### Generación de los datos de rocío recolectable
**netcdfs_dewyield_CERRA.py**: Script que genera una archivo NETCDF para cada año con los datos de rocío recolectable por hora obtenidos de los datos de CERRA cada seis horas entre 1991 y 2020.
**netcdfs_dewyield_WRF**: Script igual al anterior con los datos de WRF.
(faltan los atributos en los NETCDFs, si son necesarios se añadirán en el futuro)

### Estudio de los campos medios
**CERRA_WRF_comp.py**: (EN PROCESO) Script que genera un plot que compara la media de rocío recolectable al año en el periodo 1991-2020 de los datos CERRA y el modelo WRF. Incluye la interpolación de los datos WRF para adaptarlos a la rejilla de CERRA (de 5km a 5,5km) y la generación de archivos NETCDFs con las medias correspondientes para su futuro análisis.