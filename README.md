# Sensibilización molecular a Der p (ácaro del polvo doméstico)

**pMenor · Estadística clínica en R**

Código del análisis estadístico del artículo:

> Martín-Martínez MD, García-de la Rosa G, De Castro-Córdova J. Sensitization
> patterns to molecular allergens of house dust mites (*Dermatophagoides
> pteronyssinus*) on the Canary Islands. *Allergo J Int.* 2025.
> [doi:10.1007/s40629-025-00349-5](https://doi.org/10.1007/s40629-025-00349-5)

## La pregunta

En la alergia al ácaro *D. pteronyssinus*, Der p 1 y Der p 2 se consideran los
alérgenos mayores, y en ellos se basan buena parte del diagnóstico y de la
inmunoterapia. Pero el perfil de sensibilización cambia con el clima y la
exposición. En una población con exposición alta y constante a ácaros, ¿qué
componentes predominan, cómo cambian con la edad y hay perfiles de pacientes
diferenciados?

## El análisis

1.685 pacientes estudiados con la plataforma multiplex ALEX2, con IgE
específica frente a 9 componentes de Der p. Sensibilización: IgE ≥ 0,35 kUA/L.

| Pregunta | Método |
|---|---|
| ¿Qué componentes predominan? | Prevalencia de sensibilización por componente |
| ¿Cambia con la edad? | χ² o Fisher por componente; Kruskal-Wallis y Dunn para niveles de IgE |
| ¿Más componentes, más IgE total? | Correlación de Spearman |
| ¿Qué componentes van juntos? | Gráfico UpSet, mapa de calor y análisis de componentes principales |
| ¿Hay perfiles de pacientes? | k-means (k = 3) y su distribución por edad |

Para comprobar qué permitía la muestra se calculan la **precisión** de las
prevalencias y la **diferencia mínima detectable** entre grupos de edad, no la
potencia a posteriori, que depende del resultado obtenido.

La edad solo constaba en 582 pacientes. Los análisis por edad usan esos 582;
el resto, la cohorte completa.

### Resultados publicados

| Medida | Resultado |
|---|---|
| Sensibilizados a ≥ 1 componente | 62,3 % |
| Componentes más frecuentes | Der p 20 (46,4 %), Der p 7 (45,9 %), Der p 1 (45,7 %) |
| Der p 2 | 7,3 % (23,3 % en niños de 6-11 años; 4,6 % en adultos) |
| IgE total frente a n.º de componentes | ρ = 0,50 (p < 0,0001) |
| Grupo de edad con más sensibilización | Adolescentes (12-17 años) |
| Perfiles (k-means) | Polisensibilizados con IgE alta; mínimamente sensibilizados; intermedios. Su distribución cambia con la edad |

**Conclusión del artículo:** en esta población predominan componentes
considerados "menores" (Der p 20, 7 y 1), mientras que Der p 2 es poco
frecuente. Esto apoya un diagnóstico molecular y una inmunoterapia adaptados a
cada región.

## Datos

Este repositorio **no contiene datos de pacientes**. `R/simular_datos.R` genera
un dataset ficticio con la misma estructura (`data/der_p_simulado.csv`), con las
prevalencias por grupo de edad y el patrón de co-sensibilización del estudio.
Las cifras se parecen a las publicadas, pero no son las mismas.

## Cómo reproducirlo

1. Clona el repositorio y abre la carpeta como proyecto de RStudio.
2. Instala los paquetes con las mismas versiones que se usaron (quedan
   registradas en `renv.lock`):

   ```r
   install.packages("renv")
   renv::restore()
   ```

3. Ejecuta `rmarkdown::render("analisis.Rmd")`.

El resultado, con código, salidas y figuras, se puede leer directamente en
[`analisis.md`](analisis.md). Las figuras usan las fuentes IBM Plex Sans e IBM
Plex Mono; si no están instaladas, R usa otras y avisa.

## Estructura

```
├── analisis.Rmd          análisis completo
├── analisis.md           el mismo análisis ya ejecutado (se lee en GitHub)
├── R/
│   ├── simular_datos.R   genera los datos ficticios
│   └── theme_pmenor.R    estilo de los gráficos
├── data/                 datos simulados
└── output/               figuras
```

## Licencia

El código se publica con licencia [MIT](LICENSE): puedes usarlo y adaptarlo citando la fuente. Si lo usas en un trabajo, cita también el artículo.

---

**pMenor** · Estadística clínica en R · [pmenor.es](https://pmenor.es)
