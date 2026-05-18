# DESCARTES Exp2 microbiome

Proyecto Quarto autocontenido para revisar y publicar los resultados de microbiota 16S del experimento `descartes_exp2`.

Este repositorio contiene directamente los capitulos, figuras y tablas que usa la web. No depende de scripts, reports ni rutas externas del proyecto de analisis para renderizarse o publicarse.

## Trabajo local

Editar el contenido de la web directamente en:

```text
index.qmd
chapters/*.qmd
appendices/reproducibility.qmd
styles.css
_quarto.yml
```

Renderizar la web con:

```bash
quarto render
```

El sitio HTML se genera en `docs/`. Esa carpeta es output de render y no debe editarse a mano.

## Estructura

- `_quarto.yml`: configuracion del libro Quarto.
- `index.qmd`: pagina de inicio.
- `chapters/`: capitulos editables del libro.
- `appendices/`: material complementario y reproducibilidad.
- `assets/results/`: figuras `PNG` y tablas `CSV` usadas por los capitulos.
- `styles.css`: estilos visuales de la web.
- `.github/workflows/pages.yml`: despliegue automatico en GitHub Pages.
- `references.bib`: bibliografia citable versionada.
- `bibliography/`: PDFs locales de articulos; esta carpeta esta ignorada por Git y no se sube a GitHub.

## Bibliografia local

Guarda los PDFs de articulos en `bibliography/` solo para trabajo local. Las citas que deban formar parte del proyecto deben incorporarse a `references.bib`, que si esta versionado.

## Publicacion

Cualquier push a `main` dispara el workflow de GitHub Pages. La web publica se sirve desde:

```text
https://dcosimber.github.io/descartes_exp2/
```
