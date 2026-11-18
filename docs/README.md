# Documentation

`theory/` holds the theory notes as Markdown with LaTeX math, one file per
chapter (read them in order: 00, 01, ...). `figures/` holds the figures and
the script that generates them from ENDF/B-VIII.1 with the project's code.

```sh
make docs      # docs/build/theory.html: one self-contained file, works offline
make figures   # regenerate figures/*.svg (needs `make data`; uses cached tables)
```

The HTML uses MathML (rendered natively by current browsers, no network
needed). For MathJax instead: `make docs DOCS_MATH=--mathjax` (needs network
when viewing).
