
# Repository Guidelines

## Project Structure & Purpose <!-- {{{* -->

`main.tex` is the main manuscript, using the `revtex4-2`
document class for submission to *Physical Review A* (American Physical Society).

The repository includes:
- `main.tex`: Main LaTeX manuscript.
- `references.bib`: BibTeX bibliography database.
- `Figures/`: Figures and supporting plotting scripts.
- `jadeleon.sty`: Custom LaTeX style definitions.
- `main.pdf`: Compiled manuscript.
- `Makefile`: Existing build configuration.

The project is about the dynamics of the reduced density matrix and the quantum channel
induced by the unitary evolution of the XXZ chain with two additional spins attached at the edge with two kinds of couplings, ising coupling and exchange coupling.

<!-- *}}} -->
## LaTeX Compilation <!-- {{{* -->

Use MiKTeX on Windows or TeX Live on Linux. From the repository root,
compile the manuscript with this BibTeX-based sequence:

1. `pdflatex -synctex=1 -interaction=nonstopmode -halt-on-error main.tex`
2. `bibtex main`
3. `pdflatex -synctex=1 -interaction=nonstopmode -halt-on-error main.tex`
4. `pdflatex -synctex=1 -interaction=nonstopmode -halt-on-error main.tex`

<!-- *}}} -->
## LaTeX Style <!-- {{{* -->

Preserve the `revtex4-2` structure and use the macros defined in
`jadeleon.sty` whenever applicable. Use descriptive equation labels,
such as `\label{eq:counterexample}`, and the package's reference macros:
`\Eref`, `\Sref`, `\Fref`, and `\Tref`. Use `\ifunam` and `\unam` for 
the corresponding affiliations.

Keep collaborator annotations in their provided commands, such as
`\janote`, `\mgnote`, `\pvnote`, `\codexnote`, `\jacomment`, and
`\fxnote`; use `\DeclareCommenter` only when adding a collaborator.
Do not call the internal helper macros `\flexfxnote` or
`\fxnoteoriginal` directly.

When asked to review comments in the paper, read collaborator comments
such as `\janote`, `\mgnote`, and any other comment command except
`\codexnote`. Analyze what each comment requires and implement the
appropriate solution. Preserve each original comment command and its
content so the user can review the changes and remove the comments when
satisfied. Flag any comment that cannot be resolved without further
information or justification.

<!-- *}}} -->
## Scientific Style and Writing <!-- {{{* -->

Use the `writing-science` skill, when available, both to review or revise
existing scientific prose and to draft new prose from the user's ideas.
The requested unit may range from a sentence or paragraph to a complete
section, but not a complete paper. Read the relevant passage in `main.tex`
and its surrounding context before revising; do not ask the user to paste
text already in the repository. Improve structure before wording, and check
logical flow, precision, and readability. For diagnosis-only requests,
report findings without editing.

**Non-negotiable citation requirement:** Whenever asked to write or revise
text in the main manuscript (`main.tex`), explicitly include the relevant
LaTeX citation next to each claim drawn from an external source. Read and
verify that the cited source supports the precise claim; do not rely on
a citation elsewhere in the document, a general reference, or attribution
only in the chat response. Reuse existing bibliography entries when
available and add verified entries to `references.bib` when needed.
For results established within this manuscript, explicitly reference the
relevant section, equation, or theorem instead of inventing an external
source. If a source or justification cannot be verified, flag the gap
rather than present the claim as supported. Complete this citation check
before delivering any requested manuscript prose.

Preserve scientific meaning, evidence, uncertainty, terminology, citations,
equations, labels, and LaTeX commands unless the user explicitly requests a
substantive change. Retain the manuscript's definitions, including the
unique free state and the fixed-point condition for free channels. Flag
missing justification rather than inventing proofs or citations, and do not
introduce unsupported priority or novelty claims. Follow journal
requirements supplied by the user; do not invent new ones.

<!-- *}}} -->
## Figures & Supporting Code <!-- {{{* -->

Store manuscript figures and their supporting scripts
under `Figures/`.

Use relative paths when including figures in LaTeX.

Do not introduce absolute Windows or Linux filesystem
paths into the manuscript.

When modifying a figure:
- Preserve the scientific data and definitions.
- Keep the source script when available.
- Check that the resulting figure agrees with the text.
- Recompile the manuscript to verify integration.

Do not introduce new software dependencies without
documenting them.

<!-- *}}} -->
## Version Control <!-- {{{* -->

Before editing, inspect the Git status, preserve unrelated changes, and
limit edits to the requested scope. Do not overwrite collaborators' work
or commit, push, pull, reset, restore, discard, or force-push changes
unless explicitly instructed. When asked to commit, use a short,
descriptive commit message.
<!-- *}}} -->
