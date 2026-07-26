# PROJECT.md — Authoritative Operating Instructions for Computational Biology Repository Builds

This file captures every rule, preference, standard, and workflow decision applied when transforming a generic scaffold into a professional, research-grade, reproducible computational docking repository (demonstrated here with **Molecular Docking of Indinavir Against HIV-1 Protease Using AutoDock Vina**). It is designed to be read and reused in **any future session** that follows this project template, so that the same quality, structure, and standards are applied consistently without re-explanation.

---

## Table of Contents

1. [Project Identity](#project-identity)
2. [Audience](#audience)
3. [Writing Style & Language](#writing-style--language)
4. [Level of Professionalism](#level-of-professionalism)
5. [Absolute Rules](#absolute-rules)
6. [Repository Structure Template](#repository-structure-template)
7. [Data Provenance Policy](#data-provenance-policy)
8. [Generated vs Committed Artifacts](#generated-vs-committed-artifacts)
9. [Images & Figures Policy](#images--figures-policy)
10. [Documentation Standards](#documentation-standards)
11. [Mermaid Diagrams](#mermaid-diagrams)
12. [Reference & Citation Standards](#reference--citation-standards)
13. [Git & Commit Policy](#git--commit-policy)
14. [Workflow for Future Docking Projects](#workflow-for-future-docking-projects)
15. [Quality Checklist](#quality-checklist)
16. [Lessons Learned from This Build](#lessons-learned-from-this-build)

---

## Project Identity

A computational biology portfolio project documenting a complete molecular docking experiment. The repository is the primary artifact; it must read like a published, reproducible scientific study — not a tutorial, not a classroom assignment, not a blog post.

## Audience

Read with an undergraduate-level molecular biology background and computational fluency. Do **not** oversimplify or include unnecessary introductory biology. Trust the reader to understand concepts like protease dimers, active sites, and hydrogen bonding. Explain domain-specific *computational* decisions (e.g., why Kollman charges, why grid box size matters) but not basic biology.

## Writing Style & Language

- **Technical** — precise, domain-correct terminology (e.g., "homodimer," "Gasteiger charges," "torsion tree," "RMSD lower/upper bound").
- **Concise** — no filler, no excessive warmth, no emoji in docs.
- **Scientifically accurate** — every claim must be verifiable; never fabricate numbers, scores, or filenames. If a value is unknown, state so explicitly (e.g., placeholder).
- **Professionally formatted** — consistent heading hierarchy, tables where appropriate, DOIs on all citations, code blocks for commands and config files, callout blocks (`> **Figure X.** ...`) for screenshots.

## Level of Professionalism

The repository should be indistinguishable in quality from a portfolio piece submitted by an experienced computational biologist. Every file has a purpose. Nothing is template filler. Placeholders are honest ("Insert Screenshot Here") rather than silent gaps. Real numbers from real runs replace placeholder tables. The README front-loads the scientific objective so a recruiter or admissions committee member can assess the project in 30 seconds.

---

## Absolute Rules

These rules are **non-negotiable** in every session that builds or updates a repository following this template.

### 1. Never Commit on Behalf of the User

I will never run `git commit` or `git push`. The user commits and pushes themselves. When work is complete, I leave the working tree clean with all changes unstaged and tell the user it is ready. This rule persists across sessions (see the `no-auto-commit` memory entry).

### 2. Never Fabricate Numerical Values

No docking scores, RMSD values, binding affinities, grid coordinates, or any quantitative result shall be invented. If a real run exists, use the real values from the actual output file. If no run has been performed, the placeholder is left as text (e.g., "Docking score placeholder — do not cite until a real Vina run is performed") or omitted entirely.

### 3. Never Commit Generated Binary or Placeholder Chemistry Data

`.pdbqt` prepared files, `.map`/`.fld`/`.gpf` grid files, and binary images not captured by the user are not committed. Only real source inputs (raw PDB from RCSB, SDF from PubChem) and real output files (if a run was performed and the output exists) are committed. The `.gitignore` enforces this.

### 4. Use Authoritative, Verifiable Inputs

Protein structures come from the Protein Data Bank (RCSB PDB). Ligands come from PubChem (using the correct CID, with formula/formula verification). Grid box coordinates are derived from the actual native ligand coordinates in the PDB file, not guessed.

### 5. Cite DOIs on All References

Every referenced paper, software tool, or database has a DOI link in `docs/references.md` and `CITATION.cff`.

### 6. Be Honest About Limitations

Every docking study has limitations (rigid receptor, empirical scoring, no explicit solvent, single conformer). These are documented in `docs/discussion.md` with specific, technical explanations — not hand-waving.

### 7. Distinguish Between What Docking Can and Cannot Prove

Docking predicts a ranked binding pose and an *empirical* affinity score. It does not prove activity, calculate a thermodynamic free energy, or model water-mediated interactions at equilibrium. The discussion makes this distinction explicit.

### 8. Preserve User Identity

The repository author is the user. No AI assistant name or handle is credited as author or contributor. The CITATION.cff and README `## Author` section use the user's real name and email. The git user identity (`clawncore`) matches the GitHub account hosting the repository, and the repo URL in CITATION.cff / README points to the correct repository location.

### 9. Respect the Existing README Structure

The README follows the exact section order specified in the project brief. It is the front door of the portfolio — professional, scannable, complete.

### 10. Images Are Optional; Don't Fabricate Them

Only real screenshots captured by the user are included. If an image does not exist, it is not listed as a "placeholder" in the repository structure; the reference is either omitted or noted as "(optional)" in the figure caption with an honest instruction to add it if captured.

---

## Repository Structure Template

```
.
├── README.md                     # Front door: objective, workflow, results summary, limitations
├── LICENSE                       # CC BY 4.0 (or as specified by user)
├── .gitignore                    # Exclude .env, OS junk, AutoDock maps, raw export folders
├── CHANGELOG.md                  # Versioned change log
├── CITATION.cff                  # Machine-readable citation (CFF v1.2)
├── PROJECT.md                    # ← This file (operating instructions for future sessions)
├── docs/
│   ├── introduction.md           # Biological & methodological background
│   ├── methodology.md            # Step-by-step protocol, every parameter justified
│   ├── results.md                # Real output data, score table, RMSD interpretation
│   ├── discussion.md             # Scientific interpretation, limitations, validation needs
│   ├── conclusion.md             # Summary, lessons learned, future directions
│   ├── troubleshooting.md        # Comprehensive error/fix guide
│   ├── references.md             # Full citations with DOIs
│   └── diagrams.md               # Mermaid diagrams (workflow, prep, docking, repo, viz)
├── protein/
│   ├── protein.pdb               # Raw or prepared receptor (RCSB source)
│   └── protein.pdbqt             # Prepared receptor (AutoDockTools)
├── ligand/
│   ├── indinavir.sdf             # Ligand 3D conformer (PubChem)
│   ├── ligand.pdb                # Ligand in PDB format
│   └── ligand.pdbqt              # Prepared ligand (AutoDockTools)
├── config/
│   └── config.txt                # Vina: receptor, ligand, grid center/size, run params
├── output/
│   ├── out.pdbqt                 # Vina ranked poses (real, if run performed)
│   └── log.txt                   # Vina stdout (real, if run performed)
└── images/                       # Only real screenshots; optional figures noted as (optional)
    ├── protein.png
    └── ligand.png
```

### Notes on Structure

- `docs/weekly/`, `notebooks/`, `scripts/`, `data/`, `results/`, `configs/` — if these existed as generic scaffold, remove them. They clutter the portfolio.
- Keep `LICENSE`, `CHANGELOG.md`, `CITATION.cff`, and `PROJECT.md` at root.
- The doc directory (`docs/`) mirrors the structure of a real research paper: background → methods → results → discussion → conclusion.

---

## Data Provenance Policy

| Artifact Type | Source | Verification |
|---------------|--------|-------------|
| Protein structure | RCSB PDB | Download via `https://files.rcsb.org/download/<PDB_ID>.pdb`; verify header line |
| Ligand 3D conformer | PubChem | Download via REST API or web; verify molecular formula, IUPAC name, heavy-atom count match expected compound |
| CID disambiguation | PubChem REST | `curl -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/<CID>/property/IUPACName,MolecularFormula/JSON"` |
| Native ligand coordinates | PDB file itself | Extract `HETATM` records of co-crystallized ligand; compute geometric centroid for grid box center |

**Important:** If the user provides a folder of real data (PDB, SDF, PNG screenshots, docking output), use those files directly. Do not re-download or overwrite user's real artifacts with generic placeholder downloads. Place them according to the structure above.

---

## Generated vs Committed Artifacts

| Artifact | Committed? | Rationale |
|----------|-----------|-----------|
| Raw PDB from RCSB | ✅ Yes | Real input; reproducible download |
| Ligand SDF from PubChem | ✅ Yes | Real input; verifiable |
| `protein.pdbqt` (prepared) | ✅ Yes (if real run done) | Otherwise exclude via `.gitignore` |
| `ligand.pdbqt` (prepared) | ✅ Yes (if real run done) | Otherwise exclude |
| `output/out.pdbqt` | ✅ Yes (if real run done) | The actual docking result — valuable |
| `output/log.txt` | ✅ Yes (if real run done) | Run metadata + scores |
| `.map`, `.fld`, `.gpf` files | ❌ No | Regenerable; ` .gitignore` excludes |
| `Docking HIV 1/` (source folder) | ❌ No | ` .gitignore` excludes |
| `Docking HIV 1.zip` | ❌ No | ` .gitignore` excludes |

The `.gitignore` must be kept consistent with this policy.

---

## Images & Figures Policy

- Only the user's real screenshots are committed (`images/protein.png`, `images/ligand.png`, etc.).
- Placeholder image references in the repository structure block are **removed** if no file exists for that entry.
- Figure captions use `> **Figure X.**` callout format. Optional figures are marked `(optional)`.
- The README's `## Author` section is the last section and is factual (name + email + portfolio note).
- No AI assistant is credited as author, contributor, or generator.

---

## Documentation Standards

### `docs/methodology.md` (most important doc)

- Documents every step the user actually performed.
- Numbers are real (from real config.txt and real log.txt).
- Command lines are given for Windows (Command Prompt / PowerShell).
- Every parameter is explained (receptor path, center_x/y/z, size_x/y/z, exhaustiveness, num_modes, energy_range).
- The grid box derivation is explained (centroid of native ligand heavy atoms).
- Output files are explained (`.pdbqt` MODEL blocks, `log.txt` table).
- A reproducibility checklist is included at the end.

### `docs/results.md`

- Real score table from `output/log.txt`.
- Column definitions (affinity, RMSD l.b./u.b.).
- Headline findings (best affinity, energy spread, clustering).
- Interpretation caveats (rankings, not absolute free energies).

### `docs/discussion.md`

- Written in scientific paper style.
- Predicteds binding → biological relevance → strengths → limitations → validation importance.
- A table of potential error sources (grid placement, flexibility, scoring function, protonation state, conformer sampling, seed sensitivity).
- Explicit "can/cannot conclude" section.

### `docs/introduction.md`

- Assumes undergraduate molecular biology.
- Covers: HIV → HIV-1 protease → why it's a target → mechanism of Indinavir → role of docking → limitations of docking.
- No unnecessary biology explanations.

### `docs/troubleshooting.md`

- Organized by tool area (ADT, PDBQT export, grid, charges, Vina executable, PyMOL, file conversion, permissions, PATH).
- Each entry: symptom → cause → fix.

### `docs/references.md`

- All entries with DOIs.
- Includes software citations (Vina, ADT, PyMOL, Open Babel, PDB, PubChem) and biological literature.

### `docs/diagrams.md`

- Mermaid diagrams for: overall workflow, receptor prep, ligand prep, docking workflow, repository structure, visualization pipeline.

---

## Mermaid Diagrams

Mermaid is used for flowcharts embedded in doc files. No external image files needed for diagrams. Each diagram is a ` ```mermaid ` fenced code block. Key diagrams to include:

1. **Overall workflow** — PDB → receptor prep → ligand prep → Vina → output → visualization.
2. **Receptor preparation** — PDB → delete waters → delete ligand → add polar H → Kollman charges → PDBQT.
3. **Ligand preparation** — SDF → add H → Gasteiger charges → detect torsions → PDBQT.
4. **Docking workflow** — PDBQT inputs + config → Vina → MODEL blocks + REMARK scores.
5. **Repository structure** — tree of files.
6. **Visualization pipeline** — PDBQT → PyMOL → cartoon/sticks/color/zoom → PNG export.

---

## Reference & Citation Standards

- **CFF file** (`CITATION.cff`): schema v1.2, type `software`, version, date-released, license (CC-BY-4.0), URL pointing to the correct repository, keywords, and references to Vina paper with DOI.
- **README BibTeX**: `@software{...}` with `author = {Family, Given}`, `title`, `year`, `note`, `url`.
- **DOI links** on all references in `docs/references.md`.

---

## Git & Commit Policy

- **I never commit or push.** The user controls all git operations.
- All my work is staged but uncommitted (if any) — the user reviews and commits themselves.
- Stale scaffold files (if the user wants them removed) can be deleted at the user's explicit request. I flag what I removed but do not delete without asking unless the files are clearly template stubs that contradict the requested clean structure.
- `.gitignore` is maintained to exclude regenerable/OS/export artifacts.

---

## Workflow for Future Docking Projects (Reusable Template)

When the user asks to build a similar docking repository for a **new** protein–ligand system:

1. **Acquire real inputs** — PDB from RCSB, ligand SDF from PubChem (verify CID and molecular formula).
2. **Check for user-provided real data** — if a folder exists with prepared files, screenshots, and docking output, use those directly.
3. **Derive grid box** — extract native ligand coordinates, compute centroid (`awk` on HETATM records), choose box size (ligand span + 3–5 Å padding each direction).
4. **Write `config/config.txt`** with real center, size, and run parameters.
5. **Write `README.md`** using this template structure.
6. **Write all `docs/*.md`** — introduction, methodology, results, discussion, conclusion, troubleshooting, references, diagrams.
7. **Use real output values** — never fabricate scores, never use "placeholder" tables when real data exists.
8. **Commit author as the user** — name, email, repo URL.
9. **Place real PNGs** in `images/`; omit non-existent files.
10. **Mark the session complete** — tell the user it is ready for commit.

---

## Quality Checklist

Before marking a build session complete, verify:

- [ ] `README.md` has all required sections (overview, objective, background, workflow diagram, structure, methodology summary, software, install, execution, results summary, limitations, future work, references, citation, license, author).
- [ ] All numerical values in docs come from real outputs (config.txt, log.txt, out.pdbqt).
- [ ] `config/config.txt` center matches the native ligand centroid.
- [ ] `output/log.txt` or `output/out.pdbqt` contains the real Vina output.
- [ ] No fabricated docking scores, RMSD values, or molecule names.
- [ ] No AI assistant credited as author or contributor.
- [ ] CITATION.cff has correct repo URL, author name/email, CC BY 4.0 license.
- [ ] `docs/references.md` includes DOIs for all entries.
- [ ] `.gitignore` excludes generated maps, OS junk, and source export folder.
- [ ] Only real screenshots exist in `images/` — no placeholder image references unless real files exist.
- [ ] Mermaid diagrams render correctly (valid Mermaid syntax).
- [ ] Discussion explicitly distinguishes what docking can vs. cannot conclude.
- [ ] Limitations section is specific and technical, not vague.
- [ ] The user has **not** been asked to commit; all changes are in the working tree.

---

## Lessons Learned from This Build

1. **Always use the user's real data first.** They provided a `Docking HIV 1/` folder with prepared PDBQT, real docking output, and PNG screenshots — I used those directly rather than re-downloading generic downloads. This made the repository self-contained and authentic.

2. **Verify ligand identity carefully.** The SDF was labeled PubChem CID 5362440, not the canonical Indinavir CID 152648. REST API confirmed C₃₆H₄₇N₅O₄ and matching IUPAC name — same molecule (Indinavir free base). Document both CIDs transparently.

3. **Real docking output is gold.** The user's `ligand_out.pdbqt` contained actual `REMARK VINA RESULT` lines with real affinity and RMSD values. Populating `results.md` and `output/log.txt` with real data — not placeholders — made the portfolio credible.

4. **`receptor.pdbqt` ≠ `protein.pdbqt`.** The folder contained two receptor pdbqt files. The `config.txt` references `protein.pdbqt` (the Vina-prepared one), so that's the authoritative prepared receptor to commit. `receptor.pdbqt` is an AutoDock leftover (generated from the `.gpf` grid parameter file) — not needed.

5. **Grid box center varies slightly.** My independent computation from MK1 gave `(13.073, 22.467, 5.557)`; the actual run used `(13.043, 22.460, 5.596)`. Both are centroids of MK1 atoms; minor differences arise from which atoms are included or whether hydrogens are counted. Use the actual run values in the config.

6. **`ligand_out.pdbqt` → `output/out.pdbqt`.** Rename for cleanliness. Update README and methodology to reference `output/out.pdbqt`.

7. **RMSD interpretation matters.** The user's run did **not** pass `--ref_ligand`, so RMSD is relative to the best pose, not the crystallographic ligand. Correct the README/methodology to avoid implying reference-ligand RMSD. Document `vina --ref_ligand` as an enhancement for users who want RMSD-to-crystal.

8. **Don't over-engineage images.** The user explicitly said "no need for all those images bro." Two real screenshots are worth more than eight fabricated placeholders. Only commit real images; mention others as optional.

9. **Clean the scaffold ruthlessly.** Generic template directories (`configs/`, `data/`, `scripts/`, `notebooks/`, `results/`, `docs/weekly/`) and files (`PROGRESS.md`, `CONTRIBUTORS.md`, `requirements.txt`, `environment.yml`) are not part of a focused portfolio repo. Remove them; git retains them for recovery.

10. **PROJECT.md is the persistence layer.** This file ensures every future session inherits all rules, preferences, and workflow decisions without re-explanation. It is the single most important file for session continuity.
