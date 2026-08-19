# License Decision Required

The repository currently has **no license file**. Public GitHub code is **not** automatically
open-source licensed; without a license, others have no permission to reuse the code.

The project **author must choose a license**. Options to consider (author's decision, not made here):

- **MIT** — permissive; allows reuse with attribution.
- **Apache-2.0** — permissive; includes patent grant.
- **CC-BY-4.0** — for data/reports (not software).

Considerations specific to this project:

- The analysis consumes **public third-party data** (ClinVar, gnomAD, Ensembl, MaveDB, PubMed),
  each with its own terms; the project license applies to **our code and derived outputs**, not to
  the raw upstream data.
- The derived variant annotations are factual data (not creative expression), but the **code** and
  **reports** are copyrightable.

**Action required:** the author should add a `LICENSE` file (and update `CITATION.cff` `license`)
before any redistribution.
