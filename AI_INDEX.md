# AI index

## Identity
Start with [README](README.md), [STATUS](STATUS.md), [CLAIMS.json](CLAIMS.json), [CITATION.cff](CITATION.cff). Two-paper unrefereed collection; DOI is availability, not proof.

## Claim-to-proof map
A-fixed: [certificate paper](papers/certificates/paper.tex), thm:main/thm:nine, eta multiplier -> holomorphy -> U_M -> complete Sturm closure. B-mod8: [family paper](papers/families/paper.tex), thm:one, theta expansion plus general exclusion witnesses. B-support: thm:support, binary-form support and Frobenius parity. B-dyadic: thm:oddpart, local parity/valuation exclusion. [Combined PDF](paper.pdf) contains both papers.

## Verification
[verify.py](verify.py) runs the FLINT certificate and [standard-library reconstruction](checks/independent_bsc_audit.py), then [family checks](checks/families.py). [negative_controls.py](negative_controls.py) rejects wrong residues/moduli/levels/missing indices. [requirements](requirements.txt), [environment](ENVIRONMENT.txt), [replay scope](REPLAY_RECEIPT.md) give commands and boundaries. [check_package.py](check_package.py) verifies MANIFEST.json.

## Trust and limitations
[ASSURANCE](ASSURANCE.md) separates universal proofs and finite corroboration; no formal verification or unaffiliated reproduction. [Citation audit](CITATION_AUDIT.md) records known cases and priority uncertainty. [Revision response](reviews/revision-response.md) resolves objections. Mod16 dyadic parameter necessity and stronger subprogressions remain open here; do not infer them from samples. ell17,N21 is a required counterexample to a tempting stronger support claim.

## Reuse and provenance
[PROVENANCE](PROVENANCE.md), [licence map](LICENSES.md), [content licence](LICENSE), [code licence](LICENSE-CODE). No private reviews or third-party PDFs redistributed. Original content CC0; original code MIT; dependency licences retained. Supplied reviewer arithmetic was adapted with source provenance, not described as unaffiliated reproduction of the final bytes.
