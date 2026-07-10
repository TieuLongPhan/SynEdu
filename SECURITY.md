# Security Policy

## Supported Versions

Security fixes are applied to the latest release and the `main` branch.

## Reporting A Vulnerability

Please use GitHub's private vulnerability reporting for this repository:

1. Open the repository's **Security** tab.
2. Select **Advisories** and **Report a vulnerability**.
3. Include affected versions, reproduction details, and expected impact.

Do not open a public issue for an undisclosed vulnerability. If private
reporting is unavailable, contact the maintainer listed in `CITATION.cff` and
avoid including exploit details in the first message.

You should receive an acknowledgement within seven days. Confirmed issues will
be coordinated through a GitHub security advisory before public disclosure.

## Build Provenance

Tagged builds produce wheel and source distributions with SHA-256 checksums and
GitHub artifact attestations. After downloading a release artifact, verify it
with:

```bash
gh attestation verify <artifact> --repo TieuLongPhan/SynEdu
```
