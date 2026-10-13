# Security policy

## Reporting a vulnerability

Please do not report security problems in a public issue, pull request or
discussion.

Report them privately through GitHub instead:

1. Go to the [Security tab](https://github.com/M4i-Imaging-Mass-Spectrometry/thyra/security)
   of this repository.
2. Click **Report a vulnerability**.
3. Describe the problem, the Thyra version, and the steps or file that trigger it.
   Do not attach unpublished or patient data; describe the file instead.

Only the maintainers can see the report. We will reply as soon as we can, and we
will tell you when a fix is released.

## Supported versions

Fixes go into the next release of the latest version. Older versions do not get
security fixes, so please upgrade:

```bash
pip install --upgrade thyra
```

## Scope

Thyra reads files written by instrument vendors and by other tools. A file that makes
Thyra run code, write outside the output folder, or read files it was not given is in
scope. A file that only makes a conversion fail or run slowly is a normal bug; please
open an issue for it.
